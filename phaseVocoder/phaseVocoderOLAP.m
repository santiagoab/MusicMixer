%---------phase vocoder with overlap-add (real time buffers)

%%
dir='D:\projects\FASTQMUL\silentDisco\sounds\';
fileName='03-FreakFandango-GypsySong';
[waveIn,sr] = wavread([dir fileName '.wav'],[1000000 2000000]);

% Initialize some variables
WindowLen = 256;
hws=WindowLen/2;
AnalysisHop = 64;
SynthesisHop = round(AnalysisHop*(1/0.9));
Hopratio = SynthesisHop/AnalysisHop;
window=[hanning(WindowLen);0]; %hanning(WindowLen,'periodic');
%dsp.Window('Hanning', 'Sampling', 'Periodic');%
zp=1;
fftSize=(WindowLen)*zp;

% buffer to window and apply the FFT
%grain=zeros(WindowLen,1);

% Initialize the variables used in the processing loop.
yprevwin = zeros(WindowLen-SynthesisHop, 1);
gain = sum(window.^2)/WindowLen; %1/(WindowLen*sum(window.^2)/SynthesisHop);
unwrapdata = 2*pi*AnalysisHop*(0:WindowLen-1)'/WindowLen;
yangle = zeros(WindowLen, 1);
firsttime = true;

%% 
b = 1;
e = length(waveIn);
b = max([hws+1,b]);    
e = min([length(waveIn)-hws,e]);    
n = length(b:AnalysisHop:e);
waveOut=[];
clear olapadd;

for c=b:AnalysisHop:e
    grainIn=waveIn(c-hws:c+hws-1) .* window(1:WindowLen);
    fftIn=fft(grainIn);
    
    % Convert complex FFT data to magnitude and phase.
    ymag       = abs(fftIn);
    yprevangle = yangle;
    yangle     = angle(fftIn);
    
    % Synthesis Phase Calculation
    % The synthesis phase is calculated by computing the phase increments
    % between successive frequency transforms, unwrapping them, and scaling
    % them by the ratio between the analysis and synthesis hop sizes.
    yunwrap = (yangle - yprevangle) - unwrapdata;
    yunwrap = yunwrap - round(yunwrap/(2*pi))*2*pi;
    yunwrap = (yunwrap + unwrapdata) * Hopratio;
    if firsttime
        ysangle = yangle;
        firsttime = false;
    else
        ysangle = ysangle + yunwrap;
    end    
    % Convert magnitude and phase to complex numbers.
    fftOut = ymag .* complex(cos(ysangle), sin(ysangle));

    wGrainOut=real(ifft(fftOut));

    % Overlap-add operation
    olapadd  = [wGrainOut(1:end-SynthesisHop,:) + yprevwin; ...
                wGrainOut(end-SynthesisHop+1:end,:)];
    yistfft  = olapadd(1:SynthesisHop,:);
    yprevwin = olapadd(SynthesisHop+1:end,:);

    % Compensate for the scaling that was introduced by the overlap-add
    % operation
    yistfft = yistfft * gain;
    %grainOut=yistfft./window(1:WindowLen);
    waveOut=[waveOut ; yistfft];

end

waveOut(1:100)=zeros(100,1);


% plot(waveOut)
%wavwrite(waveIn,sr,'D:\matlab\FAST\sounds\original.wav');
wavwrite(waveOut,sr,'D:\matlab\FAST\sounds\timestrechtedPVOLA.wav');

%%
if (0)
    subplot(4,1,1)
    plot(waveIn(c-hws:c+hws-1));
    hold on
    plot(grainIn,'g')
    title('grain and windowed grain')
    subplot(4,1,2)
    plot(ymag)
    hold on
    plot(abs(fftOut),'g')
    title('ymagIn and ymagOut')
    subplot(4,1,3)
    plot(yangle)
    hold on
    plot(angle(fftOut),'g')
    title('phaseIn and phaseOut')
    subplot(4,1,4)
    plot(wGrainOut)
    hold on
    %plot(ywin,'g')
    title('grainOut and Windowed grain out')
end