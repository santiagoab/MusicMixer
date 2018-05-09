function synthSignal = pl_phaseVocoder_variable_analysis_hop(signal, tsm_factor, winSamps)
% A phase locked vocoder time-scale modification algorithm based upon Jean
% Laroche's 1999 work. The synthesis hop is fixed at one quarter the analysis window (hanning) 
% while the analysis hop is scaled by the time scale factor; this results in more FFT computations 
% than the bonada 2000 approach for slowing down - but provides smoother transitions frame to frame of 
% the magnitude spectrums and, in general, a better quality of output. This code started
% out as the code provided by Tae Hong Park, but has changed significantly
% over the years
%
% David Dorran, Audio Research Group, Dublin Institute of Technology
% david.dorran@dit.ie
% http://eleceng.dit.ie/dorran
% http://eleceng.dit.ie/arg
%
 
% make sure input is mono and transpose if necessary
[r, c] = size(signal);
if r > 1
    signal = signal';
end;
[r, c] = size(signal);
if r > 1
    disp('Note :only works on mono signals');
    synthSignal = [];
    return
end;
 
% add a few zeros to stop the algorithm failing
zpad = zeros(1, 44100/4);
signal = [signal, zpad];
if nargin < 3
    winSamps = 2048;
end
 
winSampsPow2 = winSamps;
synHopSamps = winSampsPow2/4;
anHopSamps = round(synHopSamps/tsm_factor);
 
win = hanning(winSampsPow2);
 
X = specgram(signal, winSampsPow2, 100,win, winSampsPow2 - anHopSamps);
 
moduli = abs(X);
phases = angle(X);
 
[numBins, numFrames ] = size(phases);
 
syn_phases = zeros(numBins, numFrames); % a holder for synthesis phases
 
twoPi   = 2*pi;
omega   = twoPi * anHopSamps * [0:numBins-1]'/numBins/2; %the expected phase hop per frame
 
syn_phases(:,1) = phases(:,1) .* ( synHopSamps/ anHopSamps);
 
for idx =  2: numFrames
    ddx = idx - 1;
    deltaPhi = princarg(phases(:,idx) - phases(:,ddx) -omega); %calculate priciple determination of the hetrodyned phase increment
    phaseInc = (omega+deltaPhi)/anHopSamps; % phase increment per sample
    %locate the peaks
    pk_indices = [];
    pk_indices  =  locate_peaks(moduli(:,idx));
    if(~length(pk_indices))
        pk_indices = [1 10 12]; % just in case an odd situation is encountered  e.g. a sequence of zeros
    end
    %update phase of each peak channel using the phase propagation formula
    syn_phases(pk_indices,idx)    = syn_phases(pk_indices,ddx)+synHopSamps*phaseInc(pk_indices); %synthesis phase
    %update phase of channels in region of influence
    % first calculate angle of rotation
    rotation_angles = syn_phases(pk_indices,idx) - phases(pk_indices,idx);
    start_point = 1; %initialize the starting point of the region of influence
 
    for k = 1: length(pk_indices) -1
        peak = pk_indices(k);
        angle_rotation  = rotation_angles(k);
        next_peak = pk_indices(k+1);
        end_point = round((peak + next_peak)/2);
        ri_indices = [start_point : peak-1, (peak+1) : end_point]; %indices of the region of influence
        syn_phases(ri_indices,idx) = angle_rotation + phases(ri_indices, idx);
        start_point = end_point + 1;
    end;
end;
 
%Make sure that the LHS and RHS of the DFT's of the synthesis frames are a
%complex conjuget of each other
Z = moduli.*exp(i*syn_phases);
Z = Z(1:(numBins),:);
conj_Z = conj(flipud(Z(2:size(Z,1) -1,:)));
Z = [Z;conj_Z];
 
synthSignal = zeros(round(length(signal)*tsm_factor+length(win)), 1);
 
curStart = 1;
for idx = 1:numFrames-1
    curEnd   = curStart + length(win) - 1;
    rIFFT    = real(ifft(Z(:,idx)));
    synthSignal([curStart:curEnd]) = synthSignal([curStart:curEnd]) + rIFFT.*win;
    curStart = curStart + synHopSamps;
end
 
%--------------------------------------------------------------------------
function indices = locate_peaks(ip)
%function to find peaks
%a peak is any sample which is greater than its two nearest neighbours
    index = 1;
    num = 2;
    indices = [];
    for k = 1 + num : length(ip) - num
        seg = ip(k-num:k+num);
        [val, max_index] = max(seg);
        if max_index == num + 1
            indices(index) = k;
            index = index + 1;
        end;
    end;
     
%--------------------------------------------------------------------------
function Phase = princarg(Phasein)
    two_pi = 2*pi;
    a = Phasein/two_pi;
    k = round(a);
    Phase = Phasein-k*two_pi;