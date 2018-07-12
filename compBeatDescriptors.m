function descriptors = compBeatDescriptors(waveIn, sr, startbpm) %dir, fileName)

if nargin < 3;   startbpm = 120; end

%% -------- Calculate the beat times
descriptors.beatPos = beat2(waveIn,sr); %beat2
beatDuration=descriptors.beatPos(2:end)-descriptors.beatPos(1:end-1);
bpm=1/mean(beatDuration)*60;
ratio=max(bpm,startbpm)/min(bpm,startbpm);
if(ratio>1.5) %set bpm to half
    descriptors.beatPos=descriptors.beatPos(1:2:end);
    beatDuration=descriptors.beatPos(2:end)-descriptors.beatPos(1:end-1);
%elseif(ratio<1.5) %set bpm to double?
end
descriptors.beatDuration=[beatDuration beatDuration(end)];

%% ------------- compute chroma features
%  NO VA t = tempo2(waveIn,sr);
% NO VA ERRORES descriptors.bsChromaF = chrombeatftrs(waveIn,sr, beatPos);
fftlen = 2 ^ (round(log(sr*(2048/22050))/log(2)));
nbin = 12;
f_ctr = 1000; 
f_sd = 1; 
chromaF = chromagram_IF(waveIn,sr,fftlen,nbin,f_ctr,f_sd);
ffthop = fftlen/4;
sgsrate = sr/ffthop;
descriptors.bsChromaF = beatavg(chromaF,descriptors.beatPos*sgsrate);
%imagesc(descriptors.bsChromaF)
%plot(descriptors.bsChromaF(1,:))

%% -------------- compute MFCC and energy
[cepstra,aspectrum,pspectrum] = melfcc(waveIn, sr, 'wintime', fftlen/sr, 'hoptime', ffthop/sr, 'numcep', 12); %careful here, there are two melfcc functions in my path. Use that of labrosa/rastamaa
%imagesc(cepstra);
descriptors.bsMFCC = beatavg(cepstra,descriptors.beatPos*sgsrate);

%% --- compute RMS energy
%rmsEnergy=computeRMSEnergy(waveIn, sr, fftlen, ffthop);
descriptors.energyTime=energyT(waveIn, fftlen, ffthop);
descriptors.energyStart=descriptors.energyTime(round(descriptors.beatPos*sr/ffthop));
%energyMax=descriptors.descriptors.energyTime(round(beatPos*sr/ffthop));
%[descriptors.descriptors.energyTime, energySpectrum]=energy(waveIn,sr);

%% metrical position
descriptors.metricalPos=ones(length(descriptors.energyStart),1);
for i=1:4
    descriptors.metricalPos(i:4:end)=i;
end

%% normalize features
%descriptors.bsMFCCNorm= (descriptors.bsMFCC-min(min(descriptors.bsMFCC)))./max(max(descriptors.bsMFCC));
%descriptors.bsChromaFNorm= (descriptors.bsChromaF-min(min(descriptors.bsChromaF)))./max(max(descriptors.bsChromaF));
%descriptors.energyStartMNorm= (descriptors.energyStartM-min(min(descriptors.energyStartM)))./max(max(descriptors.energyStartM));


