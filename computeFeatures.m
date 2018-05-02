function features = computeFeatures(folder,song)
%%TO DO  Add features / save all features

%% Initial variables
startbpm = 120;
dir = [folder song];
[waveIn,sr] = audioread(dir);
features.name = song;

%% Convert stereo to mono
[m,n] = size(waveIn);
if n == 2
   waveIn = sum(waveIn, 2) / size(waveIn, 2); 
end
features.wave = waveIn;  
features.sr = sr;
%% Beat times 
features.beatPos = beat2(waveIn,sr); %beat2
beatDuration=features.beatPos(2:end)-features.beatPos(1:end-1);
bpm=1/mean(beatDuration)*60;
ratio=max(bpm,startbpm)/min(bpm,startbpm);
if(ratio>1.5) %set bpm to half
    features.beatPos=features.beatPos(1:2:end);
    beatDuration=features.beatPos(2:end)-features.beatPos(1:end-1);
end
features.beatDuration=[beatDuration beatDuration(end)];
features.bpm = bpm;
%% Chroma features
fftlen = 2 ^ (round(log(sr*(2048/22050))/log(2)));
nbin = 12;
f_ctr = 1000; 
f_sd = 1; 
chromaF = chromagram_IF(waveIn,sr,fftlen,nbin,f_ctr,f_sd);
ffthop = fftlen/4;
sgsrate = sr/ffthop;
features.bsChromaF = beatavg(chromaF,features.beatPos*sgsrate);
imagesc(features.bsChromaF)
plot(features.bsChromaF(1,:))

%% MFCC and energy
%[cepstra,aspectrum,pspectrum] = melfcc(waveIn, sr, 'wintime', fftlen/sr, 'hoptime', ffthop/sr, 'numcep', 12); %careful here, there are two melfcc functions in my path. Use that of labrosa/rastamaa
%imagesc(cepstra);
%features.bsMFCC = beatavg(cepstra,features.beatPos*sgsrate);

%% compute RMS energy
features.energyTime=energyT(waveIn, fftlen, ffthop);
features.energyStart=features.energyTime(max(1,round(features.beatPos*sr/ffthop)));

%% Energy power
eng = features.energyTime;
eng = eng - mean(eng);
sv = eng.* eng;
dp = sum(sv) / length(sv);     % suma de quadrados, producto escalar
features.power = sqrt(dp);   % sqrt del producto escalar

%% metrical position
features.metricalPos=ones(length(features.energyStart),1);
for i=1:4
    features.metricalPos(i:4:end)=i;
end

end
