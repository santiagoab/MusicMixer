function [closest]= compBeatSimilarity(descriptors, transitionLimits)
doSave=0;
% %% ------------- compute chroma features
% %t = tempo(waveIn,sr);
% %bsChromaF = chrombeatftrs(waveIn,sr, beatPos);
% fftlen = 2 ^ (round(log(sr*(2048/22050))/log(2)));
% nbin = 12;
% f_ctr = 1000; 
% f_sd = 1; 
% chromaF = chromagram_IF(waveIn,sr,fftlen,nbin,f_ctr,f_sd);
% ffthop = fftlen/4;
% sgsrate = sr/ffthop;
% bsChromaF = beatavg(chromaF,beatPos*sgsrate);
% %imagesc(bsChromaF)
% %plot(bsChromaF(1,:))

% %% -------------- compute MFCC and energy
% [cepstra,aspectrum,pspectrum] = melfcc(waveIn, sr, 'wintime', fftlen/sr, 'hoptime', ffthop/sr, 'numcep', 12);
% %imagesc(cepstra);
% bsMFCC = beatavg(cepstra,beatPos*sgsrate);

% %% ---compute RMS energy
% %rmsEnergy=computeRMSEnergy(waveIn, sr, fftlen, ffthop);
% energyTime=energyT(waveIn, fftlen, ffthop);
% energyStart=energyTime(round(beatPos*sr/ffthop));
% %energyMax=energyTime(round(beatPos*sr/ffthop));
% %[energyTime, energySpectrum]=energy(waveIn,sr);

%% duration, energy, metrical position
energyStartM=ones(length(descriptors.energyStart), length(descriptors.energyStart));
durationM=ones(length(descriptors.energyStart),length(descriptors.energyStart));

for i=1:length(descriptors.energyStart)
    for j=1:length(descriptors.energyStart)
        energyStartM(i,j)=abs(descriptors.energyStart(i)- descriptors.energyStart(j));
        durationM(i,j)=abs(descriptors.beatDuration(i)- descriptors.beatDuration(j));
    end
end

%% metrical position
metricalM=ones(length(descriptors.energyStart),length(descriptors.energyStart));
%metrical=ones(length(descriptors.energyStart),1);
for i=1:4
 %   metrical(i:4:end)=i;
    metricalM(i:4:end, i:4:end)=0;
end
%% ----- weight features
timbreWeight = 1;
pitchWeight = 1;
loudStartWeight = 1;
loudMaxWeight = 1; 
durationWeight = 1;
%confidenceWeight = 1;
metricalPosWeight = 1;

%normalize matrices
bsMFCCNorm= (descriptors.bsMFCC-min(min(descriptors.bsMFCC)))./max(max(descriptors.bsMFCC));
bsChromaFNorm= (descriptors.bsChromaF-min(min(descriptors.bsChromaF)))./max(max(descriptors.bsChromaF));
energyStartMNorm= (energyStartM-min(min(energyStartM)))./max(max(energyStartM));
durationMNorm= (durationM-min(min(durationM)))./max(max(durationM));
metricalMNorm= (metricalM-min(min(metricalM)))./max(max(metricalM));
% subplot(2,2,1)
% imagesc(dist(descriptors.bsMFCC)); title('MFCC')
% subplot(2,2,2); 
% imagesc(dist(bsChromaF));title('Chroma')
% subplot(2,2,3)
% imagesc(energyStartM); title('Energy Onset')
% subplot(2,2,4)
% imagesc(durationM); title('Duration')

%wF=(timbreWeight*dist(bsMFCCNorm)+pitchWeight*dist(bsChromaFNorm)+loudStartWeight*energyStartMNorm+durationWeight*durationMNorm ).*metricalM; %+ metricalPosWeight*metricalMNorm;
wF=timbreWeight*dist(bsMFCCNorm)+pitchWeight*dist(bsChromaFNorm)+loudStartWeight*energyStartMNorm+durationWeight*durationMNorm+ metricalPosWeight*metricalMNorm;
wFNorm= (wF-min(min(wF)))./max(max(wF));

%imagesc(wF);
%wF=bsMFCC;
%
autosimilarityMatrix=wFNorm;
%dist(wF);
%autosimilarityMatrix=corrmtx(bsChromaF',bsChromaF');

% imagesc(autosimilarityMatrix)
% imagesc(autosimilarityMatrix(1:10,:))
% plot(log(autosimilarityMatrix(1,:)))
% plot(autosimilarityMatrix(1,:))


%% find closest beats. Threshold is ?
threshold=0.06;
maxnCandidates=5;
minsep=15;
clear closest;
n=length(descriptors.beatPos);
closetMat=zeros(n,n);
for i=transitionLimits(1)+1:size(autosimilarityMatrix,2)-transitionLimits(2)
    %get the n-closest below threshold and separated min of minsep
    [ASorted AIdx] = sort(autosimilarityMatrix(i,:));
    smallestNElements = ASorted(1:maxnCandidates);
    smallestNIdx = AIdx(1:maxnCandidates);
    idxInSmallest=find(smallestNElements<threshold);
    idx=sort(smallestNIdx(idxInSmallest));
    %remove neighbor same and next. 
    %idx=setxor(idx,i);
    for k=-minsep:minsep
       if(find(ismember(idx,i+k)))
           idx=setxor(idx,i+k);
       end
    end 
    %remove connections to transitions
    idxtr1=find(idx>transitionLimits(1));
    idx=idx(idxtr1);
    idxtr2=find(idx<size(autosimilarityMatrix,2)-transitionLimits(2)+1);
    idx=idx(idxtr2);
    %As js is 0-indexed, then (i-1) and (i)
    idx=idx-1;
    closest{i}=idx; 
    %closetMat(i,idx)=ones(length(idx),1);
end

if(doSave)
    % save JSON
    jsonString=savejson('',closest);
    fid=fopen([dir fileName '_closest.json'],'w');
    fprintf(fid, '%s', jsonString);
    fclose(fid);

    fid=fopen(['C:\wamp\www\tiles\' fileName '_closest.json'],'w');
    fprintf(fid, '%s', jsonString);
    fclose(fid);
end
