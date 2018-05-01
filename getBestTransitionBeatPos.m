function [trBeatIdx]=getBestTransitionBeatPos(descriptors)
    doPlots=0;
   %% few plots
    if (doPlots)
        plot(descriptors(1).energyStart)
        hold on
        plot(descriptors(2).energyStart,'g')
    end

    %% energyStart
    maxLengthTransitionBeats=min([length(descriptors(1).energyStart)/2, length(descriptors(2).energyStart)/2 , 100]);
    norm=1:maxLengthTransitionBeats*2;
    norm=[norm maxLengthTransitionBeats*2-1:-1:1];
    a=descriptors(1).energyStart(end-2*maxLengthTransitionBeats+1:end);
    b=descriptors(2).energyStart(1:maxLengthTransitionBeats*2);
    [cs,lags]=xcorr(a , b, maxLengthTransitionBeats*2-1) ; %, length(b)-1);
    %cs=cs(:)./norm(:)*max(norm)/4;
    cross.EnergyStart=cs(end-maxLengthTransitionBeats+1:end);
    %plotCorrelations(a,b,cross.EnergyStart, maxLengthTransitionBeats)
    %{
        figure
        hold on
        plot(a)
        plot(b,'g')
        plot(cross.EnergyStart,'r')
        plot(cs,'k')
    %}

    %% beatDuration
    a=descriptors(1).beatDuration(end-2*maxLengthTransitionBeats+1:end);
    b=descriptors(2).beatDuration(1:maxLengthTransitionBeats*2);
    cs=xcorr(a , b, maxLengthTransitionBeats*2-1);
    %cs=cs./norm(:)';
    cross.beatDuration=cs(end-maxLengthTransitionBeats+1:end);
    %plot(cross.beatDuration)
    %plotCorrelations(a,b,cross.beatDuration, maxLengthTransitionBeats)

    %% metrical Position
    metricalPosWeight=1/30;
    a=descriptors(1).metricalPos(end-2*maxLengthTransitionBeats+1:end);
    b=descriptors(2).metricalPos(1:maxLengthTransitionBeats*2);
    cs=xcorr(a , b, maxLengthTransitionBeats*2-1);
    cs=cs(:)./norm(:)*max(norm)*metricalPosWeight;
    cross.metricalPos=cs(end-maxLengthTransitionBeats+1:end);
    %plot(cross.metricalPos,'x-')
    %plot(cs,'x-')
    %plotCorrelations(a,b,cross.metricalPos, maxLengthTransitionBeats)

%     %% Chroma Features
%     %figure
%     %hold on
%     clear csbsChromaF;
% 
%     nFeatures=size(descriptors(1).bsChromaF,1);
%     colors=hsv(nFeatures);
%     for i=1:nFeatures
%         a=descriptors(1).bsChromaF(i,end-2*maxLengthTransitionBeats+1:end);
%         b=descriptors(2).bsChromaF(i,1:maxLengthTransitionBeats*2);
%         cs=xcorr(a , b, maxLengthTransitionBeats*2-1);
%         csbsChromaF(i,:)=cs(end-maxLengthTransitionBeats+1:end);
%         %plot(csbsChromaF(i,:), 'color', colors(i,:))
%     end
%     cross.Chroma=sum(csbsChromaF,1)/nFeatures;
%     %plot(cross.Chroma,'lineWidth',2,'color','k')
%     %%plot(prod(csbsChromaF,1)*10^(nFeatures/2),'lineWidth',2,'color','k')
% 
%     %% MFCC Features
% 
%     clear csbsMFCC;
% 
%     nFeatures=size(descriptors(1).bsMFCC,1);
%     colors=hsv(nFeatures);
%     for i=1:nFeatures
%         a=descriptors(1).bsMFCC(i,end-2*maxLengthTransitionBeats+1:end);
%         b=descriptors(2).bsMFCC(i,1:maxLengthTransitionBeats*2);
%         cs=xcorr(a , b, maxLengthTransitionBeats*2-1);
%         %cs=cs./norm*max(norm)/4;
%         csbsMFCC(i,:)=cs(end-maxLengthTransitionBeats+1:end);
%         %plot(csbsMFCC(i,:), 'color', colors(i,:))
%     end
%     cross.MFCC=sum(csbsMFCC,1)/nFeatures;
% 
%     %figure
%     %hold on
%     %plot(cross.MFCC,'lineWidth',2,'color','k')
%     %%plot(prod(csbsChromaF,1)*10^(nFeatures/2),'lineWidth',2,'color','k')

    %% compute final xcorr
    doPlot=0;
    if (doPlot)
        figure
        hold on
    end
    fieldNames=fieldnames(cross);
    colors=hsv(length(fieldnames(cross)));
    theSum=zeros(1,length(cross.metricalPos));
    for i=1:length(fieldnames(cross))
        aux=getfield(cross, fieldNames{i});
        aux=aux(:)';
        if (doPlot)
             plot(aux,'color', colors(i,:),'displayName',fieldNames{i})
        end
        theSum=theSum+aux;
    end
    theSum=theSum/length(fieldnames(cross));
    if (doPlot)
        plot(theSum,'k-x','lineWidth',2,'displayName','Weighted Corr. Sum.')
    end

    minTranDurationBeats=20;
    [val,trBeatIdx]=max(theSum(1:end-minTranDurationBeats));
    trBeatIdx=maxLengthTransitionBeats-trBeatIdx; %+minTranDurationBeats;
    
end

function plotCorrelations(a,b,coor, maxLength)
    %{
        figure
        hold on
        plot(a)
        plot(b,'g')
        corr=cross.metricalPos;
        plot(corr,'r')
        plot(cs,'k')
    %}

    
     [val,trBeatIdx]=max(corr);
     trBeatIdx=maxLength-trBeatIdx;
     figure
     plot(a,'x-')
     hold on
     plot([zeros(1, length(a)-trBeatIdx) b],'g-x')
    
end