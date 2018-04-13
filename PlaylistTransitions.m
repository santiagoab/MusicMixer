clear 
%% 
doPlots=0;
startbpm=120;

% load playlist
%dir='D:\projects\FASTQMUL\silentDisco\sounds\';
dir='D:\projects\FASTQMUL\ionicprojects\testBlank\www\';
playlistName='playlist.txt';
playlist = loadjson([dir playlistName]);
%playListOut=loadjson([dir playlistName]);

clear waveIn
% -------- Load first song of the transition
fileName{1}=playlist{1};
playListOut{1}.beat=0; %0 is the number of beats int the first transition
playListOut{1}.filename=playlist{1}; %0 is the number of beats int the first transition

%[closest, beatPos, bsChromaF, bsMFCC, energyTime, energyStart, metrical]= compBeatSimilarity(dir, fileName,0);
disp(['Anayzing ' fileName{1}]);
[waveIn{1},sr(1)] = wavread([dir fileName{1} '.wav']);
descriptors(1)=compBeatDescriptors(waveIn{1},sr(1),120);


%% for each song ...
clear lastTrack;    
for iTrack=2:length(playlist)
    if( iTrack>2)
        descriptors(1)=descriptors(2);
        waveIn{1}=waveIn{2};
        sr(1)=sr(2);
        fileName{1}=fileName{2};
        %clear descriptors(2)
    end

    %load the second track
    fileName{2}=playlist{iTrack}; %'03-FreakFandango-GypsySong'; %'03-FreakFandango-GypsySong'; %'01- Ciaocarlia'; %
    playListOut{iTrack}.filename=playlist{iTrack}; %0 is the number of beats int the first transition
    disp(['Anayzing ' fileName{2}]);
    [waveIn{2},sr(2)] = wavread([dir fileName{2} '.wav']);
    descriptors(2)=compBeatDescriptors(waveIn{2},sr(2));

    [trBeatIdx]=getBestTransitionBeatPos(descriptors);

    
    playListOut{iTrack}.beat=trBeatIdx; %0 is the number of beats int the first transition

    %% ...
    %trBeatDurs(1,:)=descriptors(1).beatDuration(end-trBeatIdx+1:end);
    %trBeatDurs(2,:)=descriptors(2).beatDuration(1:trBeatIdx);
    %plot(trBeatDurs','x-')

    %compute the stretch ratio (both tracks)
    avgBeatDur(1)=  (descriptors(1).beatPos(end-trBeatIdx) - descriptors(1).beatPos(end-trBeatIdx-10)) /10;   %avg. beat duration in the last 10 beats before the transition. (1st track)
    avgBeatDur(2)=  (descriptors(2).beatPos(trBeatIdx+10) - descriptors(2).beatPos(trBeatIdx))/10; %avg. beat dur. in the 10 beats fter the transition (2track)
    rampDuration = interp1([1,trBeatIdx],avgBeatDur,1:trBeatIdx);
    %rampPosition=cumsum(ampDuration)+descriptors(1).beatPos(end-trBeatIdx);
    if (doPlots)
        figure
        hold on
        plot(descriptors(1).beatDuration,'b-x')
        plot([nan(1, length(descriptors(1).beatDuration)-trBeatIdx+1) descriptors(2).beatDuration ],'g-x')
        plot([nan(length(descriptors(1).beatPos) - trBeatIdx +1,1)' rampDuration],'r-x' ) 
    end
    
    bBeatIdx(1)=length(descriptors(1).beatPos)-trBeatIdx+1;
    bBeatIdx(2)=1;
    
    %% Plot transition without beat-MAtching
    if (1)
        figure;
        ax(1)=subplot(3,1,1);
        b=descriptors(1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=descriptors(1).beatPos(end)*sr(1);
        plot(b:e, waveIn{1}(b:e))
        hold on
        nBeats=length(descriptors(1).beatPos)-(bBeatIdx(1)-10)+1;
        plot([descriptors(1).beatPos(bBeatIdx(1)-10:end)*sr(1) ; descriptors(1).beatPos(bBeatIdx(1)-10:end)*sr(1)], [-1*ones(1,nBeats); ones(1,nBeats)])
        ax(2)=subplot(3,1,2);
        b=descriptors(2).beatPos(1)*sr(2);
        e=descriptors(2).beatPos(trBeatIdx+10)*sr(2);
        dispOffset=descriptors(1).beatPos(bBeatIdx(1))*sr(1);
        plot((b:e)+dispOffset, waveIn{1}(b:e))       
        hold on
        plot([descriptors(2).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset ; descriptors(2).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset], [-1*ones(1,nBeats); ones(1,nBeats)])
        linkaxes(ax);
    end

    %% stretch beats
    clear waveOut e b startSamples waveTransOut
    clear lastPh1 lastPh2;

    %rampDuration=interp1([1,trBeatIdx],[descriptors(1).beatDuration(bBeatIdx(1)) descriptors(2).beatDuration(bBeatIdx(2)) ],1:trBeatIdx);
    startSamples(1)=round(descriptors(1).beatPos(bBeatIdx(1))*sr(1));
    startSamples(2)= round(descriptors(2).beatPos(bBeatIdx(2))*sr(2));
    b(1)=startSamples(1);
    b(2)=startSamples(2);
    waveTransOut{1}=[]; %[waveIn{1}(1:b(1)-1)]';
    waveTransOut{2}=[];
    fftSize=1024;
    %if ratio--> inputDuration/targetDuration
    for ibeat=1:trBeatIdx-1
        eBeatIdx(1)=length(descriptors(1).beatPos)-(trBeatIdx-ibeat)+1;
        eBeatIdx(2)=ibeat+1;
        e(1)=min(length(waveIn{1}), round(descriptors(1).beatPos(eBeatIdx(1))*sr(1)));
        e(2)=min(length(waveIn{1}), round(descriptors(2).beatPos(eBeatIdx(2))*sr(2)));
        % for each beat in the transition, stretch the beats according to the ramp
        ratio(1)=descriptors(1).beatDuration(eBeatIdx(1)-1)/rampDuration(ibeat); %ratio source/target
        %ratio(2)=descriptors(2).beatDuration(ibeat)/rampDuration(ibeat);
        % stretch waveIn{1}
        buffer=waveIn{1}(b(1):e(1));
        lengthBufferOut=length(buffer)/ratio(1); % lengthBufferOut=length(buffer)/1.5;
        if (~exist('lastPh1','var'))
            [bufferOut, lastPh1]=pvocLengthNoRatio(buffer,lengthBufferOut,fftSize); %  [bufferOut, lastPh1]=pvoc(buffer,1.5,fftSize);
        else
            [bufferOut, lastPh1]=pvocLengthNoRatio(buffer,lengthBufferOut,fftSize, lastPh1);
        end
        waveTransOut{1}=[waveTransOut{1} bufferOut'];
        b(1)=e(1)+1;
        %ratio(2)=(length(buffer)/sr(1))/descriptors(2).beatDuration(ibeat);

        % stretch waveIn{2}
        buffer=waveIn{2}(b(2):e(2));
        if (~exist('lastPh2','var'))
            [bufferOut, lastPh2]=pvocLengthNoRatio(buffer,lengthBufferOut,fftSize);
        else
            [bufferOut, lastPh2]=pvocLengthNoRatio(buffer,lengthBufferOut,fftSize, lastPh2);
        end
        waveTransOut{2}=[waveTransOut{2} bufferOut'];
        b(2)=e(2)+1;
    end
    t=(1:length(waveTransOut{1}))/length(waveTransOut{1});
    t=t*2-1; %normalize t between -1 and 1 to make the crossfades equal power (sqrt)
    highFadeIn= sqrt((1/2)*(1+t));
    highFadeOut= sqrt((1/2)*(1-t));
    %highFadeOut=1-(logspace(0,1,length(waveTransOut{1}))-1)/9;
    %highFadeIn=1-(logspace(1,0,length(waveTransOut{2}))-1)/9;
    %plot(highFadeOut)
    %plot(sqrt(highFadeOut.^2+highFadeIn.^2)) %make sure that the sum is zero??
    
    %%
    waveTransOut{1}=waveTransOut{1}.* highFadeOut;  %linspace(1,0,length(waveTransOut{1}));
    waveTransOut{2}=waveTransOut{2}.* highFadeIn;  %linspace(0,1,length(waveTransOut{2}));
    
    %save the transitions
     wavwrite(waveTransOut{1},sr(1), [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack-1) '.wav']);
     wavwrite(waveTransOut{2},44100, [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack) '.wav']);
     waveMix=[waveTransOut{1}+waveTransOut{2}];
     wavwrite(waveMix,sr(2), [dir 'transition' num2str(iTrack-1) 'mix.wav']);

    %% some plots..
    if (1)
        figure;
        ax(1)=subplot(3,1,1);
        %first plot 10 beats before the transition
        b=descriptors(1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=descriptors(1).beatPos(bBeatIdx(1))*sr(1);        
        plot(b:e, waveIn{1}(b:e)./max(waveIn{1}(b:e)),'Color',[0.8 0.8 0.8])
        hold on
        nBeats=11; %10 before plus the transition point
        t1=descriptors(1).beatPos(bBeatIdx(1)-10:bBeatIdx(1))*sr(1);
        plot([ t1; t1], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then plot the transition        
        plot((1:length(waveTransOut{1}))+e, waveTransOut{1}./max(waveTransOut{1}),'Color',[0.5 0.5 0.5])
        plot((1:length(waveTransOut{1}))+e, highFadeOut,'r')        
        t2=([0 cumsum(rampDuration(1:end-1))]+descriptors(1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %--Now the second song
        ax(2)=subplot(3,1,2);
        %plot then the transition
        plot((1:length(waveTransOut{2}))+e, waveTransOut{2}./max(waveTransOut{2}),'Color',[0.5 0.5 0.5])
        hold on
        plot((1:length(waveTransOut{1}))+e, highFadeIn,'r')        
        %t3=([0 cumsum(rampDuration(1:end))]+descriptors(1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then the next 10 beats of the song
        b2=descriptors(2).beatPos(trBeatIdx+1)*sr(2);
        e2=descriptors(2).beatPos(trBeatIdx+11)*sr(2);
        dispOffset=t2(end) - b2; %descriptors(1).beatPos(bBeatIdx(1))*sr(1);       
        plot((b2:e2)+dispOffset, waveIn{2}(b2:e2)./max(waveIn{2}(b2:e2)),'Color',[0.8 0.8 0.8])       
        t3=descriptors(2).beatPos(trBeatIdx+1:trBeatIdx+11)*sr(2)+dispOffset;
        nBeats=length(t3);
        plot([t3 ; t3], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %Finally the tempo Curve
        beatDur=[ descriptors(1).beatPos(bBeatIdx(1)-9:bBeatIdx(1))] - [ descriptors(1).beatPos(bBeatIdx(1)-10:bBeatIdx(1)-1)];
        beatDur=[beatDur rampDuration(1:end)];
        aux=[ descriptors(2).beatPos(trBeatIdx+1:trBeatIdx+10)] - [ descriptors(2).beatPos(trBeatIdx:trBeatIdx+9)];
        beatDur=[beatDur aux];
        ax(3)=subplot(3,1,3)
        plot([t1(1:end-1) t2(1:end-1) t3] ,beatDur,'kx-')
        linkaxes(ax,'x');
    end
    
    
    %%
    if (exist('lastTrack','var'))
        waveOut=[lastTrack.inTransition waveIn{1}(lastTrack.transitionLastSample+1:startSamples(1))' waveTransOut{1}];
        firstTr=lastTrack.tranNewBeatPos(1:end-1);
        middle=descriptors(1).beatPos(lastTrack.transitionLastBeat+1:end-trBeatIdx+1);
        middle=middle-descriptors(1).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
        secondTr=cumsum(rampDuration(1:end-1))+middle(end);
        beatPosOut=[firstTr middle secondTr];
    else
        waveOut=[waveIn{1}(1:startSamples(1))' waveTransOut{1}];
        idxFirstBeatInTransition=length(descriptors(1).beatPos)-trBeatIdx+1;
        trStart=descriptors(1).beatPos(idxFirstBeatInTransition); %it is the same as start of the first beat in the transition 
        beatPosOut=[descriptors(1).beatPos(1:idxFirstBeatInTransition) cumsum(rampDuration(1:end-1))+ trStart]; %remove last element of rampDuration as is the end of the last beat
    end
    
    %compute Beat similarity excluding transitions
    tr1Beat=0;
    if(exist('lastTrack','var'))
        tr1Beat=lastTrack.transitionLastBeat;
    end
    [closest]= compBeatSimilarity(descriptors(1), [tr1Beat trBeatIdx]);
    
    % save JSON
    jsonString=savejson('',closest);
    fid=fopen([dir fileName{1} '_p_closest.json'],'w');
    fprintf(fid, '%s', jsonString);
    fclose(fid);  
    
    %save processed track and beatPos
    wavwrite(waveOut,44100, [dir fileName{1} '_p.wav']);
    savetoWavesurfer(dir, [fileName{1} '_p'], beatPosOut);
    
    % save to JSON
    jsonString=savejson('',round(beatPosOut*44100));
    fid=fopen([dir fileName{1} '_p_beats.json'],'w');
    fprintf(fid, '%s', jsonString);
    fclose(fid);
    
    %% some plots
    if (doPlots)
        %first track
        firstBeatTr=length(descriptors(1).beatPos)-trBeatIdx+1;
        %beatPos(1)=[descriptors(1).beatPos(1:firstBeatTr-1) beatPosOut];
        plot(beatPosOut,'g-x') %first song
        hold on       
        %second track
        t=[1:length(descriptors(2).beatPos)]+firstBeatTr-1;
        accTime=cumsum(rampDuration(1:end));
        firstPart=[0 accTime(1:end-1)]; 
        secondPart=descriptors(2).beatPos(trBeatIdx+1:end); 
        secondPart=secondPart - descriptors(2).beatPos(trBeatIdx+1) + accTime(end); %update offset
        y=[firstPart  secondPart];  %+ descriptors(2).beatPos(1)
        y=y+descriptors(1).beatPos(firstBeatTr);
        plot(t,y,'x-') %position of the beats of the second song (respet the zero seconds at the begining of the first song)
    end

    lastTrack.inTransition=waveTransOut{2};
    lastTrack.transitionLastSample=e(2); %trPosSec(2);
    lastTrack.transitionLastBeat=trBeatIdx; 
    lastTrack.tranNewBeatPos=[0 cumsum(rampDuration(1:end))]+descriptors(2).beatPos(1); %zero and all elements of rampDuration. So the last value is the position of the start beat after the transition. 
    

end

%save last Track
% descriptors(1)=descriptors(2);
% waveIn{1}=waveIn{2};
% sr(1)=sr(2);
% fileName{1}=fileName{2};

%waveOut=[lastTrack.inTransition waveIn{2}(lastTrack.transitionLastSample+1:startSamples(2))'];
waveOut=[waveTransOut{2} waveIn{2}(e(2):end)'];
firstTr=lastTrack.tranNewBeatPos(1:end-1);
middle=descriptors(2).beatPos(lastTrack.transitionLastBeat+1:end);
middle=middle-descriptors(2).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
beatPosOut=[firstTr middle];

%save processed track and beatPos
wavwrite(waveOut,44100, [dir fileName{2} '_p.wav']);
savetoWavesurfer(dir, [fileName{2} '_p'], beatPosOut);
% save to JSON
jsonString=savejson('',round(beatPosOut*44100));
fid=fopen([dir fileName{2} '_p_beats.json'],'w');
fprintf(fid, '%s', jsonString);
fclose(fid);

%compute Beat similarity excluding transitions
[closest]= compBeatSimilarity(descriptors(2), [trBeatIdx 0]);

% save JSON
jsonString=savejson('',closest);
fid=fopen([dir fileName{2} '_p_closest.json'],'w');
fprintf(fid, '%s', jsonString);
fclose(fid);  

    
%save playlist with iformation about the length (in beats) of the
% save to JSON
jsonString=savejson('',playListOut);
fid=fopen([dir 'playlist.json'],'w');
fprintf(fid, '%s', jsonString);
fclose(fid);



