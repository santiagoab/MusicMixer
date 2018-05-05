
addpath('C:\Users\SANTI\Desktop\MusicMixer\phaseVocoder')
%% 
doPlots=0;


%% -------- Load first song of the transition
fileName{1}=MixPlaylist(1).name;
waveIn{1} = MixPlaylist(1).wave;
sr(1) = MixPlaylist(1).sr;
%playListOut{1}.beat=0; %0 is the number of beats int the first transition
%playListOut{1}.filename=MixPlaylist(1).name; %0 is the number of beats int the first transition
descriptors(1).beatPos = MixPlaylist(1).beatPos;
descriptors(1).beatDuration = MixPlaylist(1).beatDuration;



%% clear var
clear lastTrack

for iTrack=2:length(MixPlaylist)
    if( iTrack>2)
        features(1)=features(2);
        waveIn{1}=waveIn{2};
        sr(1)=sr(2);
        fileName{1}=fileName{2};
        %clear descriptors(2)
    end

    %load the second track
    fileName{2}=MixPlaylist(iTrack).name
    %playListOut{iTrack}.filename=MixPlaylist(1).name; 
    disp(['Anayzing ' fileName{2}]);
    waveIn{2} = MixPlaylist(iTrack).wave;
    sr(2) = MixPlaylist(iTrack).sr;
    %features(2)=compBeatDescriptors(waveIn{2},sr(2));
    %descriptors(2).beatPos = MixPlaylist(iTrack).beatPos;
    %descriptors(2).beatDuration = MixPlaylist(iTrack).beatDuration;
    %descriptors(2).

    [trBeatIdx]=round(getBestTransitionBeatPos(MixPlaylist(iTrack-1:iTrack)));

    
    %playListOut{iTrack}.beat=trBeatIdx; %0 is the number of beats int the first transition

    %% ...

    %compute the stretch ratio (both tracks)
    avgBeatDur(1)=  (MixPlaylist(iTrack-1).beatPos(end-trBeatIdx) - MixPlaylist(iTrack-1).beatPos(end-trBeatIdx-10)) /10;   %avg. beat duration in the last 10 beats before the transition. (1st track)
    avgBeatDur(2)=  (MixPlaylist(iTrack).beatPos(trBeatIdx+10) - MixPlaylist(iTrack).beatPos(trBeatIdx))/10; %avg. beat dur. in the 10 beats fter the transition (2track)
    rampDuration = interp1([1,trBeatIdx],avgBeatDur,1:trBeatIdx);
    %rampPosition=cumsum(ampDuration)+MixPlaylist(iTrack-1).beatPos(end-trBeatIdx);
    if (doPlots)
        figure
        hold on
        plot(MixPlaylist(iTrack-1).beatDuration,'b-x')
        plot([nan(1, length(MixPlaylist(iTrack-1).beatDuration)-trBeatIdx+1) MixPlaylist(iTrack).beatDuration ],'g-x')
        plot([nan(length(MixPlaylist(iTrack-1).beatPos) - trBeatIdx +1,1)' rampDuration],'r-x' ) 
    end
    
    bBeatIdx(1)=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
    bBeatIdx(2)=1;
    
    %% Plot transition without beat-MAtching
    if (doPlots)
        figure;
        ax(1)=subplot(3,1,1);
        b=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=MixPlaylist(iTrack-1).beatPos(end)*sr(1);
        plot(b:e, waveIn{1}(b:e))
        hold on
        nBeats=length(MixPlaylist(iTrack-1).beatPos)-(bBeatIdx(1)-10)+1;
        plot([MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:end)*sr(1) ; descriptors(1).beatPos(bBeatIdx(1)-10:end)*sr(1)], [-1*ones(1,nBeats); ones(1,nBeats)])
        ax(2)=subplot(3,1,2);
        b=MixPlaylist(iTrack).beatPos(1)*sr(2);
        e=MixPlaylist(iTrack).beatPos(trBeatIdx+10)*sr(2);
        dispOffset=descriptors(1).beatPos(bBeatIdx(1))*sr(1);
        plot((b:e)+dispOffset, waveIn{1}(b:e))       
        hold on
        plot([MixPlaylist(iTrack).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset ; MixPlaylist(iTrack).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset], [-1*ones(1,nBeats); ones(1,nBeats)])
        linkaxes(ax);
    end

    %% stretch beats
    clear waveOut e b startSamples waveTransOut
    clear lastPh1 lastPh2;

    %rampDuration=interp1([1,trBeatIdx],[MixPlaylist(iTrack-1).beatDuration(bBeatIdx(1)) MixPlaylist(iTrack).beatDuration(bBeatIdx(2)) ],1:trBeatIdx);
    startSamples(1)=round(MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1));
    startSamples(2)= round(MixPlaylist(iTrack).beatPos(bBeatIdx(2))*sr(2));
    b(1)=startSamples(1);
    b(2)=startSamples(2);
    waveTransOut{1}=[]; %[waveIn{1}(1:b(1)-1)]';
    waveTransOut{2}=[];
    fftSize=1024;
    %if ratio--> inputDuration/targetDuration
    for ibeat=1:trBeatIdx-1
        eBeatIdx(1)=length(MixPlaylist(iTrack-1).beatPos)-(trBeatIdx-ibeat)+1;
        eBeatIdx(2)=ibeat+1;
        e(1)=min(length(waveIn{1}), round(MixPlaylist(iTrack-1).beatPos(eBeatIdx(1))*sr(1)));
        e(2)=min(length(waveIn{1}), round(MixPlaylist(iTrack).beatPos(eBeatIdx(2))*sr(2)));
        % for each beat in the transition, stretch the beats according to the ramp
        ratio(1)=MixPlaylist(iTrack-1).beatDuration(eBeatIdx(1)-1)/rampDuration(ibeat); %ratio source/target
        %ratio(2)=MixPlaylist(iTrack).beatDuration(ibeat)/rampDuration(ibeat);
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
        %ratio(2)=(length(buffer)/sr(1))/MixPlaylist(iTrack).beatDuration(ibeat);

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
   %  audiowrite(waveTransOut{1},sr(1), [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack-1) '.wav']);
   %  audiowrite(waveTransOut{2},44100, [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack) '.wav']);
   %  waveMix=[waveTransOut{1}+waveTransOut{2}];
   %  audiowrite(waveMix,sr(2), [dir 'transition' num2str(iTrack-1) 'mix.wav']);

    %% some plots..
    if (doPlots)
        figure;
        ax(1)=subplot(3,1,1);
        %first plot 10 beats before the transition
        b=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1);        
        plot(b:e, waveIn{1}(b:e)./max(waveIn{1}(b:e)),'Color',[0.8 0.8 0.8])
        hold on
        nBeats=11; %10 before plus the transition point
        t1=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:bBeatIdx(1))*sr(1);
        plot([ t1; t1], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then plot the transition        
        plot((1:length(waveTransOut{1}))+e, waveTransOut{1}./max(waveTransOut{1}),'Color',[0.5 0.5 0.5])
        plot((1:length(waveTransOut{1}))+e, highFadeOut,'r')        
        t2=([0 cumsum(rampDuration(1:end-1))]+MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %--Now the second song
        ax(2)=subplot(3,1,2);
        %plot then the transition
        plot((1:length(waveTransOut{2}))+e, waveTransOut{2}./max(waveTransOut{2}),'Color',[0.5 0.5 0.5])
        hold on
        plot((1:length(waveTransOut{1}))+e, highFadeIn,'r')        
        %t3=([0 cumsum(rampDuration(1:end))]+MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then the next 10 beats of the song
        b2=MixPlaylist(iTrack).beatPos(trBeatIdx+1)*sr(2);
        e2=MixPlaylist(iTrack).beatPos(trBeatIdx+11)*sr(2);
        dispOffset=t2(end) - b2; %MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1);       
        plot((b2:e2)+dispOffset, waveIn{2}(b2:e2)./max(waveIn{2}(b2:e2)),'Color',[0.8 0.8 0.8])       
        t3=MixPlaylist(iTrack).beatPos(trBeatIdx+1:trBeatIdx+11)*sr(2)+dispOffset;
        nBeats=length(t3);
        plot([t3 ; t3], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %Finally the tempo Curve
        beatDur=[ MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-9:bBeatIdx(1))] - [ MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:bBeatIdx(1)-1)];
        beatDur=[beatDur rampDuration(1:end)];
        aux=[ MixPlaylist(iTrack).beatPos(trBeatIdx+1:trBeatIdx+10)] - [ MixPlaylist(iTrack).beatPos(trBeatIdx:trBeatIdx+9)];
        beatDur=[beatDur aux];
        ax(3)=subplot(3,1,3)
        plot([t1(1:end-1) t2(1:end-1) t3] ,beatDur,'kx-')
        linkaxes(ax,'x');
    end
    
    
    %%
    if (exist('lastTrack','var'))
        waveOut=[lastTrack.prev waveIn{1}(lastTrack.transitionLastSample+1:startSamples(1))' waveTransOut{1}+waveTransOut{2}];
        firstTr=lastTrack.tranNewBeatPos(1:end-1);
        middle=MixPlaylist(iTrack-1).beatPos(lastTrack.transitionLastBeat+1:end-trBeatIdx+1);
        middle=middle-MixPlaylist(iTrack-1).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
        secondTr=cumsum(rampDuration(1:end-1))+middle(end);
        beatPosOut=[firstTr middle secondTr];
    else %first and second songs
        waveOut=[waveIn{1}(1:startSamples(1))' waveTransOut{1}+waveTransOut{2}];
        idxFirstBeatInTransition=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
        trStart=MixPlaylist(iTrack-1).beatPos(idxFirstBeatInTransition); %it is the same as start of the first beat in the transition 
        beatPosOut=[MixPlaylist(iTrack-1).beatPos(1:idxFirstBeatInTransition) cumsum(rampDuration(1:end-1))+ trStart]; %remove last element of rampDuration as is the end of the last beat
        
        
    end
    
%     %compute Beat similarity excluding transitions
%     tr1Beat=0;
%     if(exist('lastTrack','var'))
%         tr1Beat=lastTrack.transitionLastBeat;
%     end
%     [closest]= compBeatSimilarity(MixPlaylist(iTrack-1), [tr1Beat trBeatIdx]);
%     
%     % save JSON
%     jsonString=savejson('',closest);
%     fid=fopen([dir fileName{1} '_p_closest.json'],'w');
%     fprintf(fid, '%s', jsonString);
%     fclose(fid);  
%     
     %save processed track and beatPos
     
%     savetoWavesurfer(dir, [fileName{1} '_p'], beatPosOut);
%     
%     % save to JSON
%     jsonString=savejson('',round(beatPosOut*44100));
%     fid=fopen([dir fileName{1} '_p_beats.json'],'w');
%     fprintf(fid, '%s', jsonString);
%     fclose(fid);
    
    %% some plots
    if (doPlots)
        %first track
        firstBeatTr=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
        %beatPos(1)=[MixPlaylist(iTrack-1).beatPos(1:firstBeatTr-1) beatPosOut];
        plot(beatPosOut,'g-x') %first song
        hold on       
        %second track
        t=[1:length(MixPlaylist(iTrack).beatPos)]+firstBeatTr-1;
        accTime=cumsum(rampDuration(1:end));
        firstPart=[0 accTime(1:end-1)]; 
        secondPart=MixPlaylist(iTrack).beatPos(trBeatIdx+1:end); 
        secondPart=secondPart - MixPlaylist(iTrack).beatPos(trBeatIdx+1) + accTime(end); %update offset
        y=[firstPart  secondPart];  %+ MixPlaylist(iTrack).beatPos(1)
        y=y+MixPlaylist(iTrack-1).beatPos(firstBeatTr);
        plot(t,y,'x-') %position of the beats of the second song (respet the zero seconds at the begining of the first song)
    end
    lastTrack.prev = waveOut;
    lastTrack.inTransition=waveTransOut{2};
    lastTrack.transitionLastSample=e(2); %trPosSec(2);
    lastTrack.transitionLastBeat=trBeatIdx; 
    lastTrack.tranNewBeatPos=[0 cumsum(rampDuration(1:end))]+MixPlaylist(iTrack).beatPos(1); %zero and all elements of rampDuration. So the last value is the position of the start beat after the transition. 
    

end


%waveOut=[lastTrack.inTransition waveIn{2}(lastTrack.transitionLastSample+1:startSamples(2))'];
waveOut=[lastTrack.prev waveIn{2}(e(2):end)'];

firstTr=lastTrack.tranNewBeatPos(1:end-1);
middle=MixPlaylist(iTrack).beatPos(lastTrack.transitionLastBeat+1:end);
middle=middle-MixPlaylist(iTrack).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
beatPosOut=[firstTr middle];

%save processed track and beatPos
audiowrite(['3songMix.wav'], waveOut,44100 );

addpath('C:\Users\SANTI\Desktop\MusicMixer\phaseVocoder')
%% 
doPlots=0;


%% -------- Load first song of the transition
fileName{1}=MixPlaylist(1).name;
waveIn{1} = MixPlaylist(1).wave;
sr(1) = MixPlaylist(1).sr;
%playListOut{1}.beat=0; %0 is the number of beats int the first transition
%playListOut{1}.filename=MixPlaylist(1).name; %0 is the number of beats int the first transition
descriptors(1).beatPos = MixPlaylist(1).beatPos;
descriptors(1).beatDuration = MixPlaylist(1).beatDuration;



%% clear var
clear lastTrack

for iTrack=2:length(MixPlaylist)
    if( iTrack>2)
        features(1)=features(2);
        waveIn{1}=waveIn{2};
        sr(1)=sr(2);
        fileName{1}=fileName{2};
        %clear descriptors(2)
    end

    %load the second track
    fileName{2}=MixPlaylist(iTrack).name
    %playListOut{iTrack}.filename=MixPlaylist(1).name; 
    disp(['Anayzing ' fileName{2}]);
    waveIn{2} = MixPlaylist(iTrack).wave;
    sr(2) = MixPlaylist(iTrack).sr;
    %features(2)=compBeatDescriptors(waveIn{2},sr(2));
    %descriptors(2).beatPos = MixPlaylist(iTrack).beatPos;
    %descriptors(2).beatDuration = MixPlaylist(iTrack).beatDuration;
    %descriptors(2).

    [trBeatIdx]=round(getBestTransitionBeatPos(MixPlaylist(iTrack-1:iTrack)));

    
    %playListOut{iTrack}.beat=trBeatIdx; %0 is the number of beats int the first transition

    %% ...

    %compute the stretch ratio (both tracks)
    avgBeatDur(1)=  (MixPlaylist(iTrack-1).beatPos(end-trBeatIdx) - MixPlaylist(iTrack-1).beatPos(end-trBeatIdx-10)) /10;   %avg. beat duration in the last 10 beats before the transition. (1st track)
    avgBeatDur(2)=  (MixPlaylist(iTrack).beatPos(trBeatIdx+10) - MixPlaylist(iTrack).beatPos(trBeatIdx))/10; %avg. beat dur. in the 10 beats fter the transition (2track)
    rampDuration = interp1([1,trBeatIdx],avgBeatDur,1:trBeatIdx);
    %rampPosition=cumsum(ampDuration)+MixPlaylist(iTrack-1).beatPos(end-trBeatIdx);
    if (doPlots)
        figure
        hold on
        plot(MixPlaylist(iTrack-1).beatDuration,'b-x')
        plot([nan(1, length(MixPlaylist(iTrack-1).beatDuration)-trBeatIdx+1) MixPlaylist(iTrack).beatDuration ],'g-x')
        plot([nan(length(MixPlaylist(iTrack-1).beatPos) - trBeatIdx +1,1)' rampDuration],'r-x' ) 
    end
    
    bBeatIdx(1)=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
    bBeatIdx(2)=1;
    
    %% Plot transition without beat-MAtching
    if (doPlots)
        figure;
        ax(1)=subplot(3,1,1);
        b=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=MixPlaylist(iTrack-1).beatPos(end)*sr(1);
        plot(b:e, waveIn{1}(b:e))
        hold on
        nBeats=length(MixPlaylist(iTrack-1).beatPos)-(bBeatIdx(1)-10)+1;
        plot([MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:end)*sr(1) ; descriptors(1).beatPos(bBeatIdx(1)-10:end)*sr(1)], [-1*ones(1,nBeats); ones(1,nBeats)])
        ax(2)=subplot(3,1,2);
        b=MixPlaylist(iTrack).beatPos(1)*sr(2);
        e=MixPlaylist(iTrack).beatPos(trBeatIdx+10)*sr(2);
        dispOffset=descriptors(1).beatPos(bBeatIdx(1))*sr(1);
        plot((b:e)+dispOffset, waveIn{1}(b:e))       
        hold on
        plot([MixPlaylist(iTrack).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset ; MixPlaylist(iTrack).beatPos(1:trBeatIdx+10)*sr(1)+dispOffset], [-1*ones(1,nBeats); ones(1,nBeats)])
        linkaxes(ax);
    end

    %% stretch beats
    clear waveOut e b startSamples waveTransOut
    clear lastPh1 lastPh2;

    %rampDuration=interp1([1,trBeatIdx],[MixPlaylist(iTrack-1).beatDuration(bBeatIdx(1)) MixPlaylist(iTrack).beatDuration(bBeatIdx(2)) ],1:trBeatIdx);
    startSamples(1)=round(MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1));
    startSamples(2)= round(MixPlaylist(iTrack).beatPos(bBeatIdx(2))*sr(2));
    b(1)=startSamples(1);
    b(2)=startSamples(2);
    waveTransOut{1}=[]; %[waveIn{1}(1:b(1)-1)]';
    waveTransOut{2}=[];
    fftSize=1024;
    %if ratio--> inputDuration/targetDuration
    for ibeat=1:trBeatIdx-1
        eBeatIdx(1)=length(MixPlaylist(iTrack-1).beatPos)-(trBeatIdx-ibeat)+1;
        eBeatIdx(2)=ibeat+1;
        e(1)=min(length(waveIn{1}), round(MixPlaylist(iTrack-1).beatPos(eBeatIdx(1))*sr(1)));
        e(2)=min(length(waveIn{1}), round(MixPlaylist(iTrack).beatPos(eBeatIdx(2))*sr(2)));
        % for each beat in the transition, stretch the beats according to the ramp
        ratio(1)=MixPlaylist(iTrack-1).beatDuration(eBeatIdx(1)-1)/rampDuration(ibeat); %ratio source/target
        %ratio(2)=MixPlaylist(iTrack).beatDuration(ibeat)/rampDuration(ibeat);
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
        %ratio(2)=(length(buffer)/sr(1))/MixPlaylist(iTrack).beatDuration(ibeat);

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
   %  audiowrite(waveTransOut{1},sr(1), [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack-1) '.wav']);
   %  audiowrite(waveTransOut{2},44100, [dir 'transition' num2str(iTrack-1) 'track' num2str(iTrack) '.wav']);
   %  waveMix=[waveTransOut{1}+waveTransOut{2}];
   %  audiowrite(waveMix,sr(2), [dir 'transition' num2str(iTrack-1) 'mix.wav']);

    %% some plots..
    if (doPlots)
        figure;
        ax(1)=subplot(3,1,1);
        %first plot 10 beats before the transition
        b=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10)*sr(1);
        e=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1);        
        plot(b:e, waveIn{1}(b:e)./max(waveIn{1}(b:e)),'Color',[0.8 0.8 0.8])
        hold on
        nBeats=11; %10 before plus the transition point
        t1=MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:bBeatIdx(1))*sr(1);
        plot([ t1; t1], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then plot the transition        
        plot((1:length(waveTransOut{1}))+e, waveTransOut{1}./max(waveTransOut{1}),'Color',[0.5 0.5 0.5])
        plot((1:length(waveTransOut{1}))+e, highFadeOut,'r')        
        t2=([0 cumsum(rampDuration(1:end-1))]+MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %--Now the second song
        ax(2)=subplot(3,1,2);
        %plot then the transition
        plot((1:length(waveTransOut{2}))+e, waveTransOut{2}./max(waveTransOut{2}),'Color',[0.5 0.5 0.5])
        hold on
        plot((1:length(waveTransOut{1}))+e, highFadeIn,'r')        
        %t3=([0 cumsum(rampDuration(1:end))]+MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)))*sr(1);
        nBeats=length(t2);
        plot([t2 ; t2], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        %then the next 10 beats of the song
        b2=MixPlaylist(iTrack).beatPos(trBeatIdx+1)*sr(2);
        e2=MixPlaylist(iTrack).beatPos(trBeatIdx+11)*sr(2);
        dispOffset=t2(end) - b2; %MixPlaylist(iTrack-1).beatPos(bBeatIdx(1))*sr(1);       
        plot((b2:e2)+dispOffset, waveIn{2}(b2:e2)./max(waveIn{2}(b2:e2)),'Color',[0.8 0.8 0.8])       
        t3=MixPlaylist(iTrack).beatPos(trBeatIdx+1:trBeatIdx+11)*sr(2)+dispOffset;
        nBeats=length(t3);
        plot([t3 ; t3], [-1*ones(1,nBeats); ones(1,nBeats)],'k:')
        
        %Finally the tempo Curve
        beatDur=[ MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-9:bBeatIdx(1))] - [ MixPlaylist(iTrack-1).beatPos(bBeatIdx(1)-10:bBeatIdx(1)-1)];
        beatDur=[beatDur rampDuration(1:end)];
        aux=[ MixPlaylist(iTrack).beatPos(trBeatIdx+1:trBeatIdx+10)] - [ MixPlaylist(iTrack).beatPos(trBeatIdx:trBeatIdx+9)];
        beatDur=[beatDur aux];
        ax(3)=subplot(3,1,3)
        plot([t1(1:end-1) t2(1:end-1) t3] ,beatDur,'kx-')
        linkaxes(ax,'x');
    end
    
    
    %%
    if (exist('lastTrack','var'))
        waveOut=[lastTrack.prev waveIn{1}(lastTrack.transitionLastSample+1:startSamples(1))' waveTransOut{1}+waveTransOut{2}];
        firstTr=lastTrack.tranNewBeatPos(1:end-1);
        middle=MixPlaylist(iTrack-1).beatPos(lastTrack.transitionLastBeat+1:end-trBeatIdx+1);
        middle=middle-MixPlaylist(iTrack-1).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
        secondTr=cumsum(rampDuration(1:end-1))+middle(end);
        beatPosOut=[firstTr middle secondTr];
    else %first and second songs
        waveOut=[waveIn{1}(1:startSamples(1))' waveTransOut{1}+waveTransOut{2}];
        idxFirstBeatInTransition=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
        trStart=MixPlaylist(iTrack-1).beatPos(idxFirstBeatInTransition); %it is the same as start of the first beat in the transition 
        beatPosOut=[MixPlaylist(iTrack-1).beatPos(1:idxFirstBeatInTransition) cumsum(rampDuration(1:end-1))+ trStart]; %remove last element of rampDuration as is the end of the last beat
        
        
    end
    
%     %compute Beat similarity excluding transitions
%     tr1Beat=0;
%     if(exist('lastTrack','var'))
%         tr1Beat=lastTrack.transitionLastBeat;
%     end
%     [closest]= compBeatSimilarity(MixPlaylist(iTrack-1), [tr1Beat trBeatIdx]);
%     
%     % save JSON
%     jsonString=savejson('',closest);
%     fid=fopen([dir fileName{1} '_p_closest.json'],'w');
%     fprintf(fid, '%s', jsonString);
%     fclose(fid);  
%     
     %save processed track and beatPos
     
%     savetoWavesurfer(dir, [fileName{1} '_p'], beatPosOut);
%     
%     % save to JSON
%     jsonString=savejson('',round(beatPosOut*44100));
%     fid=fopen([dir fileName{1} '_p_beats.json'],'w');
%     fprintf(fid, '%s', jsonString);
%     fclose(fid);
    
    %% some plots
    if (doPlots)
        %first track
        firstBeatTr=length(MixPlaylist(iTrack-1).beatPos)-trBeatIdx+1;
        %beatPos(1)=[MixPlaylist(iTrack-1).beatPos(1:firstBeatTr-1) beatPosOut];
        plot(beatPosOut,'g-x') %first song
        hold on       
        %second track
        t=[1:length(MixPlaylist(iTrack).beatPos)]+firstBeatTr-1;
        accTime=cumsum(rampDuration(1:end));
        firstPart=[0 accTime(1:end-1)]; 
        secondPart=MixPlaylist(iTrack).beatPos(trBeatIdx+1:end); 
        secondPart=secondPart - MixPlaylist(iTrack).beatPos(trBeatIdx+1) + accTime(end); %update offset
        y=[firstPart  secondPart];  %+ MixPlaylist(iTrack).beatPos(1)
        y=y+MixPlaylist(iTrack-1).beatPos(firstBeatTr);
        plot(t,y,'x-') %position of the beats of the second song (respet the zero seconds at the begining of the first song)
    end
    lastTrack.prev = waveOut;
    lastTrack.inTransition=waveTransOut{2};
    lastTrack.transitionLastSample=e(2); %trPosSec(2);
    lastTrack.transitionLastBeat=trBeatIdx; 
    lastTrack.tranNewBeatPos=[0 cumsum(rampDuration(1:end))]+MixPlaylist(iTrack).beatPos(1); %zero and all elements of rampDuration. So the last value is the position of the start beat after the transition. 
    

end


%waveOut=[lastTrack.inTransition waveIn{2}(lastTrack.transitionLastSample+1:startSamples(2))'];
waveOut=[lastTrack.prev waveIn{2}(e(2):end)'];

firstTr=lastTrack.tranNewBeatPos(1:end-1);
middle=MixPlaylist(iTrack).beatPos(lastTrack.transitionLastBeat+1:end);
middle=middle-MixPlaylist(iTrack).beatPos(lastTrack.transitionLastBeat+1)+lastTrack.tranNewBeatPos(end); %change offset given by last stretch
beatPosOut=[firstTr middle];

%save processed track and beatPos
audiowrite(['FinalMix.wav'], waveOut,44100 );
