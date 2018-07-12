function MixPlaylist = dataProcessing(yCoordinates, features, NumSongs)

%Resample the yCoordinates from the hand drawing to the same length of the
%number of desired songs, thus each value of the yCoordinates corresponding to the
% closest normalised energy of a song.
newY = resample(yCoordinates,NumSongs,length(yCoordinates));

TempFeatures = features; %temporary TempFeatures so the originals are not overwritten

% energy range for songs in a position. Range usually 0.1
range = 0.1;

%% weights for song parameters / next song decision
wegy = 1; %energy is normalized
wbpm = 1; %bpm is normalized
wkey = 1; %key is normalized
wmode = 1; %mode is normalized

%%

%for each position,(except the first one), in the drawing select closest songs by energy to the drawn curve within an energy range
for i=1:length(newY)
    
    if i == 1  %Condition for the first song
       
    [~, position{i}] = min(abs([TempFeatures.powerNormalized]-newY(i)));  %return for the first song index of closest song by energy to first drawing point
    
    MixPlaylist(i) = TempFeatures(position{i}); 
    TempFeatures(position{i}).powerNormalized = NaN; %Null power value to prevent using the same song again
    
    else %after the first song is saved      
    currentPos = position{i-1};
    
    
    [~,position{i}]=find(abs([TempFeatures.powerNormalized]-newY(i))<range);  % range usually 0.1
    
           if isempty(position{i}) == 0 
                   
                    currentSong = currentPos;
                    currentSongEner = TempFeatures(currentSong).powerNormalized;
                    currentSongBpm = TempFeatures(currentSong).bpmNormalized;
                    currentSongKey = TempFeatures(currentSong).keyNormalized;
                    currentSongMode = TempFeatures(currentSong).modeNormalized;
                    
                    nextPos = position{i};

                    for y = 1:length(nextPos)
                    % get next pos attributes 
                    nextSong = nextPos(y);
                    nextSongEner = TempFeatures(nextSong).powerNormalized;
                    nextSongBpm = TempFeatures(nextSong).bpmNormalized;
                    nextSongKey = TempFeatures(nextSong).keyNormalized;
                    nextSongMode = TempFeatures(nextSong).modeNormalized;

                    %% Compute weighted distances between such songs 
                    
                    dist(y) = (abs(currentSongEner - nextSongEner) * wegy) + (abs(currentSongBpm - nextSongBpm) * wbpm) + (abs(currentSongKey - nextSongKey) * wkey) + (abs(currentSongMode - nextSongMode) * wmode); 

                   
                    end
                    [~, index] = min(dist);
                    
                    MixPlaylist(i) = TempFeatures(nextPos(index));
                    TempFeatures(nextPos(index)).powerNormalized = NaN; %Null power value to prevent using the same song again
                    
                    dist = [];

            % if there are no songs within energy range, get closest song
            else isempty(position{i})
                [~, position{i}] = min(abs([TempFeatures.powerNormalized]-newY(i)));
                MixPlaylist(i) = TempFeatures(position{i}); 
                TempFeatures(position{i}).powerNormalized = NaN; %Null power value to prevent using the same song again
            end
        
    end
          
end

subplot(3,1,1);
plot(yCoordinates)
title('Your energy contour drawing');

subplot(3,1,2);
plot(newY, 'b--o')
title(['Your energy drawing for ' ,num2str(NumSongs), ' songs']);

subplot(3,1,3); 
PlaylistPower = [MixPlaylist.powerNormalized];
plot(PlaylistPower, 'r--o')
title('Created playlist energy distribution');


% subplot(4,1,3); 
% Test = [MixPlaylist.bpm];
% plot(Test)
% title('BPM');
% 
% subplot(4,1,4); 
% Test2 = [MixPlaylist.key];
% plot(Test2)
% title('Key');

end