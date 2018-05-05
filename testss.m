s = [1 1 1 2 2 1 6 2];
t = [6 3 4 5 6 6 2 7];
w = [4 56 57 87 23 500 1 2];
figure;
G = digraph(s,t,Dists);
plot(G,'EdgeLabel',G.Edges.Weight)

function drawnPlaylist = dataProcessing(yCoordinates, features, NumSongs)

%Resample the yCoordinates from the hand drawing to the same length of the
%number of desired songs, thus each value of the yCoordinates corresponding to the
% closest normalised energy of a song.
newY = resample(yCoordinates,NumSongs,length(yCoordinates));

% energy range for songs in a position. Range usually 0.1
range = 0.1;

%weights for song parameters
wbpm = 1;
wkey = 1;
wmode = 1; 

%vectors init
s = [];
t = [];
Dists =[];


%for each position in the drawing select closest songs by energy to the drawn curve within an energy range
for i=1:length(newY)
    
    [~,position{i}]=find(abs([features.powerNormalized]-newY(i))<range);  % range usually 0.1
    
    
    % if there are no songs within energy range, get closest song
    if isempty(position{i})
       [~, position{i}] = min(abs([features.powerNormalized]-newY(i)));   
    end
    
    %[~, index] = min(abs([features.powerNormalized]-newY(i)));  %return
    %index of closest song by energy to i drawing point
%     drawnPlaylist(i) = features(index);%set that song on the created playlist
%     features(index).powerNormalized = NaN; %set that song power to empty to avoid using it twice 
end


%Weights between current and next position possible songs
for j=1:length(position - 1)
currentPos = position{j};
nextPos = position{j+1};

    for k=1:length(currentPos)
        %get current pos attributes for weights (BPM, key, mode)
        currentSong = currentPos(k);
        currentSongBpm = features(currentSong).bpm;
        currentSongKey = features(currentSong).key;
        currentSongMode = features(currentSong).mode;
        
        for y = 1:length(nextPos)
           %get next pos attributes 
            nextSong = nextPos(y);
            nextSongBpm = features(y).bpm;
            nextSongKey = features(y).key;
            nextSongMode = features(y).mode;
            
           %compute weighted distances between such songs 
           dist = abs(currentSongBpm - nextSongBpm) * wbpm;  %RETOCAT
              
           %store S, T and Dists
            s = cat(2,s,currentSong);
            t = cat(2,t,nextSong);
            Dists = cat(2,Dists,dist);
            
        end
          
    end
    
   %find shortest path 
    
   % if repeating songs, remove and recompute
   
   % create final mix drawnPlaylist(i) or mixPlaylist
    
  
end


subplot(2,1,1);
plot(newY)
title(['Your energy drawing for ' ,num2str(NumSongs), ' songs.']);

subplot(2,1,2); 
PlaylistPower = [drawnPlaylist.powerNormalized];
plot(PlaylistPower)
title('Created playlist energy distribution');

end