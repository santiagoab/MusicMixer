function drawnPlaylist = drawProcessing(xCoordinates,yCoordinates, features, NumSongs)

%Resample the yCoordinates from the hand drawing to the same length of the
%number of desired songs, thus each value of the yCoordinates corresponding to the
% closest normalised energy of a song.
newY = resample(yCoordinates,NumSongs,length(yCoordinates));



for i=1:length(newY)
    
    [~, index] = min(abs([features.powerNormalized]-newY(i)));  %return index of closest song by energy to i drawing point
    
    drawnPlaylist(i) = features(index);%set that song on the created playlist
    
    features(index).powerNormalized = NaN; %set that song power to empty to avoid using it twice
    
    %Otras opciones? Interpolation? N-neighbours? ?¿?¿?¿
   
end

subplot(2,1,1);
plot(newY)
title(['Your energy drawing for ' ,num2str(NumSongs), ' songs.']);

subplot(2,1,2); 
PlaylistPower = [drawnPlaylist.powerNormalized];
plot(PlaylistPower)
title('Created playlist energy distribution');

end