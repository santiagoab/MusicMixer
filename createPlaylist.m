function playlistSongs = createPlaylist(sortedByEnergy,btnPressed)

switch btnPressed
    case 'Low'
        %Interval to select inferior 33 percentile
        X = prctile([sortedByEnergy.power],0);
        Y = prctile([sortedByEnergy.power],33);
        disp('Low energy songs')
    case 'Medium'
        %Interval to select middle 33 percentile
        X = prctile([sortedByEnergy.power],33);
        Y = prctile([sortedByEnergy.power],66);
        disp('Medium energy songs')
    case 'High'
        %Interval to select superior 33 percentile
        X = prctile([sortedByEnergy.power],66);
        Y = prctile([sortedByEnergy.power],100);
        disp('High energy songs')
    otherwise
        X = prctile([sortedByEnergy.power],0);
        Y = prctile([sortedByEnergy.power],100);
        disp('No selection made, all energy songs')
end

a = 1;

for i=1:length(sortedByEnergy) 
    if (X<=sortedByEnergy(i).power) && (sortedByEnergy(i).power<=Y)
        playlistSongs(a) = sortedByEnergy(i);
        a = a + 1;
    end
end

end