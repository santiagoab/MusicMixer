function features = doToSongs(folder,songs)

for i=1:length(songs)
    disp([folder songs(i).name]);
    features(i) = computeFeatures(folder,songs(i).name); 
    scatter([features(i).bpm],[features(i).power]);
    text([features(i).bpm],[features(i).power], features(i).name, 'horizontal','left', 'vertical','bottom','FontSize',8);
    xlabel('BPM')
    ylabel('Energy Power')
    hold on
end

bpmNormalized = [features.bpm]/norm([features.bpm]);
powerNormalized = [features.power]/norm([features.power]);

end

