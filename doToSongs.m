function features = doToSongs(folder,songs)

for i=1:length(songs)
    disp([folder songs(i).name]);
    features(i) = computeFeatures(folder,songs(i).name); 
end

disp('Song features extracted succesfully');

%Normalize features
bpmNormalized = ([features.bpm] - min([features.bpm])) / ( max([features.bpm]) - min([features.bpm]) );
powerNormalized = ([features.power] - min([features.power])) / ( max([features.power]) - min([features.power]) );
keyNormalized = ([features.key] - min([features.key])) / ( max([features.key]) - min([features.key]) );
modeNormalized = ([features.mode] - min([features.mode])) / ( max([features.mode]) - min([features.mode]) );

%save into features struct
bpmNormalized = num2cell(bpmNormalized');
powerNormalized = num2cell(powerNormalized');
keyNormalized = num2cell(keyNormalized');
modeNormalized = num2cell(modeNormalized');

[features(:).bpmNormalized] = deal(bpmNormalized{:});
[features(:).powerNormalized] = deal(powerNormalized{:});
[features(:).keyNormalized] = deal(keyNormalized{:});
[features(:).modeNormalized] = deal(modeNormalized{:});

for i=1:length(songs)
    scatter([features(i).bpm],[features(i).powerNormalized]);
    text([features(i).bpm],[features(i).powerNormalized], features(i).name, 'horizontal','left', 'vertical','bottom','FontSize',10);
    xlabel('BPM')
    ylabel('Energy Power')
    hold on
end


end

