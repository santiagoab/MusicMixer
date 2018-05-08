%{
[songs, folder] = getMusicFiles();
features = doToSongs(folder,songs);
%}
function features = doToSongs(folder,songs)

for i=1:length(songs)
    %disp([folder songs(i).name]);
    features(i) = computeFeatures(folder,songs(i).name); 
    scatter([features(i).bpm],[features(i).power]);
    text([features(i).bpm],[features(i).power], features(i).name, 'horizontal','left', 'vertical','bottom','FontSize',8);
    xlabel('BPM')
    ylabel('Energy Power')
    hold on
end

disp('Song features extracted succesfully');

%Normalize features
bpmNormalized = ([features.bpm] - min([features.bpm])) / ( max([features.bpm]) - min([features.bpm]) );
powerNormalized = ([features.power] - min([features.power])) / ( max([features.power]) - min([features.power]) );

%save into features struct
bpmNormalized = num2cell(bpmNormalized');
powerNormalized = num2cell(powerNormalized');

[features(:).bpmNormalized] = deal(bpmNormalized{:});
[features(:).powerNormalized] = deal(powerNormalized{:});

end

