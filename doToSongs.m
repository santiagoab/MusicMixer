function features = doToSongs(folder,songs)

for i=1:length(songs)
    disp([folder songs(i).name]);
    features(i) = computeFeatures(folder,songs(i).name);      
end

disp('Songs features extracted successfully!');

end
  