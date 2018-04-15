function outputs = doToSongs(folder,songs)

for i=1:length(songs)
    disp([folder songs(i).name]);
    features = computeFeatures([folder songs(i).name]);
    
    
end
end
  