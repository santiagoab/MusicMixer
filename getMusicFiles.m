function [listing,folder] = getMusicFiles()
selpath = uigetdir;
%contents(3).name;
folder = [selpath filesep]; %Add needed \ to path
%listing = dir(selpath);
listing = dir([selpath, '/*.mp3']) %just look for .mp3 to avoid loading .mat

inds = [];
n    = 0;
k    = 1;

while n < 2 && k <= length(listing)
    if any(strcmp(listing(k).name, {'.', '..'}))
        inds(end + 1) = k;
        n = n + 1;
    end
    k = k + 1;
end

listing(inds) = [];

end