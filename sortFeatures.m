function sortedByBPM = sortFeatures(features)

Afields = fieldnames(features);
Acell = struct2cell(features);
sz = size(Acell);
Acell = reshape(Acell, sz(1), []);
Acell = Acell';  
Acell = sortrows(Acell, 3);
Acell = reshape(Acell', sz);
sortedByBPM = cell2struct(Acell, Afields, 1);

end
  