function fixMicroManagerFileNames(folder)
%fixMicroManagerFileNames Removes bug-causing file suffix ('_MMStack_Default.ome')
% from micro-manager generated .tif files
%   will remove all '_MMStack_Default.ome' frome file names found in
%   recursive search of the input variable 'folder'.

d = dir([folder '\**\*.tif']);

for i =1:length(d)
    fn = fullfile(d(i).folder, d(i).name);
    newname = strrep(fn, '_MMStack_Default.ome', '')
    movefile(fn, newname)
end

end