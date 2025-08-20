function [remotedata] = copyLoadClear(file2load, tempdir)
%copyLoadClear Copies a remote file to a local directory for loading, then
%deletes local file.
%   Detailed explanation goes here

if iscell(file2load)
    file2load = file2load{:};
end

[~, name, ext] = fileparts(file2load);
copyfile(file2load, tempdir)

remotedata = load(fullfile(tempdir, name), '-mat');
localfilepath = [tempdir '\' name ext];
delete(localfilepath)

end

