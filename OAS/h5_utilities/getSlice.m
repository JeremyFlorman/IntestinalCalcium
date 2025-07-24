function [img] = getSlice(slice, indices)
%GETSLICE Summary of this function goes here
%   Detailed explanation goes here

file_index = indices.h5_file_index(slice);
path = indices.h5_path{file_index};
img_size = indices.img_size{file_index};




img = h5read(path, '/data',[1 1 indices.relative_index(slice)],[img_size(1) img_size(2) 1]);

end

