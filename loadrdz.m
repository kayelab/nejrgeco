function [raw, rz, long, long_z] = loadrdz(name, sframe) 
% Load .tif (via CaImAn-MATLAB toolbox), detrend, z-score
%
%   inputs:
%       name    :   full path to .tif file
%       sframe  :   first frame from which to start reading data
% 
%   outputs: 
%       raw     :   raw data from input .tif
%       rz      :   detrended, z-scored data from input .tif
%       long    :   output 'raw' in long (vectorized) format (i.e., space-by-time)
%       long_z  :   output 'rz' in long (vectorized) format (i.e., space-by-time)

% open file
raw = read_file(name, sframe);

% reshape red, greed, and roi data to vectorized space (x & y) by time
long = reshape(raw, [size(raw, 1) * size(raw, 2), size(raw, 3)]);

% detrend, zscore
long_z = zscore(detrend(long'))';

% reshape long_z back to 3D
rz = reshape(long_z, [size(raw, 1), size(raw, 2), size(raw, 3)]);
    
end