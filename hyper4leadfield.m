function [LeadField, coordinateTransformParameters] = hyper4leadfield(EEGarray)
% hyper4leadfield() % Calculates the leadfield using REST
% INPUT
% EEGarray                          EEGLAB array including channels locations
% OUTPUT
% LeadField                         leadfield for the electrode array
% coordinateTransformParameters     Coordinates of the coregister of the
%                                   electrodes location to the
%                                   'standard_1005'
%                                   e.g. coordinateTransformParameters = [ 0.3100  -19.2162  -43.1129    0.0823    0.0030   -1.5751  109.3809  102.7028  138.8204 ];

% Electrodes coordinates (cartesian)
Xelec = arrayfun(@(x) (x.X), EEGarray.chanlocs);
Yelec = arrayfun(@(x) (x.Y), EEGarray.chanlocs);
Zelec = arrayfun(@(x) (x.Z), EEGarray.chanlocs);
xyz_elec = cat(1, Xelec, Yelec, Zelec)';

% Load fixed dipoles coordinates
[ProgramPath, ~, ~] = fileparts(which('pop_REST_reref.m'));
xyz_dipoles = load([ProgramPath,filesep,'corti869-3000dipoles.dat']);

% Calculate the dipole orientations.
xyz_dipOri = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
xyz_dipOri ( 2601: 3000, 1 ) = 0;
xyz_dipOri ( 2601: 3000, 2 ) = 0;
xyz_dipOri ( 2601: 3000, 3 ) = 1;

% Define headmodel
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ];
headmodel.r      = [ 0.8700,0.9200,1];
headmodel.cond   = [ 1.0000,0.0125,1];
headmodel.tissue = { 'brain' 'skull' 'scalp' };

% Calculate leadfield
[G,~] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel);
LeadField = G';

[~,coordinateTransformParameters] = coregister(EEGarray.chanlocs, ...
    [fileparts(which('standard_1005.elc')) filesep 'standard_1005.elc'], ...
    'warp', 'auto', 'manual', 'off');
end
