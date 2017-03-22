clear all; 
close all; 
clc;

%% Load a single file
% To load multiple files, the first argument (filenames) should
% be a cell array of strings

bv = BATSView('2d_idl.out',{'level','rho})


% Specify the bounding box for the portion of the domain
% which will be uniformily refined
% 
$ This is specified as a (3x2) array, containing the values
% [x_min, x_max;
%  y_min, y_max;
%  z_min, z_max]
%
% Only rows up to nDim are considered, with the remainder ignored.
% Sampling points are located on the edges of this bounding box, 
% i.e. face/edge centered, and not inset (i.e. cell centered).
bndbox = [300, 1000;
            0,  700;
            0,    0];


% Specify the spacing between points in the uniformily defined grid.
%
% This is given by a vector of length 3, considering elements only up to 
% nDim.

minsizes = [1, 1, 0];


% Obtain the ndgrid-formed coordinates and resampled value of density on the
% uniform grid. The first argument (1) denotes to use the first file
% read into the BATSView object. For multi-file BATSView objects,
% this value can be iterated over up to the number of read files.
%
% By default, interpolation is performed by using the 'nearest' flag, 
% choosing the value at the data point closest to the requested point.
% This may be changed by modifying get_uniform and using 'linear'.

[Xu,Yu,Zu,Vu] = bv.get_uniform(1,'rho',bndbox,minsizes)


%% Generate a pseudocolor plot of density
figure(1);
cla; hold on;
%colormap('visit_hot_desaturated')
h=pcolor(X,Y,D);
set(h,'edgealpha',0);
shading interp

