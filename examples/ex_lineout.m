clear all; 
close all; 
clc;

%% Load a single file
bv = BATSView('2d_idl.out','level')

%% Generate some lineouts of density
ysamp = linspace(0,100,2);
figure(1);
cla; hold on;
for i=1:numel(ysamp)
    i
    start  = [   0,ysamp(i)];
    finish = [2000,ysamp(i)];
    [coords, values] = bv.lineout(1,'level',start,finish,1000);
    plot(coords(:,1),values,'k-');
end



%% Generate a pseudocolor plot of density
% figure(2);
% cla; hold on;
% colormap('visit_hot_desaturated')
% bndbox = [[300,1000];[0,700];[0,0]]
% minsizes = [0.3,0.3,0]
% [X,Y,Z,D] = bv.get_uniform(1,'rho',bndbox,minsizes);
% h=pcolor(X,Y,D);
% % set(h,'edgealpha',0);
% shading interp

