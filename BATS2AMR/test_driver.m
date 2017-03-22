clear all; close all; clc;

domainBB = [[0,2];[0,1]];
nRoot = [1,1];
NB = 8;

x1 = linspace(0,1,NB+1); x1 = 0.5*(x1(1:end-1)+x1(2:end));
x2 = linspace(0,1,NB+1); x2 = 0.5*(x2(1:end-1)+x2(2:end));

[X1,Y1] = ndgrid([x1,x1+1],[x1]);

X = [X1(:)];
Y = [Y1(:)];

dxmin = get_min_spacing(X)
dymin = get_min_spacing(Y)

nzoneX = (domainBB(1,2)-domainBB(1,1))/dxmin
nzoneY = (domainBB(2,2)-domainBB(2,1))/dymin


nrootX = nzoneX/gcd(nzoneX,nzoneY)
nrootY = nzoneY/gcd(nzoneX,nzoneY)


return

%%
figure(1);
cla; hold on
plot(X,Y,'o')
