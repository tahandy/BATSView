clear all; close all; clc;

fmtdec = @(str) fprintf('%s\n',[str(1:8),' ',str(9:16),' ',str(17:24),' ',str(25:32)]);

nx    = 4;
ny    = 4;
xbnds = single([-1,1]);
ybnds = single([-pi,exp(1)]);
x = linspace(xbnds(1),xbnds(2),nx);
y = linspace(ybnds(1),ybnds(2),ny);
[X,Y] = ndgrid(x,y);

X = (X-xbnds(1))/diff(xbnds);
Y = (Y-ybnds(1))/diff(ybnds);

X = X(:);
Y = Y(:);

% [X,Y] = ndgrid(uint32(0:nx-1),uint32(0:ny-1))


% la = X<0.5 & Y<0.5;
% X(la) = [];
% Y(la) = [];
% x = linspace(0,0.5,nx+1); x = x(1:nx);
% y = linspace(0,0.5,ny+1); y = y(1:ny);
% [Xf,Yf] = ndgrid(x,y);
% X(end+1:end+nx*ny) = Xf(:);
% Y(end+1:end+nx*ny) = Yf(:);

order = uint32(size(X));
for i=1:numel(X)
%     ix = typecast(cast(X(i),'single'),'uint32')
%     iy = typecast(cast(Y(i),'single'),'uint32')
    
    bx = single(str2num(Fr_dec2bin(X(i))));
    by = single(str2num(Fr_dec2bin(Y(i))));
    
    ix = typecast(bx,'uint32')
    iy = typecast(by,'uint32')

    
    
    order(i) = morton2(ix,iy);
end

%%
% fmtdec(dec2bin(order(13),32))
% fmtdec(dec2bin(order(14),32))
% fmtdec(dec2bin(order(15),32))
% fmtdec(dec2bin(order(16),32))
% fmtdec(dec2bin(order(17),32))
 
[~,I] = sort(order);
plot(X(I),Y(I),'k-o')