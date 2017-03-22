clear all; 
close all; 
clc;

%% Load a single file
% bv = BATSView('2d_idl.out','level');


%%
clc;

x = 0:3;
y = 0:3;

[X,Y] = ndgrid(x,y);
X = X(:);
Y = Y(:);


x = 0:10
y = 0:3


[X,Y] = ndgrid(x,y);
X = X(:);
Y = Y(:);


dx0 = 0.5;
dx1 = 0.25;

x1 = linspace(0,1,8+1); x1 = 0.5*(x1(1:end-1)+x1(2:end));
x2 = linspace(0,2,8+1); x2 = 0.5*(x2(1:end-1)+x2(2:end));

dx = x1(2)-x1(1);

[X1,Y1] = ndgrid(x1,x1);
[X2,Y2] = ndgrid(x2,x2);

X2 = X2(:);
Y2 = Y2(:);

la = X2<1 & Y2<1;
X2(la) = [];
Y2(la) = [];

X = [X2(:);X1(:)];
Y = [Y2(:);Y1(:)];


X = round(X(:)/dx);
Y = round(Y(:)/dx);




mortonOrder = uint32(size(X));

% for i=2:2
for i=1:numel(X)
    ix = typecast(cast(X(i),'single'),'uint32');
    iy = typecast(cast(Y(i),'single'),'uint32');


    x_ui32 = ix;
    x_ui32 = bitand(x_ui32,typecast(cast(hex2dec('0000ffff'),'single'),'uint32'));
    dec2bin(x_ui32)
    x_ui32 = bitand(typecast(cast(hex2dec('00ff00ff'),'single'),'uint32'),bitxor(x_ui32,bitsll(x_ui32,8)));
    x_ui32 = bitand(typecast(cast(hex2dec('0f0f0f0f'),'single'),'uint32'),bitxor(x_ui32,bitsll(x_ui32,4)));
    x_ui32 = bitand(typecast(cast(hex2dec('33333333'),'single'),'uint32'),bitxor(x_ui32,bitsll(x_ui32,2)));
    x_ui32 = bitand(typecast(cast(hex2dec('55555555'),'single'),'uint32'),bitxor(x_ui32,bitsll(x_ui32,1)));


    y_ui32 = iy;
    y_ui32 = bitand(y_ui32,typecast(cast(hex2dec('0000ffff'),'single'),'uint32'));
    y_ui32 = bitand(typecast(cast(hex2dec('00ff00ff'),'single'),'uint32'),bitxor(y_ui32,bitsll(y_ui32,8)));
    y_ui32 = bitand(typecast(cast(hex2dec('0f0f0f0f'),'single'),'uint32'),bitxor(y_ui32,bitsll(y_ui32,4)));
    y_ui32 = bitand(typecast(cast(hex2dec('33333333'),'single'),'uint32'),bitxor(y_ui32,bitsll(y_ui32,2)));
    y_ui32 = bitand(typecast(cast(hex2dec('55555555'),'single'),'uint32'),bitxor(y_ui32,bitsll(y_ui32,1)));

    zvalue = x_ui32 + bitsll(y_ui32,1);
%     zvalue = bitsll(x_ui32,-1) + y_ui32;
    
    mortonOrder(i) = zvalue;
end

%%
[tmp,IND] = sort(mortonOrder);
IND = flipud(IND);

%%
figure(1);
cla; hold on;

cmap = jet(100);

plot(X(IND),Y(IND),'ko');

for i=1:numel(IND)-1
    plot([X(i),X(i+1)],[Y(i),Y(i+1)],'color',cmap(i,:));
end


    