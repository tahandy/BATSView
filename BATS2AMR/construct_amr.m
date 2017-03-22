clear all; 
close all; 
clc;

%% Load a single file
% bv = BATSView('2d_idl.out','level');


%%
clc;

% x = linspace(intmin('uint32')/2,intmax('uint32')/2,2);
% y = linspace(intmin('uint32')/2,intmax('uint32')/2,2);
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



mortonOrder = zeros(size(X));

for i=1:1
    ix = typecast(cast(X(i),'single'),'uint32')
    iy = typecast(cast(Y(i),'single'),'uint32')
    iz = typecast(cast(0.0,'single'),'uint32')

    a_ui32 = ix;
    x_ui64 = uint64(ix)
    x_ui64 = bitand(x_ui64,typecast(hex2dec('1fffff'),'uint64'))

    x_ui64 = bitand(bitor(x_ui64,bitsll(x_ui64,32)),hex2dec('1f00000000ffff'))
    x_ui64 = bitand(bitor(x_ui64,bitsll(x_ui64,16)),hex2dec('1f0000ff0000ff'))
    x_ui64 = bitand(bitor(x_ui64,bitsll(x_ui64, 8)),hex2dec('100f00f00f00f00f'))
    x_ui64 = bitand(bitor(x_ui64,bitsll(x_ui64, 4)),hex2dec('10c30c30c30c30c3'))
    x_ui64 = bitand(bitor(x_ui64,bitsll(x_ui64, 2)),hex2dec('1249249249249249'))


    a_ui32 = iy;
    y_ui64 = uint64(a_ui32);
    y_ui64 = bitand(y_ui64,hex2dec('1fffff'));

    y_ui64 = bitand(bitor(y_ui64,bitsll(y_ui64,32)),hex2dec('1f00000000ffff'));
    y_ui64 = bitand(bitor(y_ui64,bitsll(y_ui64,16)),hex2dec('1f0000ff0000ff'));
    y_ui64 = bitand(bitor(y_ui64,bitsll(y_ui64, 8)),hex2dec('100f00f00f00f00f'));
    y_ui64 = bitand(bitor(y_ui64,bitsll(y_ui64, 4)),hex2dec('10c30c30c30c30c3'));
    y_ui64 = bitand(bitor(y_ui64,bitsll(y_ui64, 2)),hex2dec('1249249249249249'));


    a_ui32 = iz;
    z_ui64 = uint64(a_ui32);
    z_ui64 = bitand(z_ui64,hex2dec('1fffff'));

    z_ui64 = bitand(bitor(z_ui64,bitsll(z_ui64,32)),hex2dec('1f00000000ffff'));
    z_ui64 = bitand(bitor(z_ui64,bitsll(z_ui64,16)),hex2dec('1f0000ff0000ff'));
    z_ui64 = bitand(bitor(z_ui64,bitsll(z_ui64, 8)),hex2dec('100f00f00f00f00f'));
    z_ui64 = bitand(bitor(z_ui64,bitsll(z_ui64, 4)),hex2dec('10c30c30c30c30c3'));
    z_ui64 = bitand(bitor(z_ui64,bitsll(z_ui64, 2)),hex2dec('1249249249249249'));


    answer = uint64(0);
    answer = bitor(answer,x_ui64);
%     answer = bitor(answer,bitsll(y_ui64,1));
%     answer = bitor(answer,bitsll(z_ui64,2));
    
    mortonOrder(i) = answer;
end

%%
IND = 1:numel(mortonOrder);
[tmp,IND] = sort(mortonOrder);

%%
figure(1);
cla; hold on;

cmap = jet(50);

plot(X(IND),Y(IND),'ko');

for i=1:numel(IND)-1
    plot([X(i),X(i+1)],[Y(i),Y(i+1)],'color',cmap(i,:));
end

