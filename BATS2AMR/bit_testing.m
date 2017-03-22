clear all;
close all;
clc;

fmtdec = @(str) fprintf('%s\n',[str(1:8),' ',str(9:16),' ',str(17:24),' ',str(25:32)]);

a_ui32 = uint32(6);

fmtdec(dec2bin(a_ui32,32));

fmtdec(dec2bin(bitsll(a_ui32,1),32));
fmtdec(dec2bin(bitshift(a_ui32,1),32));



a = single(0.011)
a_ui32 = typecast(a,'uint32');

fprintf('a bits:    ');
fmtdec(dec2bin(a_ui32,32));


signmask = false(1,32); 
signmask(1) = true;

expmask = false(1,32);
expmask(2:9) = true;

fracmask = false(1,32);
fracmask(10:end) = true;


fprintf('signmask:  ');
fmtdec(dec2bin(binaryVectorToDecimal(signmask),32))

fprintf('expmask:   ');
fmtdec(dec2bin(binaryVectorToDecimal(expmask),32))

fprintf('fracmask:  ');
fmtdec(dec2bin(binaryVectorToDecimal(fracmask),32))


s = bitand(a_ui32,binaryVectorToDecimal(signmask))
s = bitsra(s,31)

e = bitand(a_ui32,binaryVectorToDecimal(expmask))
fmtdec(dec2bin(e,32));
e = bitsra(e,23)
fmtdec(dec2bin(e,32));
e = int32(e)-127

f = bitand(a_ui32,binaryVectorToDecimal(fracmask))
fbits = fliplr(bitget(f,1:23));
fmtdec(dec2bin(f,32));

f = 1 + sum(single(2.^-(1:numel(fbits))).*single(fbits))

f*2.0^single(e)







