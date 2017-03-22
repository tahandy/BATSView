function n_ui32 = Part1By1(a_ui32)
fmtdec = @(str) fprintf('%s\n',[str(1:8),' ',str(9:16),' ',str(17:24),' ',str(25:32)]);

n = a_ui32;

% fmtdec(dec2bin(n,32));

n = bitand(bitxor(n, bitsll(n,8)),hex2dec('00ff00ff'));
% fmtdec(dec2bin(n,32));
n = bitand(bitxor(n, bitsll(n,4)),hex2dec('0f0f0f0f'));
% fmtdec(dec2bin(n,32));
n = bitand(bitxor(n, bitsll(n,2)),hex2dec('33333333'));
% fmtdec(dec2bin(n,32));
n = bitand(bitxor(n, bitsll(n,1)),hex2dec('55555555'));
% fmtdec(dec2bin(n,32));

n_ui32 = n;