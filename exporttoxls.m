% Export the results to .xls files.
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function exporttoxls( xpos, E, B, ri1, ri2, anbpa, an345, emaxcon, emaxg,...
                    coronainit, ginit, einit_west, eginit_west, aveloss, peakloss,...
                    natpow, waveimp)
% 
xlsfile = 'ResultAC.xls';
xlssheet = 'Result';
lenxpos = length(xpos);
lencell = lenxpos + 2; % ������
widcell = 'Q'-'A' + 'Z' - 'A' + 1 + 1;
data = cell(lencell, widcell);

m1 = size(xpos,2); 
if ( m1 > 1) 
    xpos = xpos'; 
end
m1 = size(E,2); 
if ( m1 > 1) 
    E = E'; 
end
m1 = size(B,2); 
if ( m1 > 1) 
    B = B'; 
end
m1 = size(ri1,2); 
if ( m1 > 1) 
    ri1 = ri1'; 
end
m1 = size(ri2,2); 
if ( m1 > 1) 
    ri2 = ri2'; 
end
m1 = size(anbpa,2); 
if ( m1 > 1) 
    anbpa = anbpa'; 
end
m1 = size(an345,2); 
if ( m1 > 1) 
    an345 = an345'; 
end
pack1 = [xpos, E, B, ri1, ri2, anbpa, an345];
pack1col = 'K'-'A' + 1 + 'Z' - 'A' + 1;
data(3:lencell, pack1col:(pack1col+7-1)) = num2cell(pack1);

lenemaxcon = length(emaxcon);
row = 19:21;
col = 'W'-'A'+1;
pack2 = [emaxcon,  coronainit, einit_west]';
data(row, col:col+lenemaxcon-1) = num2cell(pack2);

row = 23;
col = 'Z'-'A'+1;
k = length(emaxg);
pack3 = [1:k; emaxg'];
data( row:row+1, col:(col+k-1)) = num2cell(pack3); % �糡ǿ�ȵ���

row = 26;
data(row, col) = num2cell( max(emaxg));  % ���ߵ��ǿ
pack4 = [ginit eginit_west]';
data( row+1:row+2, col:(col+k-1)) = num2cell(pack4);

row = 30;
col = 'W'-'A'+1;
k = length(aveloss);
pack5 = [aveloss, peakloss]';
data( row:row+1, col:(col+k-1)) = num2cell(pack5);

row = 32;
k = length(natpow);
pack6 = [natpow, waveimp]';
data( row:row+1, col:3:(col+3*k-3+1)) = num2cell(pack6);

xlswrite( xlsfile,data, xlssheet,'A1');