% compute the radio interference using normal implementation. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function rif = ristdexc( method, linecenter, nsubline, dsubline, Y, Z, emax, rouearth, xpos, riheight)
if ( min(dsubline) < 0.5) dsubline = dsubline .* 100; end
n = length(xpos);
rif = zeros(n,1);
ph = length(emax);
mu = 4*pi*1e-7; 
freq = 500e3; 
episilon = 8.85e-12; 

x = linecenter(:,1);
y = linecenter(:,2);

x1 = repmat( x', n, 1);
y1 = repmat( y', n, 1);

xt = repmat( xpos, 1, ph);
yt = repmat( riheight * ones(n, 1), 1, ph);

dmirror = sqrt( rouearth / (pi * freq * mu));

x1t = (x1-xt).^2;
dis = (y1-yt)./( (y1-yt).^2 + x1t) + (y1+yt+2*dmirror) ./ ( x1t + (y1+yt+2*dmirror).^2);

C = imag(Y)./(2*pi*freq);

if ( method == 1)  % 1 guo biao
    tao = 10.^((70 - 585./ emax + 35.*log10( dsubline)-10.*log10(nsubline))/20) ;
else
    tao = 10.^((37.02 + 120.*log10( emax./15) + 40*log10(dsubline./4))/20);
end
Tao = diag(tao); 
               
ZY = Z*Y;
[S,Lm] = eig(ZY);  % eig in matlab: s * lm * inv(s) = ZY. this is different from the handbook of 345 kv.
alpha = real( sqrt(Lm));

diaga = diag(alpha);

invA = C ./ (2*pi*episilon);

R = S * Tao;  % defination of eig is different.

QR = Q .* R;

QW = 30 .* invA .* inv(S) * R;
% should be 20 m outside and 2 m high.

% alpham = repmat( diaga, 1, ph) + repmat( diaga.', ph, 1);    
m1 = repmat( diaga, 1, ph);
alpham = m1+m1.';

s = ph * ph;
v = 1:s;
v = reshape(v, ph, ph)';

Q = 30 .* dis * M;

for i = 1:n
    f = dis(i,:);
    Q = 30 * f * M;
    W = diag(Q)*R;
    % ei = sigma( wmi.*exp...), 
    tp=0; tmp=0;
    wt = reshape( W, [], 1);
    wt1 = repmat( wt, 1, ph); 
    wt2 = repmat( wt', ph, 1);
    for k = 1:ph       % for each phase.
        w1 = repmat( W(:,k),1, ph);
        w2 = repmat( W(:,k).', ph, 1);
        t = w1 .* w2 ./ alpham;        
        tmp = sum(sum(t)) + tmp;    
        tp = sum(t(1:s))+tp;
    end   
    
    rif(i) = 10*log10((abs(tp*2)));
end