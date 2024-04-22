function ridb = ristd( emax, linepos, rsubline, xt, ht)
n = length(xt);
m = length( emax);   
x1 = repmat( linepos(:,1), 1, n);
y1 = repmat( linepos(:,2), 1, n);
x2 = repmat( xt.', m, 1);
y2 = ones(m,n) .* ht;
d = (x1-x2).^2 + (y1-y2).^2;
ridbtmp = 3.5 .* emax + 12.* rsubline - 30 .* log10(20./d);  
pa = 1:3:(m-2);  ria = ridbtmp(pa,:);  % A phase
pb = 2:3:(m-1);  rib = ridbtmp(pb,:);  % B phase
pc = 3:3:m;      ric = ridbtmp(pc,:);  % C phase
if ( m > 3) 
    ria = (10.^( ria ./20)).^2;   ria = sum(ria, 1);
    rib = (10.^( rib ./20)).^2;   rib = sum(rib, 1);
    ric = (10.^( ric ./20)).^2;   ric = sum(ric, 1);
    ridbtmp = 10 .* log10([ria; rib; ric]);
end
rmax = max( ridbtmp);
rmin = max( ridbtmp);
rmed = ria + rib + ric - rmax - rmin;
kt = (rmax - rmin > 3);
ridb(kt) = rmax(kt);
kt = ~kt;
ridb(kt) = (rmax(kt) + rmed(kt)) * 0.5 + 1.5;