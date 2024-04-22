% compute the power line parameters (which are frequency dependent) at Frequence
% The results are accurate in the sense that it matches with PACAD tline.exe
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% This code was initially written by Gong Youjun in 2007, and optimized by chijie@tsinghua.edu.cn around 2009.

function [Y,Z]= YZ( CurrentLinePos, rCurrentSubLine,rhoCurrentLine,EarthLinePos,rEarthLine,rhoEarthLine,Frequence,rhoG,totalroute,phase_info)
episilon = 1/36/pi*1e-9;
if ( rhoG <= 0) rhoG = 100; end

n1 = size( CurrentLinePos,1);
n2 = size(EarthLinePos,1);
n = n1 + n2;
phase = 3*totalroute;

rM = [rEarthLine;rCurrentSubLine];

LinePos = [EarthLinePos;CurrentLinePos];
P1 = zeros(n,n);
for h = 1:n
    for k = 1:n
        if ( h ==k)
            P1(h,k) = log(2*LinePos(h,2)/rM(h));
        else
            DD(h,k) = sqrt((LinePos(h,1)-LinePos(k,1))^2+(LinePos(h,2)+LinePos(k,2))^2);
            dd(h,k) = sqrt((LinePos(h,1)-LinePos(k,1))^2+(LinePos(h,2)-LinePos(k,2))^2);
            P1(h,k) = log(DD(h,k)/dd(h,k));
        end
    end
end

Pll = P1(1:n2, 1:n2);
Plt = P1(1:n2, (n2+1):n);
Ptl = P1((n2+1):n,1:n2);
Ptt = P1((n2+1):n,(n2+1):n);
Pe = Ptt - Ptl*Pll\Plt; % Pe = Ptt - Ptl*inv(Pll)*Plt;
P1 = Pe;

for i = 1:totalroute
    phaseinfo((i-1) * 3 + 1) = phase_info(i,1);
    phaseinfo((i-1) * 3 + 2) = phase_info(i,2);
    phaseinfo((i-1) * 3 + 3) = phase_info(i,3);
end
phaseinfo( totalroute * 3 + 1) = phase_info( totalroute+1,1);


for k = 1:phase
    Noi = phaseinfo(k);
    for i = (Noi+1): phaseinfo(k+1)-1
        P1(:,i) = P1(:,i) - P1(:,Noi);
    end
    for i = (Noi+1): phaseinfo(k+1)-1
        P1(i,:) = P1(i,:) - P1(Noi,:);
    end
end

for k = 1:phase
    Noi = phaseinfo(k); % 
    Pend(k,:) = P1(Noi,:);
    a1 = phase + phaseinfo(k)-(k-1);
    b1 = a1 + phaseinfo(k+1)-phaseinfo(k) - 2;
    a2 = Noi + 1;
    b2 = a2 + phaseinfo(k+1)-phaseinfo(k) - 2;
    Pend(a1:b1, :) = P1(a2:b2, :);  
end

for k = 1:phase
    Noi = phaseinfo(k); % 
    Pend_(:,k) = Pend(:,Noi);
    a1 = phase + phaseinfo(k)-(k-1);
    b1 = a1 + phaseinfo(k+1)-phaseinfo(k) - 2;
    a2 = Noi + 1;
    b2 = a2 + phaseinfo(k+1)-phaseinfo(k) - 2;
    Pend_( :,a1:b1) = Pend(:,a2:b2);  
end

Pell = Pend_(1:phase,1:phase);
Pelt = Pend_(1:phase, phase+1:size(Pend,1));
Petl = Pend_(phase+1:size(Pend,1),1:phase);
Pett = Pend_(phase+1:size(Pend,1),phase+1:size(Pend,1));
Pee =  Pell - Pelt* Pett\ Petl; % Pee =  Pell - Pelt*inv(Pett)*Petl;
C = inv(Pee).*2*pi*episilon;
Y = sqrt(-1)*2*pi*Frequence*C;

rM(1:n2,1) = rEarthLine;
rM = [rM; rCurrentSubLine];

rhoEarthLine = rhoEarthLine/1000;
rhoCurrentLine = rhoCurrentLine /1000;

dao_AL = 3.3536*10^7;
dao_Fe = 1.0256*10^7;
mui = 4*pi*10^(-7);

L = zeros(n,n);
j=sqrt(-1);
P = 1/sqrt(j*Frequence*2*pi*4*pi*10^(-7)/rhoG);
Zc = L; Z = L;
for h = 1:n
    for k = 1:n
        if ( h ==k)
            if ( h <= n2 &&  k <= n2)
%                 radius = rhoEarthLine(h);                
%             else
%                 radius = rhoCurrentLine(h-n2);
%             end
                m = sqrt(2*pi*Frequence*1.00022*mui*dao_Fe);
                delta = 1/sqrt(pi*Frequence*1.00022*mui*dao_Fe);
                r = rM(h);
                q = r-(r/delta-0.001)*delta;
                             
                k1 = (besseli(0,sqrt(j)*(m*r+0.001))-besseli(0,sqrt(j)*m*r))/0.001;
                mmc = sqrt(j)*(m*q);
                k1 = difbesseli0j(mmc);
                
                k2 = (besselk(0,sqrt(j)*(m*r+0.001))-besselk(0,sqrt(j)*m*r))/0.001;
 %               k2 = difbesselk0j(mmc);
                
                fai = -(besseli(0,sqrt(j)*(m*q+0.001))-besseli(0,sqrt(j)*m*q))/(besselk(0,sqrt(j)*(m*q+0.001))-besselk(0,sqrt(j)*m*q));
%                 fai = -difbesseli0j(sqrt(j)*(m*q)) / difbesselk0j(0,sqrt(j)*m*q);
%                 fai = -k1/k2;
                
                Zc(h,k) = rhoEarthLine(h) *j*0.5*m*r*(+fai*besselk(0,sqrt(j)*m*r))/(k1+fai*k2)*(1-(q/r)^2);
               
                Z(h,k) = j*2*pi*Frequence*4*pi*10^(-4)/2/pi*log(2*(LinePos(h,2)+P)/r)+Zc(h,k);  % ������-4�η���Ϊ��������Zc��λһ��
            else
                m = sqrt(2*pi*Frequence*1.00022*mui*dao_AL);
                delta = 1/sqrt(pi*Frequence*1.000022*mui*dao_AL);
                r = rM(h);
                q = r-(r/delta-0.1)*delta;
                mmc = sqrt(j)*m*q;
                
                k1 = (besseli(0,sqrt(j)*(m*r+0.001))-besseli(0,sqrt(j)*m*r))/0.001;
   %             k1 = difbesseli0j(mmc);
                k2 = (besselk(0,sqrt(j)*(m*r+0.001))-besselk(0,sqrt(j)*m*r))/0.001;
   %             k2 = difbesselk0j(mmc);
                fai = -(besseli(0,sqrt(j)*(m*q+0.001))-besseli(0,sqrt(j)*m*q))/(besselk(0,sqrt(j)*(m*q+0.001))-besselk(0,sqrt(j)*m*q));
   %             fai = -k1/k2;
            
                Zc(h,k) = rhoCurrentLine(h-n2)*j*0.5*m*r*(besseli(0,sqrt(j)*m*r)+fai*besselk(0,sqrt(j)*m*r))/(k1+fai*k2)*(1-(q/r)^2);
                Z(h,k) = j*2*pi*Frequence*4*pi*10^(-4)/2/pi*log(2*(LinePos(h,2)+P)/r)+Zc(h,k);
            end
        else
              Z(h,k) = j*2*pi*Frequence*4*pi*10^(-4)/2/pi*log(sqrt(((LinePos(h,1)-LinePos(k,1))^2+(LinePos(h,2)+LinePos(k,2)+2*P)^2)/((LinePos(h,1)-LinePos(k,1))^2+(LinePos(h,2)-LinePos(k,2))^2)));
        end
    end
end

Zl = Z;
Zll = Zl(1:n2, 1:n2);
Zlt = Zl(1:n2, (n2+1):n);
Ztl = Zl((n2+1):n,1:n2);
Ztt = Zl((n2+1):n,(n2+1):n);
Zee = Ztt - Ztl*inv(Zll)*Zlt;

for k = 1:phase
    Noi = phaseinfo(k);
    for i = (Noi+1): phaseinfo(k+1)-1
        Zee(:,i) = Zee(:,i) - Zee(:,Noi);
    end
end

for k = 1:phase
    Noi = phaseinfo(k); 
    for i = (Noi+1): (phaseinfo(k+1) -1)
        Zee(i,:) = Zee(i,:) - Zee(Noi,:);
    end
end

for k = 1:phase
    Noi = phaseinfo(k); % 
    Zeeend(k,:) = Zee(Noi,:);
    a1 = phase + phaseinfo(k)-(k-1);
    b1 = a1 + phaseinfo(k+1)-phaseinfo(k) - 2;
    a2 = Noi + 1;
    b2= a2 + phaseinfo(k+1)-phaseinfo(k) - 2;
    Zeeend(a1:b1, :) = Zee(a2:b2, :);  
end


for k = 1:phase
    Noi = phaseinfo(k); % 
    Zee_(:,k) = Zeeend(:,Noi);
    a1 = phase + phaseinfo(k)-(k-1);
    b1 = a1 + phaseinfo(k+1)-phaseinfo(k) - 2;
    a2 = Noi + 1;
    b2 = a2 + phaseinfo(k+1)-phaseinfo(k) - 2;
    Zee_( :,a1:b1) = Zeeend(:,a2:b2);  
end

Zell = Zee_(1:phase,1:phase);
Zelt = Zee_(1:phase, phase+1:size(Zee_,1));
Zetl = Zee_(phase+1:size(Zee_,1),1:phase);
Zett = Zee_(phase+1:size(Zee_,1),phase+1:size(Zee_,1));
Z = Zell - Zelt*inv(Zett)*Zetl;

% adjust the unit
Z = Z ./ 1000;

% http://www.uma.ac.ir/files/site1/a_akbari_994c8e8/bessel.pdf for modified
% bessel functions.
function div = difbesselj0(z)
div = besseli(-1,z);
return
function divj = difbesselj0j(z)
j = sqrt(-1);
divj = besseli(-1,z) * (1+j)/sqrt(2);
return

function div = difbesselk0(z)
div = -besselk(-1,z);
return
function divj = difbesselk0j(z)
j=sqrt(-1);
divj = -besselk(-1,z) * (1+j)/sqrt(2);
return

function div = difbesseli0(z)
div = besseli(-1,z);
return
function divj = difbesseli0j(z)
j=sqrt(-1);
divj = besseli(-1,z)*(1+j)/sqrt(2);
return