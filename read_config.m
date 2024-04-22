% compute the audible noise. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function [route_start, phase_start, ...
          exchange,diameter, om, subconductor, sag, voltage, max_voltage, current, ...
          ground_bground, ground_diameter, ground_om, ground_wire, ground_sag,...
          rouearth,...
          rough,...
          aptitude, weather, ...
          repeatdef,...
          ri_an_def,...
          selfdef] ...
         = read_config( filename)
     
fid = fopen( filename, 'r');

% 
spe_route   = fscanf( fid, '%d', 1);   % 
route_start = ones( spe_route+1, 1);   % 
phase_start = ones( spe_route+1, 3);
exchange    = zeros(spe_route, 1);

angle = [2*pi 4*pi/3 2*pi/3];
angle = [0 2/3*pi 4/3*pi];
PI = cos(angle) + sin(angle) .* sqrt(-1); % 
PI = sin(angle) + cos(angle) .* sqrt(-1);
SQRT3 = sqrt(3);
x = 0;
for j = 1: spe_route                           
    [A, cnt]  = fscanf( fid, '%f', 5);
    exchange(j) = A(1);  % �Ƿ�λ
    max1        = A(2);  % ������е�ѹ                  
    mva1        = A(3);  % ��������
    sag1        = A(4);  % ���� 
    sys1        = A(5);  % ϵͳ��ѹ���ߵ�ѹ
    
    route_start(j+1) = route_start(j);  % ��׼
    for t = 1:3
        size_of_phase = fscanf( fid, '%d', 1);                    % ���������
        if ( t < 3) % t == 1 or 2
            phase_start(j, t+1) = phase_start(j, t) + size_of_phase;   % ��¼��һ������￪ʼ
        else  % t == 3
            phase_start(j+1, 1) = phase_start(j, 3) + size_of_phase;   % ��¼��һ������￪ʼ
        end
       
        current1 = mva1 / SQRT3 / .95 / sys1;                     % �������ֵ mva/kv, ��������0.95.  MV = sqrt(3) UI 
        [A, cnt]  = fscanf( fid, '%f', 4 * size_of_phase);
        A = reshape( A, 4, size_of_phase);  A = A';
        diameter(x+1:x+size_of_phase,1)       = A(:,1);
        om (x+1:x+size_of_phase,1)            = A(:,2);
        sag(x+1:x+size_of_phase,1)            = ones(size_of_phase, 1) .* sag1;
        subconductor(x+1:x+size_of_phase,1)   = A(:,3);    
        subconductor(x+1:x+size_of_phase,2)   = A(:,4);                          % �ҵ�� 
        current(x+1:x+size_of_phase,1)        = current1 /size_of_phase * PI(t); % ÿ���ĵ���
        voltage(x+1:x+size_of_phase,1)        = sys1 /sqrt(3) * PI(t);           % ��Եص�ѹ��Чֵ
        max_voltage(x+1:x+size_of_phase,1)    = max1 * PI(t);                    % ��Чֵ       
        x                                     = x + size_of_phase;  
        route_start(j+1)                      = route_start(j+1) + size_of_phase;  
    end           
end
subconductor(:,2) = subconductor(:,2) - 2.0/3 .* sag;                           % ͳһ������ƽ����

% ������Ϣ
total_ground = fscanf( fid, '%d',1);
[gg, cnt] = fscanf( fid, '%f', 6 * total_ground);
gg = reshape( gg, 6, total_ground);
gg = gg';
ground_bground   = gg(:,1);
ground_diameter  = gg(:,2);
ground_om        = gg(:,3);
ground_wire(:,1) = gg(:,4);
ground_sag       = gg(:,6);
ground_wire(:,2) = gg(:,5) - 2.0/3.0 .* gg(:,6);                                % ͳһ������ƽ����

% ȫ�����ɱ�׼��λ
ground_diameter = ground_diameter / 1e3;  % ��
voltage = voltage * 1e3;                  % ��
max_voltage = max_voltage * 1e3;          % ��
current = current * 1e3;                  % ��
diameter = diameter /1e3;                 % ��

% ��������
% aptitude, fogy, rou, rainy, snow, sunny
[A, cnt] = fscanf( fid, '%f', 6);
A(A<0) = 0;
aptitude = A(1);
v = [2 4 5 6];
weather = A(v);
rouearth = A(3);

% �Ƿ�ѭ������
g_bool_for = fscanf( fid, '%d', 1);
if ( g_bool_for == 1)
   [A, cnt] = fscanf( fid, '%f', 3 * 5);
   A = reshape(A, 3, 5); repeatdef = A';
   t = repeatdef(:,1);   repeatdef(:,1) = repeatdef(:,2);  repeatdef(:,2) = t;
   t = repeatdef(:,3);
   t(t<0) = 1;
   repeatdef(:,3) = t;  
else
   repeatdef = repmat([0 0 1], 5, 1);          
end

[A, cnt] = fscanf( fid, '%f', 11);
ri_an_def = [A(1) A(2); A(3) A(4)];
% �Զ������AN��RI�ĸ߶������
selfdef = [A(5) A(6); A(7) A(8); A(9) A(10)];
% ����ֲڶ�
rough = A(11);

fclose(fid);