% main function of the package.
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code Written by chijie@tsinghua.edu.cn

function main(id)
 warning('off','MATLAB:dispatcher:InexactMatch');
 filename = 'config.dat';
 [route_start, phase_start, ...
          exchange,diameter, om, subconductor, sag, voltage, max_voltage, current, ...
          ground_bground, ground_diameter, ground_om, ground_wire, ground_sag,...
          rouearth,...
          rough,...
          aptitude, weather, ...
          repeatdef,...
          ri_an_def,...
          selfdef] ...
         = read_config( filename);
   
 apt_start = repeatdef(1,1); apt_end = repeatdef(1,2); apt_step = repeatdef(1,3); 
 ph_dis_start = repeatdef(2,1); ph_dis_end = repeatdef(2,2); ph_dis_step = repeatdef(2,3);    
 height_start = repeatdef(3,1); height_end = repeatdef(3,2); height_step = repeatdef(3,3);
 
 an_apt_base = ri_an_def(1,1); an_apt_addvalue = ri_an_def(1,2);
 ri_apt_base = ri_an_def(2,1); ri_apt_addvalue = ri_an_def(2,2);
 
 ef_height = selfdef(1,1); bf_height = selfdef(1,2);
 ri_height = selfdef(2,1); ri_width = selfdef(2,2);
 an_height = selfdef(3,1); an_width = selfdef(3,2); 
 
 total_route = length(route_start)-1;
 total_ground = length(ground_diameter);
 sag_total = [ground_sag; sag];
 ph = reshape( phase_start.', [], 1);
 n_div_number = ph(total_route*3+1:-1:2) - ph(total_route*3:-1:1);
 total_phase = total_route*3;
 phase_radius = zeros(total_phase,1);
 for i=1:total_phase
    phase_radius(i) = mean( diameter(ph(i):ph(i+1)-1)) * 0.5;
    phase_center(i,:) = mean( subconductor(ph(i):ph(i+1)-1,:));
 end
 
 for apt = apt_start : apt_step : apt_end
	for height = height_start:height_step:height_end
		for phase_dis = ph_dis_start: ph_dis_step: ph_dis_end
              
                u = [ zeros(total_ground,1); voltage];
                linecenter = [ ground_wire; subconductor];
                rline = [ ground_diameter; diameter] ./2;
                [p,q] = qline( u,linecenter,rline);    
            
                xpos = (-50:0.2:50)'; 
                E = eg( q, linecenter, sag_total, xpos, ef_height);
                E = E ./1e3;  % v/m-> kv/m
                Emax = max(E);                
              
                B = bg( current, subconductor, sag, xpos, bf_height);
                B = B .* 1e6;  % uT                
          
                [ef, eground, econmax, econmean] = esurf( q, linecenter, rline, length(route_start)-1, phase_start);        
                
                freq = 500e3;
                [Y,Z] = YZ( subconductor, diameter .* 0.5, om, ground_wire, ground_diameter .* 0.5, ground_om, freq, rouearth, total_route, phase_start);                         
                % RI
                if ( max(n_div_number) >= 4)
                     method = 0; rifm1 = ri( method, phase_center, n_div_number, phase_radius .* 2, Y, Z, econmean, rouearth, xpos, ri_height) - 12.5;    
                     method = 1; rifm2 = ri( method, phase_center, n_div_number, phase_radius .* 2, Y, Z, econmean, rouearth, xpos, ri_height) - 12.5;   
                else
                    rifm1 = ristd( econmax, phase_center, phase_radius, xpos, ri_height);  % ���귨
                    rifm2 = rifm1;
                end

                anbpa = noisebpa( econmean, phase_center, n_div_number, phase_radius, xpos, an_height);
                an345 = noise345( econmean, phase_center, n_div_number, phase_radius, xpos, an_height);
              
                [coronainit, ginit, einit_west, eginit_west] = einit( aptitude,  phase_radius.* 2, ground_diameter, rough);    
                
                [aveloss, peakloss] = corona( econmean, coronainit, rough, n_div_number, (phase_radius .* 100), ...
                                             weather(1), weather(2),weather(3),weather(4));
                                         
                exporttoxls( xpos, E, B, rifm1, rifm2, anbpa, an345, econmax, eground,...
                    coronainit, ginit, einit_west, eginit_west, aveloss, peakloss,...
                    0, 0);    % natpow, waveimp              
        end
    end
 end