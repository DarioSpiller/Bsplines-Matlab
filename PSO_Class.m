classdef PSO_Class
    
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mrp1
        mrp2
        mrp3
        tempo
        tt
        t1
        t2
        t3
    end
    
    methods
        
        function swarm = PSO_Class(m)
            if nargin == 1 && m > 0
                swarm.mrp1 = zeros(1,m);
                swarm.mrp2 = zeros(1,m);
                swarm.mrp3 = zeros(1,m);
                swarm.tempo = 0;
                swarm.tt = zeros(1,m);
                swarm.t1 = zeros(1,m);
                swarm.t2 = zeros(1,m);
                swarm.t3 = zeros(1,m);
            else
                error('Wrong value of m!')
            end
        end
        
        function swarm = initialize_swarm(swarm,m,mrp_l_limit,mrp_u_limit,t_min,t_max,IC,FC,POL)
            
            if POL == 1
                %%% Punti intermedi Modified Rodriguez Parameters
                swarm.mrp1(:,3:m-2) = mrp_l_limit*ones(1,m-4) + rand(1,m-4)*(mrp_u_limit - mrp_l_limit);
                swarm.mrp2(:,3:m-2) = mrp_l_limit*ones(1,m-4) + rand(1,m-4)*(mrp_u_limit - mrp_l_limit);
                swarm.mrp3(:,3:m-2) = mrp_l_limit*ones(1,m-4) + rand(1,m-4)*(mrp_u_limit - mrp_l_limit);
                
                %%% Punti iniziali e finali dello swarm, per garantire derivata nulla
                %%% iniziale
                swarm.mrp1(:,1:2) = IC(1);
                swarm.mrp2(:,1:2) = IC(2);
                swarm.mrp3(:,1:2) = IC(3);
                swarm.mrp1(:,m-1:m) = FC(1);
                swarm.mrp2(:,m-1:m) = FC(2);
                swarm.mrp3(:,m-1:m) = FC(3);
            elseif POL == 2
                swarm.mrp1(:,1:m) = mrp_l_limit*ones(1,m) + rand(1,m)*(mrp_u_limit - mrp_l_limit);
                swarm.mrp2(:,1:m) = mrp_l_limit*ones(1,m) + rand(1,m)*(mrp_u_limit - mrp_l_limit);
                swarm.mrp3(:,1:m) = mrp_l_limit*ones(1,m) + rand(1,m)*(mrp_u_limit - mrp_l_limit);
            end
            
            %%% Tempo di manovra
            swarm.tempo = t_min + rand(1)*(t_max - t_min);
            swarm.tt = linspace(0,swarm.tempo,m);
            swarm.t1 = linspace(0,1,m);
            swarm.t2 = linspace(0,1,m);
            swarm.t3 = linspace(0,1,m);
            
        end
        
        function vel = initialize_vel(vel,m,v_max_mrp,v_max_t,POL)
            
            if POL == 1
                %%% Punti intermedi Modified Rodriguez Parameters
                vel.mrp1(:,3:m-2) = v_max_mrp*(-1+2*rand(1,m-4));
                vel.mrp2(:,3:m-2) = v_max_mrp*(-1+2*rand(1,m-4));
                vel.mrp3(:,3:m-2) = v_max_mrp*(-1+2*rand(1,m-4));
                
                % Per garantire manovra rest-to-rest
                vel.mrp1(:,1:2) = 0;
                vel.mrp2(:,1:2) = 0;
                vel.mrp3(:,1:2) = 0;
                vel.mrp1(:,m-1:m) = 0;
                vel.mrp2(:,m-1:m) = 0;
                vel.mrp3(:,m-1:m) = 0;
            elseif POL == 2
                vel.mrp1(:,1:m) = v_max_mrp*(-1+2*rand(1,m));
                vel.mrp2(:,1:m) = v_max_mrp*(-1+2*rand(1,m));
                vel.mrp3(:,1:m) = v_max_mrp*(-1+2*rand(1,m));
            end
            
            %%% Tempo di manovra
            vel.tempo = v_max_t*(-1+2*rand(1));
            vel.tt(2:m) = v_max_t*(-1+2*rand(1,m-1));
            vel.tt(1) = 0;
            vel.t1 = [0 v_max_t/100*(-1+2*rand(1,m-2)) 0];
            vel.t2 = [0 v_max_t/100*(-1+2*rand(1,m-2)) 0];
            vel.t3 = [0 v_max_t/100*(-1+2*rand(1,m-2)) 0];
        end
        
        function vel = vel_evaluation(m,vel,swarm,pbest,gbest,Gbest,r,w,c_p,c_l,c_g,v_max_mrp,v_max_t,move_times,diff_times,only_int,POL)
            
            vel.mrp1 = r*(w*vel.mrp1...
                + c_p*rand(1)*(pbest.mrp1 - swarm.mrp1)...
                + c_l*rand(1)*(gbest.mrp1 - swarm.mrp1)...
                + c_g*rand(1)*(Gbest.mrp1 - swarm.mrp1));
            
            vel.mrp2 = r*(w*vel.mrp2...
                + c_p*rand(1)*(pbest.mrp2 - swarm.mrp2)...
                + c_l*rand(1)*(gbest.mrp2 - swarm.mrp2)...
                + c_g*rand(1)*(Gbest.mrp2 - swarm.mrp2));
            
            vel.mrp3 = r*(w*vel.mrp3...
                + c_p*rand(1)*(pbest.mrp3 - swarm.mrp3)...
                + c_l*rand(1)*(gbest.mrp3 - swarm.mrp3)...
                + c_g*rand(1)*(Gbest.mrp3 - swarm.mrp3));
            
            if POL == 1
                vel.mrp1(:,1:2) = 0;
                vel.mrp2(:,1:2) = 0;
                vel.mrp3(:,1:2) = 0;
                vel.mrp1(:,m-1:m) = 0;
                vel.mrp2(:,m-1:m) = 0;
                vel.mrp3(:,m-1:m) = 0;
            end
            
            if move_times == 1
                % Sub-intervals times
                if diff_times == 1
                    if only_int == 1
                        vel.tempo = zeros(size(vel.tempo));
                    else
                        vel.tempo = r*(w*vel.tempo...
                            + c_p*rand(1)*(pbest.tempo - swarm.tempo)...
                            + c_l*rand(1)*(gbest.tempo - swarm.tempo)...
                            + c_g*rand(1)*(Gbest.tempo - swarm.tempo));
                    end
                    vel.t1 = r*(w*vel.t1...
                        + c_p*rand(1)*(pbest.t1 - swarm.t1)...
                        + c_l*rand(1)*(gbest.t1 - swarm.t1)...
                        + c_g*rand(1)*(Gbest.t1 - swarm.t1));
                    vel.t2 = r*(w*vel.t2...
                        + c_p*rand(1)*(pbest.t2 - swarm.t2)...
                        + c_l*rand(1)*(gbest.t2 - swarm.t2)...
                        + c_g*rand(1)*(Gbest.t2 - swarm.t2));
                    vel.t3 = r*(w*vel.t3...
                        + c_p*rand(1)*(pbest.t3 - swarm.t3)...
                        + c_l*rand(1)*(gbest.t3 - swarm.t3)...
                        + c_g*rand(1)*(Gbest.t3 - swarm.t3));
                    vel.t1(:,1) = 0;
                    vel.t2(:,1) = 0;
                    vel.t3(:,1) = 0;
                    vel.t1(:,m) = 0;
                    vel.t2(:,m) = 0;
                    vel.t3(:,m) = 0;
                else
                    vel.tt = r*(w*vel.tt...
                        + c_p*rand(1)*(pbest.tt - swarm.tt)...
                        + c_l*rand(1)*(gbest.tt - swarm.tt)...
                        + c_g*rand(1)*(Gbest.tt - swarm.tt));
                    
                end
                
            else
                % Using only the global time
                vel.tempo = r*(w*vel.tempo...
                    + c_p*rand(1)*(pbest.tempo - swarm.tempo)...
                    + c_l*rand(1)*(gbest.tempo - swarm.tempo)...
                    + c_g*rand(1)*(Gbest.tempo - swarm.tempo));
                
            end
            
            %%% Controllo
            vel.mrp1(abs(vel.mrp1)>v_max_mrp) = v_max_mrp*sign(vel.mrp1(abs(vel.mrp1)>v_max_mrp));
            vel.mrp2(abs(vel.mrp2)>v_max_mrp) = v_max_mrp*sign(vel.mrp2(abs(vel.mrp2)>v_max_mrp));
            vel.mrp3(abs(vel.mrp3)>v_max_mrp) = v_max_mrp*sign(vel.mrp3(abs(vel.mrp3)>v_max_mrp));
            
            if move_times == 1
                if diff_times == 1
                    v_max1 = v_max_t/Gbest.tempo;
                    
                    vel.tempo(abs(vel.tempo)>v_max_t) = v_max_t*sign(vel.tempo);
                    
                    check = abs(vel.t1) > v_max1;
                    vel.t1(check) = v_max1*sign(vel.t1(check));
                    
                    check = abs(vel.t2) > v_max1;
                    vel.t2(check) = v_max1*sign(vel.t2(check));
                    
                    check = abs(vel.t3) > v_max1;
                    vel.t3(check) = v_max1*sign(vel.t3(check));
                else
                    % Sub-intervals times
                    check = abs(vel.tt)>v_max_t;
                    vel.tt(check) = v_max_t*sign(vel.tt(check));
                end
            else
                % Using only the global time
                vel.tempo(abs(vel.tempo)>v_max_t) = v_max_t*sign(vel.tempo);
            end
            
        end
        
        function swarm = displacement(vel,swarm,m,mrp_u_limit,t_max,t_min,move_times,diff_times)
            
            %%% Spostamento
            swarm.mrp1 = swarm.mrp1 + vel.mrp1;
            swarm.mrp2 = swarm.mrp2 + vel.mrp2;
            swarm.mrp3 = swarm.mrp3 + vel.mrp3;
            
            
            if move_times == 1
                if diff_times == 1
                    swarm.tempo = swarm.tempo + vel.tempo;
                    swarm.t1 = swarm.t1 + vel.t1;
                    swarm.t2 = swarm.t2 + vel.t2;
                    swarm.t3 = swarm.t3 + vel.t3;
                else
                    % Sub-intervals times
                    swarm.tt = swarm.tt + vel.tt;
                end
            else
                % Using only the global time
                swarm.tempo = swarm.tempo + vel.tempo;
            end
            
            %%% Controllo
            swarm.mrp1(abs(swarm.mrp1)>mrp_u_limit) = mrp_u_limit*sign(swarm.mrp1(abs(swarm.mrp1)>mrp_u_limit));
            swarm.mrp2(abs(swarm.mrp2)>mrp_u_limit) = mrp_u_limit*sign(swarm.mrp2(abs(swarm.mrp2)>mrp_u_limit));
            swarm.mrp3(abs(swarm.mrp3)>mrp_u_limit) = mrp_u_limit*sign(swarm.mrp3(abs(swarm.mrp3)>mrp_u_limit));
            
            
            if move_times == 1
                if diff_times == 1
                    % t1
                    swarm.t1(swarm.t1>1) = 1;
                    if swarm.t1(end) < 0;
                        swarm.t1(end) = 0;
                    end
                    swarm.t1(swarm.t1<0) = 0;
                    
                    if ~all(diff(swarm.t1)>0)
                        for ii = 2:m-1
                            if swarm.t1(ii) <= swarm.t1(ii-1)
%                                 swarm.t1(ii) = swarm.t1(ii-1) + 1e-10;
                                swarm.t1(ii) = 0.5*(swarm.t1(ii-1) + swarm.t1(ii+1));
                            end
                        end
                    end
                    
                    % t2
                    swarm.t2(swarm.t2>1) = 1;
                    if swarm.t2(end) < 0;
                        swarm.t2(end) = 0;
                    end
                    swarm.t2(swarm.t2<0) = 0;
                    
                    if ~all(diff(swarm.t2)>0)
                        for ii = 2:m-1
                            if swarm.t2(ii) <= swarm.t2(ii-1)
%                                 swarm.t2(ii) = swarm.t2(ii-1) + 1e-10;
                                swarm.t2(ii) = 0.5*(swarm.t2(ii-1) + swarm.t2(ii+1));
                            end
                        end
                    end
                    
                    % t3
                    swarm.t3(swarm.t3>1) = 1;
                    if swarm.t3(end) < 0;
                        swarm.t3(end) = 0;
                    end
                    swarm.t3(swarm.t3<0) = 0;
                    
                    if ~all(diff(swarm.t3)>0)
                        for ii = 2:m-1
                            if swarm.t3(ii) <= swarm.t3(ii-1)
%                                 swarm.t3(ii) = swarm.t3(ii-1) + 1e-10;
                                swarm.t3(ii) = 0.5*(swarm.t3(ii-1) + swarm.t3(ii+1));
                            end
                        end
                    end
                    
                    swarm.tempo(swarm.tempo>t_max) = t_max;
                    swarm.tempo(swarm.tempo<t_min) = t_min;
                    
                else
                    % Sub-intervals times
                    swarm.tt(swarm.tt>t_max) = t_max;
                    if swarm.tt(end) < t_min;
                        swarm.tt(end) = t_min;
                    end
                    swarm.tt(swarm.tt<0) = 0;
                    
                    if ~all(diff(swarm.tt)>0)
                        for ii = 2:m-1
                            if swarm.tt(ii) <= swarm.tt(ii-1)
                                %                                 swarm.tt(ii) = swarm.tt(ii-1) + 1e-10;
                                swarm.tt(ii) = 0.5*(swarm.tt(ii-1) + swarm.tt(ii+1));
                            end
                        end
                    end
                    
                    swarm.tempo = swarm.tt(end);
                end
            else
                swarm.tempo(swarm.tempo>t_max) = t_max;
                swarm.tempo(swarm.tempo<t_min) = t_min;
                swarm.tt = linspace(0,swarm.tempo,m);
            end
            
        end
        
        function new_swarm = new_bsplines_struct(k,m,t,swarm,tt,ni)
            
            % m e t sono già aggiornati per l'aggiunta del nuovo punto
            
            a = zeros(m,1);
            
            new_swarm(1:length(swarm)) = PSO_Class(m);
            
            % Ricalcolo tutti i punti per l'interpolazione
            for ind = 1:m
                
                if ind <= ni-k
                    a(ind) = 1;
                elseif (ind >= ni-k +1) && (ind<=ni-1)
                    a(ind) = (tt-t(ind))/(t(ind+k)-t(ind));
                else
                    a(ind) = 0;
                end
                
                for iii = 1:length(swarm)
                    
                    if ind < m
                        mrp1_a = swarm(iii).mrp1(ind);
                        mrp2_a = swarm(iii).mrp2(ind);
                        mrp3_a = swarm(iii).mrp3(ind);
                        tt_a = swarm(iii).tt(ind);
                        t1_a = swarm(iii).t1(ind);
                        t2_a = swarm(iii).t2(ind);
                        t3_a = swarm(iii).t3(ind);
                    else
                        mrp1_a = 0;
                        mrp2_a = 0;
                        mrp3_a = 0;
                        tt_a = 0;
                        t1_a = 0;
                        t2_a = 0;
                        t3_a = 0;
                    end
                    
                    if ind > 1
                        mrp1_b = swarm(iii).mrp1(ind-1);
                        mrp2_b = swarm(iii).mrp2(ind-1);
                        mrp3_b = swarm(iii).mrp3(ind-1);
                        tt_b = swarm(iii).tt(ind-1);
                        t1_b = swarm(iii).t1(ind-1);
                        t2_b = swarm(iii).t2(ind-1);
                        t3_b = swarm(iii).t3(ind-1);
                    else
                        mrp1_b = 0;
                        mrp2_b = 0;
                        mrp3_b = 0;
                        tt_b = 0;
                        t1_b = 0;
                        t2_b = 0;
                        t3_b = 0;
                    end
                    
                    new_swarm(iii).mrp1(ind) = a(ind)*mrp1_a + (1-a(ind))*mrp1_b;
                    new_swarm(iii).mrp2(ind) = a(ind)*mrp2_a + (1-a(ind))*mrp2_b;
                    new_swarm(iii).mrp3(ind) = a(ind)*mrp3_a + (1-a(ind))*mrp3_b;
                    new_swarm(iii).tt(ind)   = a(ind)*tt_a + (1-a(ind))*tt_b;
                    new_swarm(iii).tempo     = swarm(iii).tempo;
                    new_swarm(iii).t1(ind)   = a(ind)*t1_a + (1-a(ind))*t1_b;
                    new_swarm(iii).t2(ind)   = a(ind)*t2_a + (1-a(ind))*t2_b;
                    new_swarm(iii).t3(ind)   = a(ind)*t3_a + (1-a(ind))*t3_b;
                    
                end
                
            end
            
        end
        
    end
    
end

