function Comparison_Caso_Ideale
% clear all; 
close all; clc;

%% CONSTANTS

NN = 100;                   % Punti/2 usati per la creazione delle curve
n = 50;                     % Numero particelle
phi_f = 60*pi/180;          % Angolo finale
u = 1;                      % Massimo controllo

%% ISTANTE DI SWITCH

% known_switch = 1 --> Lo switch da u = 1 a u = -1 avviene esattamente a
% metà intervallo temporale
% known_switch = 0 --> Lo switch da u = 1 a u = -1 viene trovato con il PSO

known_switch = 1;

%% SCELTA DELLA FUNZIONE DA APPROSSIMARE

% a. CASE = 1 --> approssimo l'angolo
% b. CASE = 2 --> approssimo il Parametro Modificato di Rodriguez

CASE = 1;

if CASE == 1
    str = sprintf('Si sta approssimando l''angolo');
elseif CASE == 2
    str = sprintf('Si sta approssimando l''MRP');
end

disp(str);

%% SCELTA DEL POLINOMIO APPROSSIMANTE

% a. POL = 1 --> approssimo con BSplines
% b. POL = 2 --> approssimo con Polinomi di Chebyschev

POL = 1;

if POL == 1
    str = sprintf('Si sta approssimando BSplines');
    m = 15;                     % Numero punti per particella
    k = 7;                      % (Ordine-1) dei polinomi interpolanti
elseif POL == 2
    str = sprintf('Si sta approssimando con Polinomi di Chebyshev');
    m = 3;                     % Ordine massimo dei polinomi di Chebyshev
    check = 0;
end

disp(str);

%% IDEAL CASE - OPTIMAL STRUCTURE

% Final time
TF = 2*sqrt(phi_f/u);
T1 = TF/2;

% Time vectors
tt1 = linspace(0,T1,NN);
tt2 = linspace(T1,TF,NN);
tt = [tt1 tt2(2:end)];

% Angle History
phi1 = 0.5*u*tt1.^2;
phi2 = 2*u*T1*tt2 - u*T1^2 - 0.5*u*tt2.^2;
phi = [phi1 phi2(2:end)];

% Optimal MRP 
MRP1 = tan(phi/4);

% Velocity History
v1 = u*tt1;
v2 = 2*u*T1-u*tt2;
v = [v1 v2(2:end)];
  
% Optimal control
u1 = ones(1,NN);
optimal_control = [u1 -u1(2:end)];

figure
plot(tt,phi)
title('Optimal angle')

figure
plot(tt,MRP1)
title('Optimal MRP1')

figure
plot(tt,v)
title('Optimal Velocity')

%% Definizione condizioni iniziali e finali

if CASE == 1
    IC = [0,0,0]';
    FC = [phi_f;0;0];
elseif CASE == 2
    IC = [0,0,0]';
    FC = [tan(phi_f/4);0;0];
end

%% Creazione delle strutture per le particelle

% Le seguenti variabili sono create tramite la classe PSO_Class che
% riassume tutte le caratteristiche fondamentali dello swarm PSO.ù

swarm(1:n) = PSO_Class(m);
vel(1:n) = PSO_Class(m);
pbest = swarm;
gbest(1:n) = PSO_Class(m);
Gbest = PSO_Class(m);

%% INIZIALIZZAZIONE PARTICELLE

% Attenzione: velocità troppo elevate non portano a convergenza!!
% Valori per i parametri di rodriguez
if CASE == 1
    mrp_l_limit = 0; 
    mrp_u_limit = 2*phi_f; 
elseif CASE == 2
    mrp_l_limit = -1.5*tan(phi_f/4); 
    mrp_u_limit =  1.5*tan(phi_f/4); 
end

v_max_mrp = 0.01*(mrp_u_limit - mrp_l_limit); %0.01

% Valori per il tempo
if known_switch == 1
    t_max = TF;
    t_min = TF;
else
    % Uso il campo 'tempo' come costante di switch
    t_max = 1;
    t_min = 0;
end

k_time = 0.001; %0.01

v_max_t = k_time*(2*t_max);

% Initialization of the swarm and velocities
for i = 1:n
    if known_switch == 1
        swarm(i) = initialize_swarm(swarm(i),m,mrp_l_limit,mrp_u_limit,TF,TF,IC,FC,POL);
    else
        swarm(i) = initialize_swarm(swarm(i),m,mrp_l_limit,mrp_u_limit,0,1,IC,FC,POL);
    end
    vel(i) = initialize_vel(vel(i),m,v_max_mrp,v_max_t,POL);
    swarm(i).mrp3 = zeros(1,m);
    vel(i).mrp3 = zeros(1,m);
end

%%  STRUTTURA DELLE NON-UNIFORM CLAMPED B-SPLINE

%%% Fattore moltiplicativo della lunghezza dei vettori in output
L = 2*NN;
if POL == 1
    [N,t] = define_bsplines(L,k,m);
    
    N0 = N(k).valori(:,:);
    N1 = N(k-1).valori(:,2:end-1);
    N2 = N(k-2).valori(:,3:end-2);
elseif POL == 2    
    T = chebyshev_initialization(m,NN,check);
    tau = linspace(-1,1,NN);
    
    T0 = 0;    
    tempo1 = T0 + 0.5*(tau+1)*(TF/2-T0);
    tempo2 = TF/2 + 0.5*(tau+1)*(TF-TF/2);
    tempo = [tempo1 tempo2(2:end)];
end

goal = 1e-14;
move_times = 0;
diff_times = 0;
only_int = 1;

%%% Numero di iterate interne ed esterne

N_int = 20; %nel paper = 20

N_ext = 100; %50

%%% Variabili PSO

J = 1e10*ones(n,1);
J_g = 1e10*ones(n,1);
J_G = 1e10;

%%% Moltiplicatori della velocità

w0 = 1.2;%1.2;          1
wf = 0.6;%0.4; Più è alto più aumenta il carattere random.
%     Rischio di perdere tempo e di perdere la soluzione ottima
%
r0 = 0.9;
rf = 0.9;

c_p = 1.5;%1.5

%%%% Aggiunto il 25-03-2015
c_l0 = 2;
c_g0 = 0;

c_lf = 0;
c_gf = 2;

%%% Salvataggio precisione
errore = 1e10;
errore_best = 1e10;
i = 0;

%% Avvio computazione

while errore_best>goal && i <= 100
    
    i = i+1;
    
    c_l = c_l0 - (c_l0-c_lf)*(i-1)/N_ext;
    c_g = c_g0 - (c_g0-c_gf)*(i-1)/N_ext;
    
    if c_l < c_lf
        c_l = c_lf;
    end
    
    if c_g > c_gf
        c_g = c_gf;
    end
    
    %%%%
    
    disp(' ')
    disp('================================================================')
    disp(' ')
    st1=sprintf('%s','CICLO ESTERNO, ITERATA NUMERO: ', num2str(i));
    disp(st1)
    disp(' ')
    if i == 1
        st=sprintf('Primo ciclo');
        disp(st)
    else
        st2=sprintf('Errore rispetto al tempo ottimo: %g', errore_best);
        disp(st2)
    end
    
    if i > 1
        
        close all
        
        figure()
        
%         subplot(4,1,1)
%         hold on
%         plot(tt,phi)
%         plot(tt,angle_best,'r')
%         title('Angle Bspline Approximation')
%         
%         subplot(4,1,2)
%         hold on
%         plot(tt,MRP1)
%         plot(tt,mrp_best,'r')
%         title('MRP Approximation')
%         
%         subplot(4,1,3)
%         hold on
%         plot(tt,v)
%         plot(tt,omega_best,'r')
%         title('Velocity Approximation')
%         
%         subplot(4,1,4)
%         hold on
%         plot(tt,optimal_control)
%         plot(tt,u_best,'r')
%         title('Control Approximation')
        
        figure
        hold on
        plot(tt,phi,'o')
        plot(tt,angle_best,'ro')
        
        pause(0.4)
    end
    
    disp(' ')
    disp('================================================================')
    
    %% INIZIO DEL CICLO PER L'EVOLUZIONE DELLO SWARM
    
    for iii = 1:N_int
        
        % Valutazinoe dei coefficienti di restrizione e di inerzia
        r = r0 - (r0-rf)*(i-1)/N_ext - (r0-rf)/N_ext*(iii-1)/N_int;
        w = w0 - (w0-wf)*(i-1)/N_ext - (w0-wf)/N_ext*(iii-1)/N_int;
        
        w(w<wf) = wf;
        r(r<rf) = rf;
        
        for j = 1:n
            
            if CASE == 1
                if POL == 1
                    [tempo,angle,omega,control] = bspline_interpolation_struct_trad(swarm(j).tt,swarm(j).t1,swarm(j).t2,swarm(j).t3,...
                        swarm(j).tempo,swarm(j).mrp1,swarm(j).mrp2,swarm(j).mrp3,N0,N1,N2,t,k,m,diff_times);
                    angle1 = nakeinterp1(tempo(:,1),angle(1,:)',tt');
                    errore = sum((phi-angle1').^2);
                elseif POL == 2
                    angle1   = chebyshev_approximation(T,swarm(j).mrp1,check);   
                    angle2   = chebyshev_approximation(T,swarm(j).mrp2,check);   
                    angle = [angle1 angle2(2:end)];
                    
                    if known_switch == 0
                        ttt1 = linspace(T0,TF*swarm(j).tempo,NN);
                        ttt2 = linspace(TF*swarm(j).tempo,TF,NN);
                        ttt = [ttt1 ttt2(2:end)];
                        angle = nakeinterp1(ttt',angle',tempo');
                        angle = angle';
                    end
                    
                    mrp1 = tan(angle/4);
                    
                    coeff_omega1 = zeros(m+1,1);
                    for nn = (m-1):-1:2
                        coeff_omega1(nn) = coeff_omega1(nn + 2) + 2*(nn)*swarm(j).mrp1(nn+1);
                    end
                    coeff_omega1(1) = 0.5*coeff_omega1(3) + swarm(j).mrp1(2);
                    
                    omega1   = (2/(TF/2))*chebyshev_approximation(T,coeff_omega1(1:m),check);
                    
                    coeff_omega2 = zeros(m+1,1);
                    for nn = (m-1):-1:2
                        coeff_omega2(nn) = coeff_omega2(nn + 2) + 2*(nn)*swarm(j).mrp2(nn+1);
                    end
                    coeff_omega2(1) = 0.5*coeff_omega2(3) + swarm(j).mrp2(2);
                    
                    omega2   = (2/(TF/2))*chebyshev_approximation(T,coeff_omega2(1:m),check);
                    omega = [omega1 omega2(2:end)];
                    
                    if known_switch == 0
                        omega = nakeinterp1(ttt',omega',tempo');
                        omega = omega';
                    end
                    
                    coeff_control1 = zeros(m,1);
                    for nn = (m-2):-1:2
                        coeff_control1(nn) = coeff_control1(nn + 2) + 2*(nn)*coeff_omega1(nn+1);
                    end
                    coeff_control1(1) = 0.5*coeff_control1(3) + coeff_omega1(2);
                    control1   = (2/(TF/2))^2*chebyshev_approximation(T,coeff_control1(1:m),check);
                    
                    coeff_control2 = zeros(m,1);
                    for nn = (m-2):-1:2
                        coeff_control2(nn) = coeff_control2(nn + 2) + 2*(nn)*coeff_omega2(nn+1);
                    end
                    coeff_control2(1) = 0.5*coeff_control2(3) + coeff_omega2(2);
                    control2   = (2/(TF/2))^2*chebyshev_approximation(T,coeff_control2(1:m),check);
                    control = [control1 control2(2:end)];
                    
                    if known_switch == 0
                        control = nakeinterp1(ttt',control',tempo');
                        control = control';
                    end
                    
                    errore = sum((phi-angle).^2) + sum(abs(control)>1);
                
                end
                
            elseif CASE == 2
                
                if POL == 1
                    [tempo,mrp,~,~] = bspline_interpolation_struct_trad_mex(swarm(j).tt,swarm(j).t1,swarm(j).t2,swarm(j).t3,...
                        swarm(j).tempo,swarm(j).mrp1,swarm(j).mrp2,swarm(j).mrp3,N0,N1,N2,t,k,m,diff_times);
                    mrp = nakeinterp1(tempo(:,1),mrp(1,:)',tt');
                    errore = sum((MRP1-mrp').^2);
                elseif POL == 2
                    mrp11   = chebyshev_approximation(T,swarm(j).mrp1,check);   
                    mrp12   = chebyshev_approximation(T,swarm(j).mrp2,check);   
                    mrp1 = [mrp11 mrp12(2:end)];
                    
                    angle1_coeff = 4*atan(swarm(j).mrp1);
                    angle2_coeff = 4*atan(swarm(j).mrp2);
                    
                    angle = 4*atan(mrp1);
                    
                    coeff_omega1 = zeros(m+1,1);
                    for nn = (m-1):-1:2
                        coeff_omega1(nn) = coeff_omega1(nn + 2) + 2*(nn)*angle1_coeff(nn+1);
                    end
                    coeff_omega1(1) = 0.5*coeff_omega1(3) + angle1_coeff(2);
                    
                    omega1   = (2/(TF/2))*chebyshev_approximation(T,coeff_omega1(1:m),check);
                    
                    coeff_omega2 = zeros(m+1,1);
                    for nn = (m-1):-1:2
                        coeff_omega2(nn) = coeff_omega2(nn + 2) + 2*(nn)*angle2_coeff(nn+1);
                    end
                    coeff_omega2(1) = 0.5*coeff_omega2(3) + angle2_coeff(2);
                    
                    omega2   = (2/(TF/2))*chebyshev_approximation(T,coeff_omega2(1:m),check);
                    omega = [omega1 omega2(2:end)];
                    
                    coeff_control1 = zeros(m,1);
                    for nn = (m-2):-1:2
                        coeff_control1(nn) = coeff_control1(nn + 2) + 2*(nn)*coeff_omega1(nn+1);
                    end
                    coeff_control1(1) = 0.5*coeff_control1(3) + coeff_omega1(2);
                    control1   = (2/(TF/2))^2*chebyshev_approximation(T,coeff_control1(1:m),check);
                    
                    coeff_control2 = zeros(m,1);
                    for nn = (m-2):-1:2
                        coeff_control2(nn) = coeff_control2(nn + 2) + 2*(nn)*coeff_omega2(nn+1);
                    end
                    coeff_control2(1) = 0.5*coeff_control2(3) + coeff_omega2(2);
                    control2   = (2/(TF/2))^2*chebyshev_approximation(T,coeff_control2(1:m),check);
                    control = [control1 control2(2:end)];
                    
                    errore = sum((MRP1-mrp1).^2);
                    
                end
            end
            
            % Valutazione FUNZIONE DI FITNESS
            fitness = errore;
            
            % Valutazione del Personal Best
            if fitness < J(j)
                J(j) = fitness;
                pbest(j) = swarm(j);
            end
            
            % Valutazione del Global Best
            if  J(j) < J_G
                J_G = J(j);
                Gbest = pbest(j);
                errore_best = fitness;
                if POL == 1
                    angle_best  = nakeinterp1(tempo(:,1),angle(1,:)',tt');
                    mrp_best    = nakeinterp1(tempo(:,1),tan(0.25*angle(1,:))',tt');
                    omega_best  = nakeinterp1(tempo(:,1),omega(1,:)',tt');
                    u_best      = nakeinterp1(tempo(:,1),control(1,:)',tt');
                elseif POL == 2
                    angle_best  = angle;
                    mrp_best    = mrp1;
                    omega_best  = omega;
                    u_best      = control;
                end
                
            end
        end
        
        
        %% SPOSTAMENTO DELLO SWARM + LOCAL BEST
        
        for j=1:n
            
            % Valutazione del Local Best
            gbest(j) = local_best_search(j,n,J,J_g,pbest); 
            
            % Valutazione delle velocità
            vel(j) = vel_evaluation(m,vel(j),swarm(j),pbest(j),gbest(j),Gbest,r,w,c_p,c_l,c_g,v_max_mrp,v_max_t,move_times,diff_times,only_int,POL);
            
            % Valutazione degli spostamenti
            swarm(j) = displacement(vel(j),swarm(j),m,mrp_u_limit,t_max,t_min,move_times,diff_times); 
            
        end
        
    end
end

end

function local_best = local_best_search(j,n,J,J_g,pbest)

a = j-2;
a(a==-1)=n-1;
a(a==0)=n;
b = j-1;
b(b==0)=n;
c = j+1;
c(c==n+1)=1;
d = j + 2;
d(d==n+1)=1;
d(d==n+2)=2;

J_g(j) = min([J_g(j),J(j),J(a),J(b),J(c),J(d)]);

if sum(J == J_g(j)) == 1
    local_best = pbest(J == J_g(j));
else
    d = (1:length(J))';
    local_best = pbest( min(d(J == J_g(j))));
end

end