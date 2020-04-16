% QUESTA FUNZIONE SERVE CONFRONTARE LA GENERAZIONE DELLE B-SPLINE TRAMITE
% ALGORITMO RICORSIVO DI DE-BOOR E ALGORITMO TRIANGOLARE DI CUI AL SITO
% http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-coef.html
% http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/

% Nota bene che il Knot-Vector è un vettore non-decrescente, il che vuol
% dire che knot successivi possono al più essere uguali, ma il knot i+1 non
% può essere minore del knot i.

clc; close all;
% clear all;

%% DICHIARAZIONI VARIABILI GLOBALI

example = 4;

L = 10000;
p = 4;             % p - 1 è l'ordine dei polinomi interpolanti
% Nota bene: per avere le basis-function uguali a 1 agli estremi
% dell'intervallo, devo avere multiple-knots all'inizio e alla fine in
% numero pari a p!
if example == 1
    m = 8;
    t = linspace(0,1,m+1);
elseif example == 2
    m = 9;             % m + 1 è il numero di knot points
    t = [ 0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1];
elseif example == 3
    m = 11;             % m + 1 è il numero di knot points
    t = [ 0, 0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1, 1];
elseif example == 4
    % Per simulare il paper di Guadalajara
    p = 8;  % IAC 8
    m = 15; % IAC 15
    
    t = 1:1:m+1;
    n = m - p + 1;
    
    t(t<=p) = 0;
    t(t>p & t<=n) = (1:n-p)/(n-p+1);
    t(t>n) = 1;
end
n = m - p + 1;     % n è il numero di polinomi di grado p-1 da combinare
                   % linearmente che posso ottenere da m + 1 knot points

% Rispetto al sito citato sopra, n è differente (il valore di n qui definito
% è pari a n+1 del sito). Ugualmente per p, quindi il p qui definito è pari
% a p+1 del sito. m è lo stesso.
                   
% Verifica
if n < 1
    error('attenzione! Non puoi ottenere polinomi di grado p-1')
end

% Verifica
check = m == (n - 1) + (p - 1) + 1;
if check ~= 1
    error('attenzione! il vincolo che definisce le b-splines non è rispettato')
end

diff_times = 1;
TF = 10;
tt = linspace(0,TF,n);
t1 = tt/TF;
t2 = tt/TF.*(1 + [1,TF/50*(1-2*rand(1,n-2)), 1]); 
t3 = tt/TF.*(1 + [1,TF/50*(1-2*rand(1,n-2)), 1]);
tempo_ref = linspace(0,TF,L);
% angle1 = [0,linspace(0,3,n-2),3];
angle1 = [10,10, 20*(1-2*rand(1,n-4)),-5,0];
angle2 = [10,10,11.1676595543435,-6.48070381406463,-5.13633918216165,-1.97484927044590,-5,0]; 
angle3 = t1;
phi_f = pi/3;

if mod(n,2) == 0
    t_central = mean(tt(n/2:n/2+1));
    theta1 = 0.5*phi_f/(t_central^2) * [tt(1:n/2) t_central].^2;
    theta2 = fliplr(theta1);
    theta2 = phi_f-theta2;
    ref_traj = tan([theta1(1:end-1) theta2(2:end)]/4);
else
    theta1 = 0.5*phi_f/(tt(floor(n/2)+1)^2) * tt(1:floor(n/2+1)).^2;
    theta2 = fliplr(theta1);
    theta2 = phi_f-theta2;
    ref_traj = tan([theta1 theta2(2:end)]/4);
end
angle01 = ref_traj;
angle01(2) = angle01(1);  angle01(end-1) = angle01(end);

%% DE-BOOR RECURSIVE ALGORITHM

% m+1 = numero di knots
% p = degree of the basis functions
% n + 1 = numero di basis-function di grado p
% m = n + p + 1;

tic
[N,t,u] = define_bsplines(L,p,n,t);
DE_BOOR_time = toc;
disp(' ')
disp('==== ALGORITMO DE BOOR =====')
disp(' ')
st = sprintf('Total time: %f',DE_BOOR_time);
disp(st);

N0 = N(p).valori(:,:);
N1 = N(p-1).valori(:,2:end-1);
N2 = N(p-2).valori(:,3:end-2);

%%% PLOT BASIS FUNCTIONS
for ind = 1:p
    figure
    plot(u,N(ind).valori); grid on
    hold on
    plot(t(t~=0),0,'o')
    text_title = sprintf('Basis-function ordine %u',ind-1);
    title(text_title)
end

%%% PLOT INTERPOLATIONS
%
% ESEMPIO CON VETTORE DEI TEMPI NON UNIFORME
% [tempo,angle,omega,control,qx,qy] = bspline_interpolation_struct_trad(tt,t1,t2,t3,...
%     TF,angle1,angle2,angle3,N0,N1,N2,t,p,n,diff_times);
% 
% figure
% axes('FontSize',16);
% hold on
% plot(t2*TF,angle2,'MarkerSize',30,'Marker','.','LineWidth',2,...
%     'Color',[0.5 0.5 0.5])
% plot(tempo(:,2),angle(2,:),'LineWidth',3,'Color',[0 0 0])
% grid on
% xlabel('\lambda','FontSize',18)
% ylabel('f(\lambda)','FontSize',16)
% title('B-Spline approximation','FontWeight','bold','FontSize',18)

[tempo,angle,omega,control] = bspline_interpolation_simple(u,...
    TF,angle1,angle2,angle3,N0,N1,N2,t,p,n);


figure
set(gcf,'units','normalized','outerposition',[0.5 0 0.4 0.7])
axes('Position',[0.13 0.285393258426966 0.775 0.608988764044945],...
    'FontSize',16);
hold on
plot(tt,angle2,'MarkerSize',30,'Marker','.','LineWidth',2,...
    'Color',[0.5 0.5 0.5])
plot(tempo,angle(2,:),'LineWidth',3,'Color',[0 0 0])
grid on
xlabel('\lambda','FontSize',18)
ylabel('f(\lambda)','FontSize',16)
title('B-Spline curve','FontWeight','bold','FontSize',18)
% Create textbox
for i= 1:D
    s{i,:} = sprintf('a_%u = %6.2f',i-1,angle2(i));
end
S = {[s{1},'    ',s{2},'    ',s{3},'    ',s{4}];
     [s{5},'    ',s{6},'    ',s{7},'    ',s{8}]};
annotation(gcf,'textbox',...
    [0.186748035673793 0.0134831460674158 0.669323928344198 0.142670499085446],...
    'String',S,...
    'FontSize',16,...
    'FitBoxToText','off',...
    'LineWidth',1,...
    'BackgroundColor',[1 1 1]);

D = 8;
NN=200; 
tempo = linspace(0,TF,NN);
tau=2*tempo/TF - 1;
T(1,:)=ones(1,NN);
T(2,:)=tau;
dTdtau(1,:)=zeros(1,NN);
dTdtau(2,:)=ones(1,NN);
ddTddtau(1,:)=zeros(1,NN);
ddTddtau(2,:)=zeros(1,NN);
ddTddtau(3,:)=4*ones(1,NN);
ddTddtau(4,:)=24*tau;
dtaudt=2/(TF);

AA = angle2(3:6);
[a4,a2,a3,a1]=coeff_a_i(AA,[10 0/dtaudt],[0 3.5/dtaudt],D);
AAA = [a1 a2 a3 a4 AA];
[x,xdot,xddot]=Traiettoria_dot_off(AAA,D,T,dTdtau,ddTddtau,tau,dtaudt);

figure
set(gcf,'units','normalized','outerposition',[0.5 0 0.4 0.7])
axes('Position',[0.13 0.285393258426966 0.775 0.608988764044945],...
    'FontSize',16);
plot(tempo,x,'LineWidth',3,'Color',[0 0 0])
grid on
xlabel('\lambda','FontSize',18)
ylabel('f(\lambda)','FontSize',16)
title('Chebyshev approximation','FontWeight','bold','FontSize',18)
% Create textbox
for i= 1:D
    s{i,:} = sprintf('a_%u = %6.2f',i,AAA(i));
end
S = {[s{1},'    ',s{2},'    ',s{3},'    ',s{4}];
     [s{5},'    ',s{6},'    ',s{7},'    ',s{8}]};
annotation(gcf,'textbox',...
    [0.186748035673793 0.0134831460674158 0.669323928344198 0.142670499085446],...
    'String',S,...
    'FontSize',16,...
    'FitBoxToText','off',...
    'LineWidth',2,...
    'BackgroundColor',[1 1 1]);


angle11 = angle1;
angle11(3) = 0.18;
[tempo01,angle01,omega01,control01,qx,qy] = bspline_interpolation_struct_trad(tt,t1,t2,t3,...
    tempo,angle11,angle2,angle3,N0,N1,N2,t,p,n,diff_times);

diff1 = (angle(1,:)-angle01(1,:));
diff2 = (omega(1,:)-omega01(1,:));
diff3 = (control(1,:)-control01(1,:));

start_diff1 = 0;
end_diff1 = 0;
start_diff2 = 0;
end_diff2 = 0;
start_diff3 = 0;
end_diff3 = 0;
for ind = 1:1000
    if diff1(ind)~=0 && start_diff1 == 0
        start_diff1 = ind-1;
    end
    if diff1(ind) == 0 && start_diff1 ~=0 && end_diff1 == 0
        end_diff1 = ind;
    end
    
    if diff2(ind)~=0 && start_diff2 == 0
        start_diff2 = ind-1;
    end
    if diff2(ind) == 0 && start_diff2 ~=0 && end_diff2 == 0
        end_diff2 = ind;
    end
    
    if diff3(ind)~=0 && start_diff3 == 0
        start_diff3 = ind-1;
    end
    if diff3(ind) == 0 && start_diff3 ~=0 && end_diff3 == 0
        end_diff3 = ind;
    end
end

if end_diff3 == 0
    end_diff3 = 1000;
end

if start_diff3 == 0
    start_diff3 = 1;
end

% Comparison of the results changing only one control point
figure
title('Interpolation Results')
subplot(3,1,1)
hold on
plot(tempo,angle(1,:)); grid on
plot(tempo01,angle01(1,:),'r'); grid on
plot(tempo(start_diff1),angle(1,start_diff1),'o')
plot(tempo(end_diff1),angle(1,end_diff1),'o')

subplot(3,1,2)
hold on
plot(tempo,omega(1,:)); grid on
plot(tempo01,omega01(1,:),'r'); grid on
plot(tempo(start_diff2),omega(1,start_diff2),'o')
plot(tempo(end_diff2),omega(1,end_diff2),'o')

subplot(3,1,3)
hold on
plot(tempo,control(1,:)); grid on
plot(tempo01,control01(1,:),'r'); grid on
plot(tempo(start_diff3),control(1,start_diff3),'o')
plot(tempo(end_diff3),control(1,end_diff3),'o')



%% ALGORITMO TRIANGOLARE

N00 = zeros(L,n);
N11 = zeros(L,n-1);
N22 = zeros(L,n-2);
req_time = zeros(L,1);

%%% SENZA RAFFINAMENTO DEL VETTORE U. IL TEMPO ESCE NEGLI ISTANTI GENERATI
%%% DAL VETTORE U

for ind = 1:L
    
    tic
    [N00(ind,:),N11(ind,:),N22(ind,:)] = triangular_algorithm_TBM_mex(int8(n),int8(p),u(ind),int8(m),t);
    req_time(ind) = toc;
    
end

%%% computational times
disp(' ')
disp('==== ALGORITMO TRIANGOLARE NODI ORIGINALI =====')
disp(' ')
st = sprintf('Total time: %f',sum(req_time));
disp(st);
st = sprintf('Mean time per time instant: %f',mean(req_time));
disp(st);

% CHECK BASIS-FUNCTIONS
err0 = sum(sum(N00-N0));
if err0 ~= 0
    disp('Errore! err0 = %f',err0);
end
err1 = sum(sum(N11-N1));
if err1 ~= 0
    disp('Errore! err1 = %f',err1);
end
err2 = sum(sum(N22-N2));
if err2 ~= 0
    disp('Errore! err2 = %f',err2);
end

if abs(err0) + abs(err1) + abs(err2) == 0
    disp('Algoritmo triangolare sui nodi orginali: OK')
else
    disp('Algoritmo triangolare sui nodi orginali: NOT OK')
end

tempo_int = N00*tt';
angle_int = N00*angle01';

[qx1,qy1,qxx1,qyy1] = derivative_control_points(n,p,t,tt,angle01);

% first derivative
ff_dot = N11*qy1;
gg_dot = N11*qx1;
f_dot = ff_dot./gg_dot;
figure
plot(tempo,omega(1,:))
hold on
plot(tempo_int,f_dot,'r');

% Second derivative
ff_ddot = N22*qyy1;
gg_ddot  = N22*qxx1;
den1 = gg_dot.*gg_dot.*gg_dot;
f_ddot = (ff_ddot.*gg_dot - ff_dot.*gg_ddot)./den1;
figure
plot(tempo,control(1,:))
hold on
plot(tempo_int,f_ddot,'r');

%%% CON RAFFINAMENTO DEL VETTORE U. IL TEMPO ESCE EQUISPAZIATO
%%% DAL VETTORE U

NN00 = zeros(L,n);
NN11 = zeros(L,n-1);
NN22 = zeros(L,n-2);
req_time = zeros(L,1);
Tempo = zeros(L,1);

u1 = u;
for ind = 1:L
    err = 1e10;
    delta = 0;
    U = u1(ind);
    while abs(err) > 1e-12
        U = U - delta;
        tic
        [NN00(ind,:),NN11(ind,:),NN22(ind,:)] = triangular_algorithm_TBM_mex(int8(n),int8(p),U,int8(m),t);
        req_time(ind) = toc;
        Tempo(ind) = NN00(ind,:)*tt';
        err = Tempo(ind) - tempo_ref(ind);
        if delta == 0
            % Inizializzo la procedura di convergenza
            delta = 0.1*err/TF;
        else
            % Arrivo a convergenza con il metodo delle secanti
            delta = (U - U_prev)/(Tempo(ind)-tempo_prev)*err;
        end
        U_prev = U;
        tempo_prev = Tempo(ind);
    end
    u1(ind) = U;
end

disp(' ')
disp('==== ALGORITMO TRIANGOLARE NODI MODIFICATI =====')
disp(' ')
%%% computational times
st = sprintf('Total time: %f',sum(req_time));
disp(st);
st = sprintf('Mean time per time instant: %f',mean(req_time));
disp(st);

tempo_int1 = NN00*tt';
angle_int1 = NN00*angle01';

[qx1,qy1,qxx1,qyy1] = derivative_control_points(n,p,t,tt,angle01);
% first derivative
ff_dot1 = NN11*qy1;
gg_dot1 = NN11*qx1;
f_dot1 = ff_dot1./gg_dot1;

figure
plot(tempo_int,ff_dot)
hold on
plot(tempo_int1,ff_dot1,'r')

figure
plot(tempo_int,gg_dot)
hold on
plot(tempo_int1,gg_dot1,'r')

figure
plot(tempo_int,f_dot)
hold on
plot(tempo_int1,f_dot1,'r')

figure
plot(tempo,omega(1,:))
hold on
plot(tempo_int1,f_dot1,'r');

% Second derivative
ff_ddot1 = NN22*qyy1;
gg_ddot1  = NN22*qxx1;
den1 = gg_dot1.*gg_dot1.*gg_dot1;
f_ddot = (ff_ddot1.*gg_dot1 - ff_dot1.*gg_ddot1)./den1;
figure
plot(tempo,control(1,:))
hold on
plot(tempo_int1,f_ddot,'r');


%%% CHECK BASIS-FUNCTIONS
figure
hold on
plot(u,N00); 
plot(u1,NN00); 
grid on
text_title = sprintf('Basis-function ordine 0');
title(text_title)

figure
hold on
plot(u,N11); 
plot(u1,NN11); 
grid on
text_title = sprintf('Basis-function ordine 0');
title(text_title)

figure
hold on
plot(u,N22); 
plot(u1,NN22); 
grid on
text_title = sprintf('Basis-function ordine 0');
title(text_title)

%%% CHECK INTERPOLATION
figure
hold on
grid on
plot(tempo_int,angle_int)
plot(tempo_int1,angle_int1,'r')

%% PROVE DI INTEGRAZIONE

tt = linspace(0,pi,n);
f1 = exp(tt); %-cos(tt); exp(tt); %(1/3)*tt.^3;
f2 = zeros(size(f1));
f3 = f2;

[tempo,f,f_dot,f_2dot,qx,qy1,f_int] = bspline_interpolation_struct_trad(tt,0,0,0,0,...
    f1,f2,f3,N0,N1,N2,t,p,n,0);

close all

figure
plot(tempo,f(1,:))
figure
plot(tempo,f_dot(1,:))
I = f(1,end)-f(1,1)
Q = trapz(tempo,f_dot(1,:))

[tempo,f,f_dot,f_2dot,qx,qy1] = bspline_interpolation_struct_trad(tt,0,0,0,0,...
    f_int,f2,f3,N0,N1,N2,t,p,n,0);

figure
plot(tempo,f(1,:))
figure
plot(tempo,f_dot(1,:))
I = f(1,end)-f(1,1)
Q = trapz(tempo,f_dot(1,:))
