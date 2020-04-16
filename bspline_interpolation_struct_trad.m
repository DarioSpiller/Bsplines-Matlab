function [tempo,f,f_dot,f_2dot,qx,qy1,f_int] = bspline_interpolation_struct_trad(tt,t1,t2,t3,tf,...
    f1,f2,f3,N0,N1,N2,t,k,m,diff_times)

% Questo file è il sorgente per il .mex

nn = m-1;
f = zeros(3,size(N0,1));
f_dot = zeros(3,size(N0,1));
f_2dot = zeros(3,size(N0,1));
Tempo = 0;
Tempo1 = 0;
Tempo2 = 0;
Tempo3 = 0;
qx1 = zeros(nn,1);
qx2 = zeros(nn,1);
qx3 = zeros(nn,1);
qx = zeros(nn,1);
g_dot1  = zeros(nn,1);
g_dot2  = zeros(nn,1);
g_dot3  = zeros(nn,1);
g_dot  = zeros(nn,1);
g_2dot1  = zeros(nn,1);
g_2dot2  = zeros(nn,1);
g_2dot3  = zeros(nn,1);
g_2dot  = zeros(nn,1);

%%% COSTRUZIONE DELLE CURVE CON LE B-SPLINES
if diff_times == 1
    Tempo1 = t1*tf;
    Tempo2 = t2*tf;
    Tempo3 = t3*tf;
else
    Tempo = tt;
end

f(1,:) = N0*f1';
f(2,:) = N0*f2';
f(3,:) = N0*f3';
if diff_times == 1
    tempo1    = N0*Tempo1';
    tempo2    = N0*Tempo2';
    tempo3    = N0*Tempo3';
    tempo = [tempo1 tempo2 tempo3];
else
    tempo    = N0*Tempo';
end

%%% Derivata prima

if isempty(N1)
    % Per calcolare la derivata prima devo aver combinato almeno due
    % polinomi precedentemente (ovvero f1, f2 ed f3 devono avere almeno due
    % componenti)
    warning('Non posso calcolare la derivata prima!')
    f_dot = NaN;
else
    qy1 = zeros(nn,1);
    qy2 = zeros(nn,1);
    qy3 = zeros(nn,1);
    if diff_times == 1
        qx1 = zeros(nn,1);
        qx2 = zeros(nn,1);
        qx3 = zeros(nn,1);
    else
        qx = zeros(nn,1);
    end
    
    tt = t(2:end-1);
    
    for ii = 1:nn
        if diff_times == 1
            qx1(ii) = (k-1)*(Tempo1(ii+1) - Tempo1(ii))/(tt(k+ii-1)-tt(ii));
            qx2(ii) = (k-1)*(Tempo2(ii+1) - Tempo2(ii))/(tt(k+ii-1)-tt(ii));
            qx3(ii) = (k-1)*(Tempo3(ii+1) - Tempo3(ii))/(tt(k+ii-1)-tt(ii));
        else
            qx(ii) = (k-1)*(Tempo(ii+1) - Tempo(ii))/(tt(k+ii-1)-tt(ii));
        end
        qy1(ii) = (k-1)*(f1(ii+1) - f1(ii))/(tt(k+ii-1)-tt(ii));
        qy2(ii) = (k-1)*(f2(ii+1) - f2(ii))/(tt(k+ii-1)-tt(ii));
        qy3(ii) = (k-1)*(f3(ii+1) - f3(ii))/(tt(k+ii-1)-tt(ii));
        
    end
    
    f_int(1) = 0;
    for ii = 1:nn
       f_int(ii+1) = f_int(ii) + qy1(ii)*(tt(k+ii-1)-tt(ii))/(k-1);
    end
    
    
    
    f1_dot = N1*qy1;
    f2_dot = N1*qy2;
    f3_dot = N1*qy3;
    if diff_times == 1
        g_dot1  = N1*qx1;
        g_dot2  = N1*qx2;
        g_dot3  = N1*qx3;
    else
        g_dot  = N1*qx;
    end
    
    if diff_times == 1
        f_dot(1,:) = f1_dot./g_dot1;
        f_dot(2,:) = f2_dot./g_dot2;
        f_dot(3,:) = f3_dot./g_dot3;
    else
        f_dot(1,:) = f1_dot./g_dot;
        f_dot(2,:) = f2_dot./g_dot;
        f_dot(3,:) = f3_dot./g_dot;
    end
end

%%% Derivata seconda
if isempty(N2)
    % Per calcolare la derivata seconda devo aver combinato almeno tre
    % polinomi precedentemente (ovvero f1, f2 ed f3 devono avere almeno tre
    % componenti)
    warning('Non posso calcolare la derivata seconda!')
    f_2dot = NaN;
else
    qyy1 = zeros(nn-1,1);
    qyy2 = zeros(nn-1,1);
    qyy3 = zeros(nn-1,1);
    qxx1 = zeros(nn-1,1);
    qxx2 = zeros(nn-1,1);
    qxx3 = zeros(nn-1,1);
    qxx = zeros(nn-1,1);
    
    ttt = t(3:end-2);
    
    for ii = 1:nn-1
        if diff_times == 1
            qxx1(ii)  = (k-2)*(qx1(ii+1) - qx1(ii))/(ttt(k+ii-2)-ttt(ii));
            qxx2(ii)  = (k-2)*(qx2(ii+1) - qx2(ii))/(ttt(k+ii-2)-ttt(ii));
            qxx3(ii)  = (k-2)*(qx3(ii+1) - qx3(ii))/(ttt(k+ii-2)-ttt(ii));
        else
            qxx(ii)  = (k-2)*(qx(ii+1) - qx(ii))/(ttt(k+ii-2)-ttt(ii));
        end
        
        qyy1(ii) = (k-2)*(qy1(ii+1) - qy1(ii))/(ttt(k+ii-2)-ttt(ii));
        qyy2(ii) = (k-2)*(qy2(ii+1) - qy2(ii))/(ttt(k+ii-2)-ttt(ii));
        qyy3(ii) = (k-2)*(qy3(ii+1) - qy3(ii))/(ttt(k+ii-2)-ttt(ii));
    end
    
    f1_2dot = N2*qyy1;
    f2_2dot = N2*qyy2;
    f3_2dot = N2*qyy3;
    
    if diff_times == 1
        g_2dot1  = N2*qxx1;
        g_2dot2  = N2*qxx2;
        g_2dot3  = N2*qxx3;
    else
        g_2dot  = N2*qxx;
    end
    
    if diff_times == 1
        den1 = g_dot1.*g_dot1.*g_dot1;
        den2 = g_dot2.*g_dot2.*g_dot2;
        den3 = g_dot3.*g_dot3.*g_dot3;
        
        f_2dot(1,:) =(f1_2dot.*g_dot1-f1_dot.*g_2dot1)./den1;
        f_2dot(2,:) =(f2_2dot.*g_dot2-f2_dot.*g_2dot2)./den2;
        f_2dot(3,:) =(f3_2dot.*g_dot3-f3_dot.*g_2dot3)./den3;
    else
        den = g_dot.*g_dot.*g_dot;
        f_2dot(1,:) =(f1_2dot.*g_dot-f1_dot.*g_2dot)./den;
        f_2dot(2,:) =(f2_2dot.*g_dot-f2_dot.*g_2dot)./den;
        f_2dot(3,:) =(f3_2dot.*g_dot-f3_dot.*g_2dot)./den;
    end
    
end
end