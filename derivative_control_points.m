function [qx1,qy1,qxx1,qyy1] = derivative_control_points(N,P,t,t1,f1)
%%% Derivata prima
qy1 = zeros(N-1,1);
qx1 = zeros(N-1,1);

% Definisco il nuovo vettore tt togliendo primo e ultimo elemento
tt = t(2:end-1);

for ii = 1:N-1
    
    % aggiorno il valore di p (ricorda che il grado del polinomio prima di
    % derivare è P-1. Il grado della derivata prima è quindi P-2.
    p = P-1;
    
    den = (tt(p+ii)-tt(ii));
    
    % La p davanti è simile al coefficiente della derivata del polinomio:
    % f = x^n --> f' = n*x^(n-1)
    qx1(ii) = p*(t1(ii+1) - t1(ii))/den;
    qy1(ii) = p*(f1(ii+1) - f1(ii))/den;
    
end

%%% Derivata seconda

qyy1 = zeros(N-2,1);
qxx1 = zeros(N-2,1);

% Definisco il nuovo vettore tt togliendo primo e ultimo elemento
ttt = t(3:end-2);

for ii = 1:N-2
    
    % aggiorno il valore di p (ricorda che il grado del polinomio prima di
    % derivare è P-1. Il grado della derivata seconda è quindi P-3.
    p = P-2;
    
    den = (ttt(p+ii)-ttt(ii));
    
    % La p davanti è simile al coefficiente della derivata del polinomio:
    % f = x^n --> f' = n*x^(n-1) --> f'' = n*(n-1)*x^(n-2)
    qxx1(ii) = p*(qx1(ii+1) - qx1(ii))/den;    
    qyy1(ii) = p*(qy1(ii+1) - qy1(ii))/den;
    
end

end