function [N,t,u] = define_bsplines(L,k,m,t)


N(1:k,1) = struct('valori',[]);

% Computazione

if nargin == 3
    %%% Knot points adimensionalizzati (tra 0 e 1)
    t = 1:1:m+k;
    
    t(t<=k) = 0;
    t(t>k & t<=m) = (1:m-k)/(m-k+1);
    t(t>m) = 1;
    
end

%%% Base functions
u = 0:t(end)/(L-1):t(end);

%% Questo modulo serve a evitare singolarità, che avvengono quando un 
% elemento di u è uguale a uno di t
err = 0;
for ind = 2:m+k-1
    err = err + sum(u(2:end-1)==t(ind));
end

if err ~= 0
    u(2:end-1) = u(2:end-1) + 1e-14;
end

%% Valutazione di N

N(1).valori = zeros(L,k+m-1);

for i = 1:m+k-1
    
    N(1).valori((u>= t(i) & u<=t(i+1)),i) = 1;
    
end

valori = zeros(L,k+m-2);

for kk = 2:k
    
    y = N(kk-1).valori;
        
    for i = 1:m+k-kk
        
        A = (u(:)-t(i))./(t(i+kk-1)-t(i));
        B = (t(i+kk)-u(:))./(t(i+kk)-t(i+1));
        
        A(isnan(A)) = 0;
        B(isnan(B)) = 0;
        
        A(isinf(A)) = 0; 
        B(isinf(B)) = 0; 
        
        valori(:,i) = A.*y(:,i) + B.*y(:,i+1);
        
    end
    
    N(kk).valori = valori(:,1:k+m-kk);
    
end

end
