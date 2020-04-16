function [N_0,N_1,N_2] = triangular_algorithm_TBM(n,p,u,m,knot_points)
% Input: n, p, m, u, and m+1 clamped knots { u0, ..., um }
% Output: Coefficients N0,p(u), N1,p(u), ..., Nn,p(u) in N[0], N[1], ..., N[n]

% TO-BE-MEXED VERSION!

% Algorithm

%% degree 0 coefficient

% rule out special cases
if u == knot_points(1)
    % u corrisponde al primo knot point
    k = int8(1);
elseif u == knot_points(m+1)
    % u corrisponde all'ultimo knot point
    k = n+p-1;
else
    % now u is between the first and the last knot point
    % Let u be in knot span [uk,uk+1);
    
    % Original MATLAB code
    % k = int8(sum(u > knot_points));
    
    % C-like code
    k = int8(0);
    while u>knot_points(k+1)
        k = k + 1;
    end
end

max_length = n+p-1;
L_N_1 = 2:n;
L_N_2 = 3:n;
N = zeros(1,max_length);
N_0 = zeros(1,n);    
N_1 = N_0(L_N_1);   % derivata prima
N_2 = N_0(L_N_2);   % derivata seconda

num = 0;            %#ok, for mexing
den = 0;            %#ok, for mexing
coeff = 0;          %#ok, for mexing
coeff1 = 0;         %#ok, for mexing
coeff2 = 0;         %#ok, for mexing

if k == 0
    disp(knot_points)
end
N(k) = 1.0;

%% degree d goes from 1 to p

% globally used quantities
KP_k = knot_points(k);
KP_kplus1 = knot_points(k+1);
num1 = (KP_kplus1 - u);
num2 = (u - KP_k);
zero_num1 = num1 == 0;
zero_num2 = num2 == 0;

for d = 1:p-1
    
    if k - d > 0
        % right (south-west corner) term only
        den = (KP_kplus1 - knot_points(k-d+1));
        
        coeff = num1/den;
        % to deal with singularities
        if zero_num1 && den == 0
            coeff = 1;
        end
        
        N(k-d) = coeff * N((k-d)+1);
        
    end
    
    
    if  d >= 2 && k>=2  && k < n + p - 1
        % compute internal terms
        for i = k-d+1:k-1
            if i > 0 && i + d < m + 1
                coeff1 = (u - knot_points(i))/(knot_points(i+d) - knot_points(i));
                coeff2 = (knot_points(i+d+1) - u)/(knot_points(i+d+1) - knot_points(i+1));
                N(i) = coeff1 * N(i) + coeff2 * N(i+1);
            end
        end
    end
    
    if k <= n + (p-1) - d
        % let (north-west corner) term only
        den = (knot_points(k+d) - KP_k);
        
        coeff = num2/den;
        % to deal with singularities
        if zero_num2 && den == 0
            coeff = 1;
        end
        
        N(k) = coeff * N(k);
        
    end    
    
    if d == p-2
        
        % First derivative
        if N(1) ~= 1 && N(end) ~= 1
            % Common case
            N_1 = N(L_N_1);
        elseif N(1) == 1 
            % to deal with singularities
            N_1 = N(L_N_1 - 1);
        elseif N(end) == 1
            % to deal with singularities
            N_1(end) = 1;
        end
        
    elseif d == p-3
        
        % Second derivative
        if N(1) ~= 1 && N(end) ~= 1
            % Common case
            N_2 = N(L_N_2);
        elseif N(1) == 1
            % to deal with singularities
            N_2 = N(L_N_2 - 2);
        elseif N(end) == 1
            % to deal with singularities
            N_2(end) = 1;
        end
        
    end
    
    
end

% array N[0..n] has the coefficients.
N_0 = N(1:n);

end