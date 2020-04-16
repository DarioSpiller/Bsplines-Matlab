function N = triangular_algorithm(n,p,u,knot_points)
% Input: n, p, m, u, and m+1 clamped knots { u0, ..., um }
% Output: Coefficients N0,p(u), N1,p(u), ..., Nn,p(u) in N[0], N[1], ..., N[n]

% Algorithm

m = length(knot_points);

%% degree 0 coefficient

% rule out special cases
if u == knot_points(1)
    % u corrisponde al primo knot point
    k = 1;
    N = zeros(1,n);
elseif u == knot_points(end)
    % u corrisponde all'ultimo knot point
    k = n+p-1;
    N = zeros(1,k);
else
    % now u is between the first and the last knot point
    % Let u be in knot span [uk,uk+1);
    k = sum(u > knot_points);
    N = zeros(1,max([k,n]));
end

N(k) = 1.0;

%% degree d goes from 1 to p

for d =1:p-1
    
    if k - d > 0
        % right (south-west corner) term only
        coeff_1 = (knot_points(k+1) - u)/(knot_points(k+1) - knot_points(k-d+1));
        
        % to deal with singularities
        if isnan(coeff_1)
            coeff_1 = 1;
        end
        
        N(k-d) = coeff_1 * N((k-d)+1);
        
    end
    
    
    if  d >= 2 && k>=2  && k < n + p - 1
        % compute internal terms
        for i = k-d+1:k-1
            if i > 0 && i + d < m
                coeff_1 = (u - knot_points(i))/(knot_points(i+d) - knot_points(i));
                coeff_2 = (knot_points(i+d+1) - u)/(knot_points(i+d+1) - knot_points(i+1));
                N(i) = coeff_1 * N(i) + coeff_2 * N(i+1);
            end
        end
    end
    
    if k <= n + (p-1) - d
        % let (north-west corner) term only
        coeff_2 = (u - knot_points(k))/(knot_points(k+d) - knot_points(k));
        
        % to deal with singularities
        if isnan(coeff_2)
            coeff_2 = 1;
        end
        
        N(k) = coeff_2 * N(k);
        
    end
end

% array N[0..n] has the coefficients.
if length(N) > n
    N(n+1:end) = [];
end

end