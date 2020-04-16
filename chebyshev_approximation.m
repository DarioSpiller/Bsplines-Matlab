function f_approx = chebyshev_approximation(T,coeff,check)

if size(coeff,2) ~= size(T,1)
    
    coeff = coeff';

end

f_approx = coeff*T;

if check == 1
    figure
    plot(x,f_approx)
end

end