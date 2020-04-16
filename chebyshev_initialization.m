function T = chebyshev_initialization(order,N,check)

x = linspace(-1,1,N);

T = zeros(order,N);

T(1,:) = ones(size(x));
T(2,:) = x;

for n = 3:order
    T(n,:) = 2*x.*T(n-1,:) - T(n-2,:);
end

if check == 1
    figure
    hold on
    for n = 1:order
        plot(x,T(n,:))
    end
end

end