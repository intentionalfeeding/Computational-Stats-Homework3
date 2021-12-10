%MAGY6973 Computational Statistics Homework 3 Problem 4(b)
%The subroutine is to sampling from the beta distribution using inverse CDF method
function x = Sampling_beta(N,alpha,beta)
x = zeros(N,1);
u = rand(N,1);
constant = gamma(alpha+beta)/gamma(alpha)/gamma(beta);
for i =1:N
    y = u(i);
    t = 0.5;
    CDF = beta_CDF(t,alpha,beta);
    PDF = constant*t^(alpha-1)*(1-t)^(beta-1); 
    error = abs(CDF-y);
    while error > 10^-2
    t = t-(CDF-y)/PDF;
    CDF = beta_CDF(t,alpha,beta);
    PDF = gamma(alpha+beta)/gamma(alpha)/gamma(beta)*t^(alpha-1)*(1-t)^(beta-1);
    error = abs(CDF-y);
    end
    x(i) = t;
end
end




