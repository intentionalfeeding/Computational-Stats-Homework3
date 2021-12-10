%MAGY6973 Computational Statistics Homework 3 Problem 4(a)
%The subroutine is to compute the beta CDF using adaptive 8th-oder Gauss-Legendre quadrature
%12-03-2021
function CDF = beta_CDF(x,alpha,beta)
%compute the roots and weights of 8th order of legendre polynomial
n = 8;
[y_roots,w_values] = get_legendre(n);
%
%adaptive procedures
I = myfunction(0,x,y_roots,w_values,alpha,beta,100);
CDF = I*gamma(alpha+beta)/gamma(alpha)/gamma(beta);
%%%%%        
    function I = myfunction(a_top,b_top,y_roots,w_values,alpha,beta,I_top)
    a = a_top; b = a_top + (b_top-a_top)/2;    
    f_values = (b-a)/2*((b-a)/2*(y_roots+1)+a).^(alpha-1).*(1-(b-a)/2*(y_roots+1)-a).^(beta-1);
    I_bottom_1 = sum(f_values.*w_values);
    a = a_top + (b_top-a_top)/2 ; b = b_top;   
    f_values = (b-a)/2*((b-a)/2*(y_roots+1)+a).^(alpha-1).*(1-(b-a)/2*(y_roots+1)-a).^(beta-1);
    I_bottom_2 = sum(f_values.*w_values);
    %
    error = abs(I_top - I_bottom_1 - I_bottom_2);
    if error > 10^-8
       I1 = myfunction(a_top,a_top + (b_top-a_top)/2,y_roots,w_values,alpha,beta,I_bottom_1);
       I2 = myfunction(a_top + (b_top-a_top)/2,b_top,y_roots,w_values,alpha,beta,I_bottom_2);
       I = I1+I2;
    else
       I = I_top; 
    end
    end
%%%%%
end




