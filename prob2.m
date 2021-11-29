clear all;
syms X;
f(X) = cos(X);
n = 2.^[1 2 3 4];
N = length(n);
errI = zeros(1,N);
errP = zeros(1,N);
for i=1:1:N
    h = pi/(n(i)-1);
    x = 0:h:pi;
    y = cos(x);
    p = polyfit(x,y,n(i)-1);
    
    %by degree n-1 polynomial:
    PP(X) = 0*X;
    for j=1:1:n(i)
        PP(X) = PP(X)+p(j)*X^(n(i)-j);
    end
    errI(i) = sqrt(double(int((PP(X)-f(X))^2,X,0,pi)));
    
    %by projecting onto nodes:
    P = zeros(n(i),n(i));
    for j = 1:1:n(i)
        for k = 1:1:n(i)
            P(j,k) = legendreP(j-1,x(k)/pi);
        end
    end
    C = inv(P')*y';
    f_hat(X) = 0*X;
    for j=1:1:n(i)
        f_hat(X) = f_hat(X)+C(j)*legendreP(j-1,X/pi);
    end
    errP(i) = sqrt(double(int((f_hat(X)-f(X))^2,X,0,pi)));
    %plot
    xx = 0:0.001:pi;
    figure(i);
    plot(xx,f(xx));
    hold on;
    plot(xx,PP(xx));
    hold on;
    plot(xx,f_hat(xx));
    title(['n = ',num2str(n(i))]);
    xlabel('x');
    ylabel('y');
    legend('f(x)=cos(x)','polynomial','legendre');
end
disp('  n         errI                errP');
disp('---------------------------------------------');
for i=1:1:N
    fprintf('%3d',n(i));
    fprintf('%20.16f',errI(i));
    fprintf('%20.16f\n',errP(i));
end
disp('---------------------------------------------');