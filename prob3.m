clear all;
x = 2.^[0 1 2 3 4 5 6 7 8];
N = length(x);
n = zeros(1,N);
J0x = zeros(1,N);
tol = 10^-12;
for i = 1:1:N
    ii = 2;
    h = 2*pi/ii;
    t = 0:h:2*pi;
    %Bessel function, change to periodic form
    y = 1/(2*pi)*cos(x(i)*cos(t));
    J0 = trapz(t,y);
    while(1)
        ii = ii+1;
        h = 2*pi/ii;
        t = 0:h:2*pi;
        y = 1/(2*pi)*cos(x(i)*cos(t));
        if(abs(J0-trapz(t,y))<tol)
            break;
        end
        J0 = trapz(t,y);
    end
    n(i) = ii;
    J0x(i) = J0;
end
%output
disp('  x            J0(x)        n');
disp('--------------------------------');
for i=1:1:N
    fprintf('%3d  ',x(i));
    fprintf('%20.16f  ', J0x(i));
    fprintf('%3d\n', n(i));
end
disp('--------------------------------');
