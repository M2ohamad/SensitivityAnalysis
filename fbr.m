tau = 500 * 50 * 10^(-6);
t = 0:0.001:(tau);
v1 = 230*sqrt(2)*cos(100*pi*t);
t1 = arctan(4/(10*pi))*(1/(100*pi));
A = sqrt(2)*230*cos(100*pi*t1)*exp(t1/tau);
v2 = A*exp(-t/tau);