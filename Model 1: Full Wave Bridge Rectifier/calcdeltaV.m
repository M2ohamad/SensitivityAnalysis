function deltaV = calcdeltaV(R,C)
%This function calculates the change in voltage given a capacitance and
%resistance value.
% Constant Definitions
tau = R * C;
freq = 50;
period = 1/(freq);
w = (2 * pi * freq);
t = 0:0.00001:(period);

% Equation Definitions
v1 =@(t) abs(230 * sqrt(2) * cos(w*t));
t1 = (1/w) * atan(1/(w*tau));
A = 230 * sqrt(2) * cos(w*t1) * exp(t1/tau);
v2 =@(t) A * exp(-t/tau);

% Calculate vMax using t1
vMax = v1(t1);

% Calculate vMin and t2 using Fsolve
myfunc = @(t) A*exp(-t/tau) + sqrt(2)*230*cos(w*t);
opts = optimset('Diagnostics','off', 'Display','off');
t2 = fsolve(myfunc, period/2, opts);
vMin = v2(t2);

deltaV = vMax - vMin;

end

