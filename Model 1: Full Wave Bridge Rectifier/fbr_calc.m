%% STEP 1: Define the Model: Additive Model
close all; clc;
% Constant Definitions
N = 1000; % sample size
R = 50; % R:50 - 200
C = 1 * 10^(-4); %C: 10^-4  -  5*10^-4
tau = R * C;
freq = 50;
period = 1/(2*freq);
w = (2 * pi * freq);
t = 0:0.00001:(period);

% Model Equation Definitions
v1 =@(t) abs(230 * sqrt(2) * cos(w*t));
t1 = (1/w) * atan(1/(w*tau));
A = 230 * sqrt(2) * cos(w*t1) * exp(t1/tau);
v2 =@(t) A * exp(-t/tau);

% Calculate Maximum Voltage
vMax = 230*sqrt(2);

% Calculate Minimum Voltage
myfunc = @(t2) A*exp(-t2/tau) + sqrt(2)*230*cos(w*t2);
opts = optimset('Diagnostics','off', 'Display','off');
t2 = fsolve(myfunc, period/2, opts);
vMin = v2(t2);

% Calculate deltaV
deltaV = vMax - vMin;

%% STEP 2: Sample Input
% Generate Uniform Random Sample Matricies A and B
% r = a + (b-a).*rand(N,1);
CA = 0.0001 + (0.0005-0.0001).*rand(N,1);
CB = 0.0001 + (0.0005-0.0001).*rand(N,1);

RA = 50 + (200-50).*rand(N,1);
RB = 50 + (200-50).*rand(N,1);

arrA = [CA,RA];
arrB = [CB,RB];

% Construct C1 and C2
C_1 = [CA, RB];
C_2 = [CB, RA];

% Define the deltaV arrays all at once using the deal command. Each array
% is a column vector of N values.
[dV_A, dV_B, dV_C1, dV_C2] = deal((1:N)');

%% STEP 3: Evaluating Model Against Sampled Input Combinations
for rows = 1:N
    for col = 1:2
        Ri = arrA(rows,2);
        Ci = arrA(rows,1);
        dV_A(rows) = calcdeltaV(Ri,Ci);
        
        Ri = arrB(rows,2);
        Ci = arrB(rows,1);
        dV_B(rows) = calcdeltaV(Ri,Ci);
        
        Ri = C_1(rows,2);
        Ci = C_1(rows,1);
        dV_C1(rows) = calcdeltaV(Ri,Ci);
        
        Ri = C_2(rows,2);
        Ci = C_2(rows,1);
        dV_C2(rows) = calcdeltaV(Ri,Ci);
    end
end

%% STEP 4: Post-processing input/output samples to compute Sensitivity Indices

% Total (Unconditional) Variance: Total Average Deviation
E_dV = (sum(dV_A) + sum(dV_B)) / (2 * N);
E2_dV = ((sum(dV_A.^2)) + (sum(dV_B.^2))) / (2 * N);
TotalVar = E2_dV - (E_dV)^2;

% Conditional Variance
% X1 and C1
E_X1 = (sum(dV_C1) + sum(dV_A)) / (2 * N);
E2_X1 = (sum(dV_C1 .* dV_A)) / N;
Var_X1 = E2_X1 - (E_X1)^2;
% X2 and C2
E_X2 = (sum(dV_C2) + sum(dV_A)) / (2 * N);
E2_X2 = (sum(dV_C2 .* dV_A)) / N;
Var_X2 = E2_X2 - (E_X2)^2;

% First-Order Sobol Indices  (Main Effect Index)
S1 = Var_X1 / TotalVar;
S2 = Var_X2 / TotalVar;
SSum = S1+S2;

% Used to collect indicies vs sample size data
% j = j + 1;
% for i=j
%     ssnn(i,:) = [SSum, S1, S2, N];
% end

%% STEP 5: Visualization
% Plot of Model
figure (1);
plot(t,v1(t),'linewidth',2)
set(gca,'xlim',[min(t) max(t)],'ylim',[0 350])
hold on
plot(t,v2(t),'g','linewidth',2)
plot(t2,vMin,'r.','markersize',18)
plot(t1,v1(t1),'r.','markersize',18)
plot(0,vMax,'r.','markersize',18)

% Plotting Histograms of Delta V
figure (2);
subplot(4,1,1), histogram(dV_A), title('(A) Change in Voltage');
subplot(4,1,2), histogram(dV_B), title('(B) Change in Voltage');
subplot(4,1,3), histogram(dV_C1), title('(C1) Change in Voltage');
subplot(4,1,4), histogram(dV_C2), title('(C2) Change in Voltage');

% Plotting Sensitivity Indicies vs Sample Sizes
% figure (3);
% %load ssnn; % Loads the indicies vs sample data  (SS contains (S1+S2, S1,
% %S2) data) % NEW data must be saved first!
% labelinputs = {'S1+S2', 'S1', 'S2'};
% plot_convergence(ssnn(:,1:3),ssnn(:,4)',[],[],[],[],[],labelinputs);

% Plot of DeltaV vs Input Parameters
figure (4);
scatter(arrA(:,1),dV_A); 
xlabel('Sample A Capacitance (F)'),ylabel('vMax - vMin (V)');

figure (5);
scatter(arrA(:,2),dV_A);
xlabel('Sample A Resistance (Ohm)'),ylabel('vMax - vMin (V)');


