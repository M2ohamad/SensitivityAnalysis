%% STEP 1: Sample Input 
close all;
% Parameters
N = 10*1e3; % sample size
U = 100; % constant parameter [V]
M = 3; % number of input parameters

amin = 0.2; amax = 0.8;
Lmin = 1e-3; Lmax = 1e-2;
fmin = 5e3; fmax = 5e4;  %fmin = 10e3; fmax = 1e3;

% Values for Sample A
X1_a = amin + (amax-amin).*rand(N,1); % alfa, uniform distribution between 0.2 and 0.8
X2_L = Lmin + (Lmax-Lmin).*rand(N,1); % L, uniform distribution between 1e-3 and 1e-2 [H]
X3_f = fmin + (fmax-fmin).*rand(N,1); % f, uniform distribution between 1e3 and 10e3 [Hz]
sampA = [X1_a, X2_L, X3_f];

% Values for Sample B
X1_a = amin + (amax-amin).*rand(N,1);
X2_L = Lmin + (Lmax-Lmin).*rand(N,1);
X3_f = fmin + (fmax-fmin).*rand(N,1);
sampB = [X1_a, X2_L, X3_f];

% Constructing C1, C2, C3;
sampC1 = [sampA(:,1), sampB(:,2), sampB(:,3)];
sampC2 = [sampB(:,1), sampA(:,2), sampB(:,3)];
sampC3 = [sampB(:,1), sampB(:,2), sampA(:,3)];

%% STEP 2: Evaluating Model Against Sampled Data
D_AB123 = zeros(N,5);
D_AB123(:,1) = modelfunction(U,sampA);
D_AB123(:,2) = modelfunction(U,sampB);
D_AB123(:,3) = modelfunction(U,sampC1);
D_AB123(:,4) = modelfunction(U,sampC2);
D_AB123(:,5) = modelfunction(U,sampC3);
%outputs for inputs A, B, C1, C2, and C3 respectively in [Amps]

%% STEP 3: Post-processing
% Total (Unconditional) Variance: Total Average Deviation
E_D(1) = (sum(D_AB123(:,1)) + sum(D_AB123(:,2))) / (2 * N);
E_D(2) = ((sum(D_AB123(:,1).^2)) + (sum(D_AB123(:,2).^2))) / (2 * N);
TotalVar = E_D(2) - (E_D(1))^2;

% Conditional Variance
% A & C1, keeping X1 fixed
E_X1(1) = (sum(D_AB123(:,1)) + sum(D_AB123(:,3))) / (2 * N);   % f0
E_X1(2) = (sum( D_AB123(:,1).*D_AB123(:,3) )) / N;
Var(1) = E_X1(2) - (E_X1(1))^2;

% A & C2, keeping X2 fixed
E_X2(1) = (sum(D_AB123(:,1)) + sum(D_AB123(:,4))) / (2 * N);
E_X2(2) = (sum( D_AB123(:,1).*D_AB123(:,4) )) / N;
Var(2) = E_X2(2) - (E_X2(1))^2;

% A & C3, keeping X3 fixed
E_X3(1) = (sum(D_AB123(:,1)) + sum(D_AB123(:,5))) / (2 * N);
E_X3(2) = (sum( D_AB123(:,1).*D_AB123(:,5) )) / N;
Var(3) = E_X3(2) - (E_X3(1))^2;

% First Order Sensitivity Indices
Sindices(1) = Var(1) / TotalVar;   % Signficance of Alfa
Sindices(2) = Var(2) / TotalVar;   % Sensitivity of Inductance [H]
Sindices(3) = Var(3) / TotalVar;   % Sensitivity of frequency [Hz]
IndicesSum = sum(Sindices);

% Non-additive model (interaction between parameters exists),
% therefore we expect the sum of first order indices to be less than one.

%% STEP 4: Visualization
% Plotting Sensitivity Indicies vs Sample Sizes
figure;
labelinputs = {'S1: alfa', 'S2: L', 'S3: f','S1+S2+S3'};
k = N/10:N/10:N;
SampledIndices = zeros(length(k),M+1);  % M+1 is including sum of indices
for i = 1:length(k)
    SampledIndices(i,:) = converge(D_AB123,k(i));
end
plot_convergence(SampledIndices(:,1:4),k,[],[],[],[],[],labelinputs);

% Histogram of Output A
VA = D_AB123(:,1);
figure; hist(VA); title('Current Output from Sample A'); 
xlabel('Current [A]'); ylabel('Number of Points');



% Scatter Plot of Input vs Output for Sample A
figure; scatter(sampA(:,1),D_AB123(:,1),3); 
title('Alfa vs Output for Sample A'); xlabel('Alfa'); ylabel('Current [A]');
figure; scatter(sampA(:,2),D_AB123(:,1),3); 
title('L vs Output for Sample A'); xlabel('Inductance [H]'); ylabel('Current [A]');
figure; scatter(sampA(:,3),D_AB123(:,1),3); 
title('Freq vs Output for Sample A'); xlabel('Frequency [Hz]'); ylabel('Current [A]');
hold on


%% Regional Sensitivity Analysis 
%Ripple of current output should not be more that 0.5 A
X = sampA;
Y = D_AB123(:,1);
threshold = 0.5;

% Finding Behavioral Samples
[N,M] = size(X);
[N,P] = size(Y);

% Creates a matrix of size (N,1) of the threshold value. Then logically
% compares if Y is less than that value
behav_indx = (Y<repmat(threshold,N,1))==P;
mvd = nan(1,M) ;  % preallocation of max vertical distance

% Factor Mapping
Xb = X(behav_indx,:);    % Collects all the behavioral inputs in X
Nb = size(Xb,1);         % Number of behavioral inputs
Xnb = X(~behav_indx,:);  % Collects all non-behavioral inputs in X
Nnb = size(Xnb,1);       % Number of non-behavioral inputs

% Computing Indices
for i=1:M
   xx = unique(sort(X(:,1)));
   CDFb = empiricalcdf(Xb(:,1),xx);
   CDFnb = empiricalcdf(Xnb(:,1),xx);
   mvd(i) = max(abs(CDFb-CDFnb));
end

plot(xx,CDFb,'Color','r','LineWidth',2)
hold on
plot(xx,CDFnb,'Color','b','LineWidth',2)




