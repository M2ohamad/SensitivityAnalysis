function output = converge(Y,Nk)
% Repeat computations using an increasing number of samples so to assess
% convergence: 
% Method: Take a large sample (50,000) and randomly select an increasing
% fixed size to use as actual sample.
sY = datasample(Y,Nk);

% Total (Unconditional) Variance: Total Average Deviation
E_D(1) = (sum(sY(:,1)) + sum(sY(:,2))) / (2 * Nk);
E_D(2) = ((sum(sY(:,1).^2)) + (sum(sY(:,2).^2))) / (2 * Nk);
TotalVar = E_D(2) - (E_D(1))^2;

% Conditional Variance
% A & C1, keeping X1 fixed
E_X1(1) = (sum(sY(:,1)) + sum(sY(:,3))) / (2 * Nk);   % f0
E_X1(2) = (sum( sY(:,1).*sY(:,3) )) / Nk;
Var(1) = E_X1(2) - (E_X1(1))^2; 

% A & C2, keeping X2 fixed
E_X2(1) = (sum(sY(:,1)) + sum(sY(:,4))) / (2 * Nk);
E_X2(2) = (sum( sY(:,1).*sY(:,4) )) / Nk;
Var(2) = E_X2(2) - (E_X2(1))^2;

% A & C3, keeping X3 fixed 
E_X3(1) = (sum(sY(:,1)) + sum(sY(:,5))) / (2 * Nk);
E_X3(2) = (sum( sY(:,1).*sY(:,5) )) / Nk;
Var(3) = E_X3(2) - (E_X3(1))^2;

% First Order Sensitivity Indices
S1 = Var(1) / TotalVar;   % Signficance of Alfa
S2 = Var(2) / TotalVar;   % Sensitivity of Inductance [H]
S3 = Var(3) / TotalVar;   % Sensitivity of frequency [Hz]
Sum = S1 + S2 + S3;

output = [S1 S2 S3 Sum];
   
end

