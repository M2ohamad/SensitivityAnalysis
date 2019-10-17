function [totalInfluence,ED_intervals] = calc_totalInfluence(X_LFI, Y , nIntervals, intervalSize)
% Calculates the total influence of a parameter by summing the weighted
% energy distances of each interval in the input space
%
% Usage:
% [variable_Influence, variable_ED_intervals] = calc_totalInfluence(X, Y, nInt, intSize);
%
% Input:
% X = matrix of samples                           - matrix (N, M)
% Y = matrix of samples                           - matrix (N, M)
%                      M = number of parameters   N = sample size
% nInt = number of intervals                      - scalar
% intSize = size of each interval                 - scalar
%
% Output:
%     totalInfluence = summation of weighted E.D. for each interval - scalar
%     ED_intervals = row vector of E.D. for each interval   - matrix (1, nInt)
%
% Use reshape to extract corresponding Y points in each interval
% Each page in multidimensional array is an interval
L = reshape( X_LFI(:,1) , intervalSize , 1 , nIntervals ) ;
F = reshape( X_LFI(:,2) , intervalSize , 1 , nIntervals ) ;
I = reshape( X_LFI(:,3) , intervalSize , 1 , nIntervals ) ;
%LF = cat( 2 , L , F ) ;

% Calcualte the energy distance for each interval against the total sample space
ED_intervals = 1:nIntervals ;
for i = 1:nIntervals
    ED_intervals(i) = calculateED( Y , I(:,:,i) ) ;
end
% Calculate the probabilty that X lies in each interval (same for unif
% dist)
probability_int = intervalSize / length(X_LFI) ;

% Calculate the total influence
totalInfluence = 0 ;
for i = 1:nIntervals
    infl = probability_int * ED_intervals(i) ;
    totalInfluence = totalInfluence + infl ;
end

end

