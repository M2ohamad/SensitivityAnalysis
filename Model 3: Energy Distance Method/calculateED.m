function [D, Dnormalized] = calculateED(X,Y)
% Calculates the squared energy distance value between two
% independent functions.
%
% Usage:
% D = edist_function(X, Y, N);
%
% Input:
% X = matrix of samples                           - matrix (N, M)
% Y = matrix of samples                           - matrix (N, M)
%                      M = number of parameters   N = sample size
%
% Output:
%       D = square of energy distance                    - scalar
%               D^2(X,Y) = 2A - B - C
%               A = ||X-Y||   B = ||X-X'||  C = ||Y-Y'||
%       Dnormalized = D normalized between 0 and 1       - scalar
%               Dnormalized = D / (2*A)
%               Equals 0 if and only if X and Y are identically distributed

[N1,M] = size(X) ;
[N2,P] = size(Y) ;

if M>2 || P>2; error('Number of parameters must be 2 or less.');
elseif M==1 && P==1
    % Approximate A
    A = 0 ;
    for i = 1:N1
        for j = 1:N2
            dij = abs(X(i,1) - Y(j,1)) ;
            A = A + dij ;
        end
    end
    A = A / N1 / N2 ;
    
    % Approximate B
    B = 0;
    for i = 1:N1
        for j = 1:N1
            dij = abs(X(i,1) - X(j,1)) ;
            B = B + dij ;
        end
    end
    B = B/N1/N1 ;
    
    % Approximate C
    C = 0;
    for i = 1:N2
        for j = 1:N2
            dij = abs(Y(i,1) - Y(j,1)) ;
            C = C + dij ;
        end
    end
    C = C/N2/N2 ;
    
elseif M==2 && P==2
    % Approximate A
    A = 0;
    for i= 1:N1
        for j = 1:N2
            dij = sqrt( ( X(i,1)-Y(j,1) ).^2 + ( X(i,2)-Y(j,2) ).^2 )  ;
            A = A + dij ;
        end
    end
    A = A/N1/N2 ;
    
    % Approximate B
    B = 0;
    for i = 1:N1
        for j = 1:N2
            dij = sqrt( ( X(i,1)-X(j,1) ).^2 + ( X(i,2)-X(j,2) ).^2 )  ;
            B = B + dij ;
        end
    end
    B = B/N1/N1 ;
    
    % Approximate C
    C = 0;
    for i = 1:N1
        for j = 1:N2
            dij = sqrt( ( Y(i,1)-Y(j,1) ).^2 + ( Y(i,2)-Y(j,2) ).^2 )  ;
            C = C + dij ;
        end
    end
    C = C/N2/N2 ;
end

% Calculate Squared Energy Distance
D = 2*A - B - C ;
    
% Calculate Normalized Sqrd Energy Dist
Dnormalized = D / (2*A) ;


