function DX = modelfunction(X,alpha,U)
%MODELFUNCTION D = alfa * (1-alfa) * U / L / F;   
%   Outputs the current of the DC-DC converter [A]
DX =  alpha .* (1 - alpha) .* U ./ X(:,2) ./ X(:,1);
end

