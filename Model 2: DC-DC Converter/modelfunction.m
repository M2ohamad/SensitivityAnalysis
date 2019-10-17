function DX = modelfunction(U,X)
%MODELFUNCTION D = alfa * (1-alfa) * U / L / F;   
%   Outputs the current of the DC-DC converter [A]
DX =  X(:,1) .* (1 - X(:,1)) .* U ./ X(:,2) ./ X(:,3);
end

