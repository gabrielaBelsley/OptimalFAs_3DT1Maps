function NumbPixelsAveraged = sizeROI(desiredCoV,achievedCoV)

%   Inputs:
        % desired CoV
        % achieved CoV
        
%   Outputs:
        %NumbPixelsAveraged: number pixels to average to achieve the
        %desired CoV
        
% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022

%desiredCoV = achievedCoV/sqrt(NumbPixelsAveraged);
%<=>NumbPixelsAveraged = (achievedCoV/desiredCoV)^2
NumbPixelsAveraged = ceil((achievedCoV/desiredCoV)^2);

end

