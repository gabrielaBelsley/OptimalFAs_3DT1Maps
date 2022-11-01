function [CoVT1_Fixed4FAComb] = varianceSSSPGR_Fixed4FAComb(Fixed4FAComb,CoVT1_4Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst)

% what is the maximum CoVT1 across the T1, B1+ range for a specific FA
% combination of 4 FAs?

%     Inputs: 
%       Fixed4FAComb - Values of 4 FAs to find CoV T1
%       CoVT1_4Angles - CoV T1 vector
%       iTR - TR value
%       istdB1 - B1+ noise value
%       T1Range - range of T1 values
%       B1TrueRange - range of B1+ values
%       SNREst - SNR value

%     Outputs: 
%       CoVT1_Fixed4FAComb - CoV T1 for the 4 FAs input into variable Fixed4FAComb

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


FA1 = Fixed4FAComb(1);
FA2 = Fixed4FAComb(2);
FA3 = Fixed4FAComb(3);
FA4 = Fixed4FAComb(4);

CoVT1_Fixed4FAComb = zeros(1,length(SNREst));
for indexNoise = 1:length(SNREst)
    [~,index_maxstdT1] = max(squeeze(CoVT1_4Angles(FA1,FA2,FA3,FA4,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxstdT1);
    
    CoVT1_Fixed4FAComb(1,indexNoise) = CoVT1_4Angles(FA1,FA2,FA3,FA4,indexT1,indexB1,indexNoise,iTR,istdB1);
end

end

