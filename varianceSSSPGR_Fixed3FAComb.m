function [CoVT1_Fixed3FAComb] = varianceSSSPGR_Fixed3FAComb(Fixed3FAComb,CoVT1_3Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst)

% what is the maximum CoVT1 across the T1, B1+ range for a specific FA
% combination of 3 FAs?

%     Inputs: 
%       Fixed3FAComb - Values of 3 FAs to find CoV T1
%       CoVT1_3Angles - CoV T1 vector
%       iTR - TR value
%       istdB1 - B1+ noise value
%       T1Range - range of T1 values
%       B1TrueRange - range of B1+ values
%       SNREst - SNR value

%     Outputs: 
%       CoVT1_Fixed3FAComb - CoV T1 for the 3 FAs input into variable Fixed3FAComb

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


FA1 = Fixed3FAComb(1);
FA2 = Fixed3FAComb(2);
FA3 = Fixed3FAComb(3);


CoVT1_Fixed3FAComb = zeros(1,length(SNREst));
for indexNoise = 1:length(SNREst)
    [~,index_maxstdT1] = max(squeeze(CoVT1_3Angles(FA1,FA2,FA3,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxstdT1);
    
    CoVT1_Fixed3FAComb(1,indexNoise) = CoVT1_3Angles(FA1,FA2,FA3,indexT1,indexB1,indexNoise,iTR,istdB1);
end

end

