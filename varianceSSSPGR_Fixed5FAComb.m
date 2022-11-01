function [CoVT1_Fixed5FAComb] = varianceSSSPGR_Fixed5FAComb(Fixed5FAComb,CoVT1_5Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst)

% what is the maximum CoVT1 across the T1, B1+ range for a specific FA
% combination of 5 FAs?

%     Inputs: 
%       Fixed5FAComb - Values of 5 FAs to find CoV T1
%       CoVT1_5Angles - CoV T1 vector
%       iTR - TR value
%       istdB1 - B1+ noise value
%       T1Range - range of T1 values
%       B1TrueRange - range of B1+ values
%       SNREst - SNR value

%     Outputs: 
%       CoVT1_Fixed5FAComb - CoV T1 for the 5 FAs input into variable Fixed5FAComb


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


FA1 = Fixed5FAComb(1);
FA2 = Fixed5FAComb(2);
FA3 = Fixed5FAComb(3);
FA4 = Fixed5FAComb(4);
FA5 = Fixed5FAComb(5);

CoVT1_Fixed5FAComb = zeros(1,length(SNREst));
for indexNoise = 1:length(SNREst)
    [~,index_maxstdT1] = max(squeeze(CoVT1_5Angles(FA1,FA2,FA3,FA4,FA5,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxstdT1);
    
    CoVT1_Fixed5FAComb(1,indexNoise) = CoVT1_5Angles(FA1,FA2,FA3,FA4,FA5,indexT1,indexB1,indexNoise,iTR,istdB1);
end

end

