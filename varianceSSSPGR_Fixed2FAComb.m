function [CoVT1_Fixed2FAComb] = varianceSSSPGR_Fixed2FAComb(Fixed2FAComb,CoVT1_2Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst)

% what is the maximum CoVT1 across the T1, B1+ range for a specific FA
% combination of 2 FAs?

%     Inputs: 
%       Fixed2FAComb - Values of 2 FAs to find CoV T1
%       CoVT1_2Angles - CoV T1 vector
%       iTR - TR value
%       istdB1 - B1+ noise value
%       T1Range - range of T1 values
%       B1TrueRange - range of B1+ values
%       SNREst - SNR value

%     Outputs: 
%       CoVT1_Fixed2FAComb - CoV T1 for the 2 FAs input into variable Fixed2FAComb

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


FA1 = Fixed2FAComb(1);
FA2 = Fixed2FAComb(2);


CoVT1_Fixed2FAComb = zeros(1,length(SNREst));
for indexNoise = 1:length(SNREst)
    [~,index_maxstdT1] = max(squeeze(CoVT1_2Angles(FA1,FA2,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxstdT1);
    
    CoVT1_Fixed2FAComb(1,indexNoise) = CoVT1_2Angles(FA1,FA2,indexT1,indexB1,indexNoise,iTR,istdB1);
end

end

