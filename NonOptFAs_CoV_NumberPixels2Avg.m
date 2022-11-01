function [achievedCoV_Standard5FAs,achievedCoV_Standard4FAs,achievedCoV_Standard3FAs,achievedCoV_Standard2FAs,NPixels2Avg_Standard5FAs,NPixels2Avg_Standard4FAs,NPixels2Avg_Standard3FAs,NPixels2Avg_Standard2FAs,achievedCoV_5UniformFAs,NPixels2Avg_5UniformFAs] = NonOptFAs_CoV_NumberPixels2Avg(StdB1Range,T1Range,B1TrueRange,SNREst,CoVT1_5Angles,CoVT1_4Angles,CoVT1_3Angles,CoVT1_2Angles)

%Number of pixels to average together to achieve the desired CoV using Non optimal FAs

%   Inputs
%       StdB1Range - B1+ factor noise range
%       T1Range - Range of T1 values
%       B1TrueRange - Range of B1+ factor values
%       SNREst - SNR values
%       CoVT1_5Angles - CoV T1 vector for 5 FAs
%       CoVT1_4Angles - CoV T1 vector for 4 FAs
%       CoVT1_3Angles - CoV T1 vector for 3 FAs
%       CoVT1_2Angles - CoV T1 vector for 2 FAs

%   Outputs
%       achievedCoV_Opt5FAs: CoV achieved using 5 FAs
%       achievedCoV_Opt4FAs: CoV achieved using 4 FAs
%       achievedCoV_Opt3FAs: CoV achieved using 3 FAs
%       achievedCoV_Opt2FAs: CoV achieved using 2 FAs
%       NPixels2Avg_Opt5FAs: Number of pixels to average to achieve the desired CoV T1 using 5 FAs
%       NPixels2Avg_Opt4FAs: Number of pixels to average to achieve the desired CoV T1 using 4 FAs
%       NPixels2Avg_Opt3FAs: Number of pixels to average to achieve the desired CoV T1 using 3 FAs
%       NPixels2Avg_Opt2FAs: Number of pixels to average to achieve the desired CoV T1 using 2 FAs


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022

iTR = 1;
istdB1 = size(StdB1Range,1);%worst case stdB1
desiredCoV = 34/1025; %34ms VFA/meanT1 fibrosis stage 5-6: Table 3 repeatability chapter
nSNRs = length(SNREst);
achievedCoV_5UniformFAs = zeros(1,nSNRs);
achievedCoV_Standard5FAs = zeros(1,nSNRs);
achievedCoV_Standard4FAs = zeros(1,nSNRs);
achievedCoV_Standard3FAs = zeros(1,nSNRs);
achievedCoV_Standard2FAs = zeros(1,nSNRs);

NPixels2Avg_5UniformFAs = zeros(1,nSNRs);
NPixels2Avg_Standard5FAs = zeros(1,nSNRs);
NPixels2Avg_Standard4FAs = zeros(1,nSNRs);
NPixels2Avg_Standard3FAs = zeros(1,nSNRs);
NPixels2Avg_Standard2FAs = zeros(1,nSNRs);

%Fixed FAs
[CoVT1_3691215] = varianceSSSPGR_Fixed5FAComb([3,6,9,12,15],CoVT1_5Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_2221515] = varianceSSSPGR_Fixed5FAComb([2,2,2,15,15],CoVT1_5Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_221515] = varianceSSSPGR_Fixed4FAComb([2,2,15,15],CoVT1_4Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_2215] = varianceSSSPGR_Fixed3FAComb([2,2,15],CoVT1_3Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_215] = varianceSSSPGR_Fixed2FAComb([2,15],CoVT1_2Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);

for iNoiseVFA = 1:length(SNREst)
    

    achievedCoV_5UniformFAs(1,iNoiseVFA) = CoVT1_3691215(1,iNoiseVFA);
    NPixels2Avg_5UniformFAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_5UniformFAs(1,iNoiseVFA));
    
    achievedCoV_Standard5FAs(1,iNoiseVFA) = CoVT1_2221515(1,iNoiseVFA);
    NPixels2Avg_Standard5FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Standard5FAs(1,iNoiseVFA));
    
    achievedCoV_Standard4FAs(1,iNoiseVFA) = CoVT1_221515(1,iNoiseVFA);
    NPixels2Avg_Standard4FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Standard4FAs(1,iNoiseVFA));
    
    achievedCoV_Standard3FAs(1,iNoiseVFA) = CoVT1_2215(1,iNoiseVFA);
    NPixels2Avg_Standard3FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Standard3FAs(1,iNoiseVFA));
    
    achievedCoV_Standard2FAs(1,iNoiseVFA) = CoVT1_215(1,iNoiseVFA);
    NPixels2Avg_Standard2FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Standard2FAs(1,iNoiseVFA));
end


end

