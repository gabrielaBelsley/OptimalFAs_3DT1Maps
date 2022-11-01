function [achievedCoV_Opt5FAs,achievedCoV_Opt4FAs,achievedCoV_Opt3FAs,achievedCoV_Opt2FAs,NPixels2Avg_Opt5FAs,NPixels2Avg_Opt4FAs,NPixels2Avg_Opt3FAs,NPixels2Avg_Opt2FAs] = OptFAs_CoV_NumberPixels2Avg(OptFA_stdT1_T1B1Range_5FA,OptFA_stdT1_T1B1Range_4FA,OptFA_stdT1_T1B1Range_3FA,OptFA_stdT1_T1B1Range_2FA,SNREst,StdB1Range)

%Number of pixels to average together to achieve the desired CoV using optimal FAs

%   Inputs
%       OptFA_stdT1_T1B1Range_5FA: cell with optimal 5 FAs that minimize the
%       CoV T1 for the B1+ and T1 value that results in the max CoV T1 from the range of B1+ and T1s simulated. 
%       OptFA_stdT1_T1B1Range_4FA: cell with optimal 4 FAs that minimize the
%       CoV T1 for the B1+ and T1 value that results in the max CoV T1 from the range of B1+ and T1s simulated. 
%       OptFA_stdT1_T1B1Range_3FA: cell with optimal 3 FAs that minimize the
%       CoV T1 for the B1+ and T1 value that results in the max CoV T1 from the range of B1+ and T1s simulated. 
%       OptFA_stdT1_T1B1Range_2FA: cell with optimal 2 FAs that minimize the
%       CoV T1 for the B1+ and T1 value that results in the max CoV T1 from the range of B1+ and T1s simulated. 
%       SNREst - SNR value
%       StdB1Range - B1+ factor noise range

%   Outputs
%       achievedCoV_Opt5FAs: CoV achieved using 5 Optimal FAs
%       achievedCoV_Opt4FAs: CoV achieved using 4 Optimal FAs
%       achievedCoV_Opt3FAs: CoV achieved using 3 Optimal FAs
%       achievedCoV_Opt2FAs: CoV achieved using 2 Optimal FAs
%       NPixels2Avg_Opt5FAs: Number of pixels to average to achieve the desired CoV T1 using 5 Optimal FAs
%       NPixels2Avg_Opt4FAs: Number of pixels to average to achieve the desired CoV T1 using 4 Optimal FAs
%       NPixels2Avg_Opt3FAs: Number of pixels to average to achieve the desired CoV T1 using 3 Optimal FAs
%       NPixels2Avg_Opt2FAs: Number of pixels to average to achieve the desired CoV T1 using 2 Optimal FAs


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


iTR = 1;
indexStdB1 = size(StdB1Range,1);%worst case stdB1
desiredCoV = 34/1025; %34ms VFA/meanT1 fibrosis stage 5-6: Table 3 repeatability chapter
nSNRs = length(SNREst);
achievedCoV_Opt5FAs = zeros(1,nSNRs);
achievedCoV_Opt4FAs = zeros(1,nSNRs);
achievedCoV_Opt3FAs = zeros(1,nSNRs);
achievedCoV_Opt2FAs = zeros(1,nSNRs);
NPixels2Avg_Opt5FAs = zeros(1,nSNRs);
NPixels2Avg_Opt4FAs = zeros(1,nSNRs);
NPixels2Avg_Opt3FAs = zeros(1,nSNRs);
NPixels2Avg_Opt2FAs = zeros(1,nSNRs);


for iNoiseVFA = 1:length(SNREst)
    
    
    Opt5FAs = OptFA_stdT1_T1B1Range_5FA{iNoiseVFA,iTR,indexStdB1};
    achievedCoV_Opt5FAs(1,iNoiseVFA) = Opt5FAs(1,6);
    NPixels2Avg_Opt5FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Opt5FAs(1,iNoiseVFA));
    
    Opt4FAs = OptFA_stdT1_T1B1Range_4FA{iNoiseVFA,iTR,indexStdB1};
    achievedCoV_Opt4FAs(1,iNoiseVFA) = Opt4FAs(1,5);
    NPixels2Avg_Opt4FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Opt4FAs(1,iNoiseVFA));
    
    Opt3FAs = OptFA_stdT1_T1B1Range_3FA{iNoiseVFA,iTR,indexStdB1};
    achievedCoV_Opt3FAs(1,iNoiseVFA) = Opt3FAs(1,4);
    NPixels2Avg_Opt3FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Opt3FAs(1,iNoiseVFA));
    
    Opt2FAs = OptFA_stdT1_T1B1Range_2FA{iNoiseVFA,iTR,indexStdB1};
    achievedCoV_Opt2FAs(1,iNoiseVFA) = Opt2FAs(1,3);
    NPixels2Avg_Opt2FAs(1,iNoiseVFA)=sizeROI(desiredCoV,achievedCoV_Opt2FAs(1,iNoiseVFA));
end


end

