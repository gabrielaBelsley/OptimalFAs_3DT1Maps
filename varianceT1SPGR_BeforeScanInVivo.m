%% Calculate the Variance in T1 from Steady-State SPGR equation using the Cramér-Rao-Lower-Bound 

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022

clearvars; clc; close all

%Steady-State SPGR equation
f = @(M0,B1,Tau,alphaNom)  M0.*((sin(alphaNom.*B1).*(1-Tau))./(1-Tau.*cos(alphaNom.*B1)));
%Initialize parameters
maxFA = 15;
alpha1 = 1:maxFA;
%Calculating CoVT1 for the worst case scenario SNR:
%lowest SNR measured in a PSC Patient was 26 at 3 degrees
SNREst = [12.5 25 50];
FA_SNRMeasured = 2;%degrees
M0=5000;
T1_SNREst = 800;
TR = 4.1;
Tau_SNREst = exp(-TR./T1_SNREst);
Signal = f(M0,1,Tau_SNREst,deg2rad(FA_SNRMeasured));
NoiseVFARange = Signal./SNREst;

%Without Contrast
T1Range = 700:100:1200;
%With Contrast: T1 shortening effect
%T1Range = 300:100:600; 
B1TrueRange = 0.59:0.05:1.15;% Roberts et al. paper 
StdB1Exp=4.6;%maximum stdB1 experimental across 10 volunteers
StdB1Range = ones(1,length(B1TrueRange)).*([StdB1Exp/sqrt(4); StdB1Exp/sqrt(3); StdB1Exp/sqrt(2); StdB1Exp])./100;
nparam = 2; %unkowns: M0 and T1
TRRange = 4.1;

% initialize output variables
stdT15Angles = zeros(length(alpha1),length(alpha1),length(alpha1),length(alpha1),length(alpha1),length(T1Range),length(B1TrueRange),length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
stdT14Angles = zeros(length(alpha1),length(alpha1),length(alpha1),length(alpha1),length(T1Range),length(B1TrueRange),length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
stdT13Angles = zeros(length(alpha1),length(alpha1),length(alpha1),length(T1Range),length(B1TrueRange),length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
stdT12Angles = zeros(length(alpha1),length(alpha1),length(T1Range),length(B1TrueRange),length(NoiseVFARange),length(TRRange),size(StdB1Range,1));


%calculate the T1 variance for different number of total VFA SPGR FAs: 2, 3, 4, 5 FAs
for ialpha1 = 1:maxFA
    for ialpha2 = (ialpha1):maxFA  %FA2 > (or =) FA1 %FA1<FA2: (ialpha1+1):30
        
        if (ialpha1 == ialpha2)
            stdT12Angles(ialpha1,ialpha2)= NaN;
        else
            alphaVectorNominal_2angles = deg2rad([ialpha1,ialpha2]);
            [varT1_T1B1Range] = varianceSSSPGR_T1B1Range_VarystdB1(alphaVectorNominal_2angles,M0,nparam,TRRange,T1Range,B1TrueRange,StdB1Range,NoiseVFARange);
            stdT12Angles(ialpha1,ialpha2,:,:,:,:,:)= sqrt(varT1_T1B1Range);
            %FA,FA,T1,B1,SNR,TR,StdB1
        end
        
        for ialpha3 = (ialpha2):maxFA %FA3 > (or =) FA2 %FA2<FA3: (ialpha2+1):30
            
            if (ialpha1 == ialpha2) && (ialpha1 == ialpha3)
                stdT13Angles(ialpha1,ialpha2,ialpha3)= NaN;
            else
                alphaVectorNominal_3angles = deg2rad([ialpha1,ialpha2,ialpha3]);
                [varT1_T1B1Range] = varianceSSSPGR_T1B1Range_VarystdB1(alphaVectorNominal_3angles,M0,nparam,TRRange,T1Range,B1TrueRange,StdB1Range,NoiseVFARange);
                stdT13Angles(ialpha1,ialpha2,ialpha3,:,:,:,:,:)= sqrt(varT1_T1B1Range);
            end
            
            for ialpha4 = (ialpha3):maxFA %FA4 > (or =) FA3
                
                if (ialpha1 == ialpha2) && (ialpha1 == ialpha3) && (ialpha1 == ialpha4)
                    stdT14Angles(ialpha1,ialpha2,ialpha3,ialpha4)= NaN;
                else
                    alphaVectorNominal_4angles = deg2rad([ialpha1,ialpha2,ialpha3,ialpha4]);
                    [varT1_T1B1Range] = varianceSSSPGR_T1B1Range_VarystdB1(alphaVectorNominal_4angles,M0,nparam,TRRange,T1Range,B1TrueRange,StdB1Range,NoiseVFARange);
                    stdT14Angles(ialpha1,ialpha2,ialpha3,ialpha4,:,:,:,:,:)= sqrt(varT1_T1B1Range);
                end
                
                for ialpha5 = (ialpha4):maxFA %FA5 > (or =) FA4
                    
                    % when FA1=FA2=FA3 we do not want to estimate T1
                    if (ialpha1 == ialpha2) && (ialpha1 == ialpha3) && (ialpha1 == ialpha4) && (ialpha1 == ialpha5) 
                        stdT15Angles(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5)= NaN;
                    else
                        alphaVectorNominal_5angles = deg2rad([ialpha1,ialpha2,ialpha3,ialpha4,ialpha5]);
                        [varT1_T1B1Range] = varianceSSSPGR_T1B1Range_VarystdB1(alphaVectorNominal_5angles,M0,nparam,TRRange,T1Range,B1TrueRange,StdB1Range,NoiseVFARange);
                        stdT15Angles(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5,:,:,:,:,:)= sqrt(varT1_T1B1Range);
                        
                    end
                end               
            end
        end
    end
end





