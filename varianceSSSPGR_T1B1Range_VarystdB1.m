function [varT1_T1B1Range] = varianceSSSPGR_T1B1Range_VarystdB1(alphaVectorNominal,M0,nparam,TRRange,T1Range,B1TrueRange,StdB1Range,NoiseVFARange)

% Calculate the variance in T1 using the Steady-State SPGR equation for a
% range of T1s, a range of B1+s, a range of noise values in the B1+ map (stdB1) and
% a range of noise values in the VFA SPGR

%     Inputs:
%       alphaVectorNominal - vector with nominal FAs
%       M0 - M0 value
%       nparam - number of parameters to estimate: T1 and M0
%       TRRange - range of TR values
%       T1Range - range of T1 values
%       B1TrueRange - range of B1+ factor values
%       StdB1Range - range of noise values from the B1+ map
%       NoiseVFARange - range of noise values from the VFA SPGR acquisition

%     Outputs:
%       varT1_T1B1Range - T1 variance


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


varT1_T1B1Range = zeros(length(T1Range),length(B1TrueRange),length(NoiseVFARange),length(TRRange),size(StdB1Range,1));


cntTR = 0;
for iTR = TRRange
    cntTR = cntTR+1;
    
    cntT1 = 0;
    for iT1True = T1Range %loop over the T1 range
        cntT1 = cntT1+1;
        Tau = exp(-iTR/iT1True);
        
        cntB1True = 0;
        for iB1True = B1TrueRange %loop over the B1+ range
            cntB1True = cntB1True+1;
            %cntStdB1 = 0;
            for cntStdB1 = 1:size(StdB1Range,1) %loop over the stdB1+ range
                %    = cntStdB1+1;
                istdB1 = StdB1Range(cntStdB1,cntB1True);
                
                cntNoiseVFA = 0;
                for inoiseVFA =  NoiseVFARange
                    cntNoiseVFA = cntNoiseVFA+1;
                    
                    [varT1] = varianceT1(M0,iB1True,Tau,alphaVectorNominal,inoiseVFA,nparam,istdB1,iTR);
                    
                    varT1_T1B1Range(cntT1,cntB1True,cntNoiseVFA,cntTR,cntStdB1) = varT1;
                end
            end
        end
        
    end
end
end

