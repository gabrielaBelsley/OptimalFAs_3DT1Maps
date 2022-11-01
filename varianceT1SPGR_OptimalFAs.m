%% Exploring T1 and B1+ Multidimensional Parameter space to find Optimal VFA FAs

% The optimal FA combination is the one that minimizes the worst (maximum)
% variance in the estimated T1 over the whole T1 range and B1+ range that
% make up our parameter space.

% The variance in the T1 estimate is affected by the SNR of the VFA data
% collected, the T1 range we wish to estimate, the range and uncertainty of
% the B1+ calculated.

% Run script after running varianceT1SPGR_BeforeScanInVivo.m
% uses the stdT12Angles, stdT13Angles, stdT14Angles, stdT15Angles output
% from varianceT1SPGR_BeforeScanInVivo.m

%This script generates Figures 1,2,3 and Table 1 in PAPER

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022


clc; close all

CoVT1_2Angles = zeros(size(stdT12Angles));
CoVT1_3Angles = zeros(size(stdT13Angles));
CoVT1_4Angles = zeros(size(stdT14Angles));
CoVT1_5Angles = zeros(size(stdT15Angles));

maxstdT12FA_T1B1Range = zeros(maxFA,maxFA);
maxstdT13FA_T1B1Range = zeros(maxFA,maxFA,maxFA);
maxstdT14FA_T1B1Range = zeros(maxFA,maxFA,maxFA,maxFA);
maxstdT15FA_T1B1Range = zeros(maxFA,maxFA,maxFA,maxFA,maxFA);

index_maxstdT12FA_T1B1Range= zeros(maxFA,maxFA);
index_maxstdT13FA_T1B1Range= zeros(maxFA,maxFA,maxFA);
index_maxstdT14FA_T1B1Range= zeros(maxFA,maxFA,maxFA,maxFA);
index_maxstdT15FA_T1B1Range= zeros(maxFA,maxFA,maxFA,maxFA,maxFA);

OptFA_stdT1_T1B1Range_2FA = cell(length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
OptFA_stdT1_T1B1Range_3FA = cell(length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
OptFA_stdT1_T1B1Range_4FA = cell(length(NoiseVFARange),length(TRRange),size(StdB1Range,1));
OptFA_stdT1_T1B1Range_5FA = cell(length(NoiseVFARange),length(TRRange),size(StdB1Range,1));

%Implements Equation 9 in Paper
% Step 1. find the combination of B1+ and T1 values that give the maximum CoV T1 
% for a total of 2, 3, 4, 5 SPGR FAs. 
% Step 2. find what FA values minimize the max CoV T1 found in step 1. 

for indexNoise = 1:length(NoiseVFARange)
    for iTR = 1:length(TRRange)
        for istdB1 = 1:size(StdB1Range,1)
            
            for ialpha1 = 1:maxFA
                for ialpha2 =(ialpha1):maxFA
                    
                    for iT1 = 1:length(T1Range)
                        CoVT1_2Angles(ialpha1,ialpha2,iT1,:,indexNoise,iTR,istdB1)=stdT12Angles(ialpha1,ialpha2,iT1,:,indexNoise,iTR,istdB1)./T1Range(iT1);
                    end
                    
                    %Step 1 : find the combination of B1+ and T1 values that give the maximum CoV T1 
                    [maxstdT12FA_T1B1Range(ialpha1,ialpha2),index_maxstdT12FA_T1B1Range(ialpha1,ialpha2)] = max(squeeze(CoVT1_2Angles(ialpha1,ialpha2,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
                    
                    
                    for ialpha3 = (ialpha2):maxFA
                        
                        %we want to minimize the worst case scenario for CoV
                        %(i.e. normalized by T1) and not stdT1
                        for iT1 = 1:length(T1Range)
                            CoVT1_3Angles(ialpha1,ialpha2,ialpha3,iT1,:,indexNoise,iTR,istdB1)=stdT13Angles(ialpha1,ialpha2,ialpha3,iT1,:,indexNoise,iTR,istdB1)./T1Range(iT1);
                        end
                        
                        [maxstdT13FA_T1B1Range(ialpha1,ialpha2,ialpha3),index_maxstdT13FA_T1B1Range(ialpha1,ialpha2,ialpha3)] = max(squeeze(CoVT1_3Angles(ialpha1,ialpha2,ialpha3,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
                        
                        for ialpha4 = (ialpha3):maxFA
                            
                            for iT1 = 1:length(T1Range)
                                CoVT1_4Angles(ialpha1,ialpha2,ialpha3,ialpha4,iT1,:,indexNoise,iTR,istdB1)=stdT14Angles(ialpha1,ialpha2,ialpha3,ialpha4,iT1,:,indexNoise,iTR,istdB1)./T1Range(iT1);
                            end
                            
                            [maxstdT14FA_T1B1Range(ialpha1,ialpha2,ialpha3,ialpha4),index_maxstdT14FA_T1B1Range(ialpha1,ialpha2,ialpha3,ialpha4)] = max(squeeze(CoVT1_4Angles(ialpha1,ialpha2,ialpha3,ialpha4,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
                            
                            
                            for ialpha5 = (ialpha4):maxFA
                                
                                for iT1 = 1:length(T1Range)
                                    CoVT1_5Angles(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5,iT1,:,indexNoise,iTR,istdB1)=stdT15Angles(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5,iT1,:,indexNoise,iTR,istdB1)./T1Range(iT1);
                                end
                                % maximum CoVT1 across the B1TrueRange and T1Range for each FA
                                %IMP: min max approach over the CoV and not
                                %stdT1: stdT1=100 for T1=1000 has the same CoV as stdT1=50ms for T1=500ms
                                %although stdT1=80ms>50ms, the 80/1000(=T1)<50/500(=T1) so
                                %we don't want to search for the maximum stdT1,
                                %but for the maximum CoV within the T1/B1 parameter range
                                [maxstdT15FA_T1B1Range(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5),index_maxstdT15FA_T1B1Range(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5)] = max(squeeze(CoVT1_5Angles(ialpha1,ialpha2,ialpha3,ialpha4,ialpha5,:,:,indexNoise,iTR,istdB1)),[],'all','omitnan','linear');
                            end
                            
                            
                        end
                    end
                end
            end
            
            %Step 2: find what FA values minimize the max CoV T1 found in step 1. 
            %2FA
            [Opt2FA] = varianceSSSPGR_MinMaxStdT1(T1Range,B1TrueRange,maxstdT12FA_T1B1Range,index_maxstdT12FA_T1B1Range);
            OptFA_stdT1_T1B1Range_2FA{indexNoise,iTR,istdB1} = Opt2FA;
            %3FA
            [Opt3FA] = varianceSSSPGR_MinMaxStdT1(T1Range,B1TrueRange,maxstdT13FA_T1B1Range,index_maxstdT13FA_T1B1Range);
            OptFA_stdT1_T1B1Range_3FA{indexNoise,iTR,istdB1} = Opt3FA;
            %4 FA
            [Opt4FA] = varianceSSSPGR_MinMaxStdT1(T1Range,B1TrueRange,maxstdT14FA_T1B1Range,index_maxstdT14FA_T1B1Range);
            OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,istdB1} = Opt4FA;
            %5FA
            [Opt5FA] = varianceSSSPGR_MinMaxStdT1(T1Range,B1TrueRange,maxstdT15FA_T1B1Range,index_maxstdT15FA_T1B1Range);
            OptFA_stdT1_T1B1Range_5FA{indexNoise,iTR,istdB1} = Opt5FA;
            
        end
    end
end

% OptFA_stdT1_T1B1Range Structure in each cell
%Column1: OptFA1
%Column2: OptFA2
%Column3: OptFA3
%Column4: OptFA4
%Column5: OptFA5
%Column6: CoVT1
%Column7: T1_OptFA
%Column8: B1_OptFA


%% Table 1: Three, Four and Five optimal FA set for stdB1 = [1 2 3], T1 = 700-1200ms, B1+=0.59-1.13
clc
indexStdB1 = 1;
% TR1 = 4.1ms
iTR = 1;
Table1Opt4FA_TR4_stdB11Perc = zeros(length(NoiseVFARange),7);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt4FA_TR4_stdB11Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_4FA{iNoiseVFA,iTR,indexStdB1};
end

Table1Opt2FA_TR4_stdB11Perc = zeros(length(NoiseVFARange),5);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt2FA_TR4_stdB11Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_2FA{iNoiseVFA,iTR,indexStdB1};
end

Table1Opt3FA_TR4_stdB11Perc = zeros(length(NoiseVFARange),6);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt3FA_TR4_stdB11Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_3FA{iNoiseVFA,iTR,indexStdB1};
end

Table1Opt5FA_TR4_stdB11Perc = zeros(length(NoiseVFARange),8);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt5FA_TR4_stdB11Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_5FA{iNoiseVFA,iTR,indexStdB1};
end
%-------------------------------

indexStdB1 = 2;
% TR1 = 4.1ms
iTR = 1;
Table1Opt3FA_TR4_stdB12Perc = zeros(length(NoiseVFARange),6);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt3FA_TR4_stdB12Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_3FA{iNoiseVFA,iTR,indexStdB1};
end

Table1Opt4FA_TR4_stdB12Perc = zeros(length(NoiseVFARange),7);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt4FA_TR4_stdB12Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_4FA{iNoiseVFA,iTR,indexStdB1};
end

Table1Opt5FA_TR4_stdB12Perc = zeros(length(NoiseVFARange),8);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt5FA_TR4_stdB12Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_5FA{iNoiseVFA,iTR,indexStdB1};
end

clc

%=========TABLE 1 PAPER============
%worst case stdB1: used in the thesis to be consistent with min-max approach
indexStdB1 = size(StdB1Range,1);%worst case stdB1
iTR = 1;
Table1Opt3FA_TR4_stdB13Perc = zeros(length(NoiseVFARange),6);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt3FA_TR4_stdB13Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_3FA{iNoiseVFA,iTR,indexStdB1};
end
Table1Opt4FA_TR4_stdB13Perc = zeros(length(NoiseVFARange),7);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt4FA_TR4_stdB13Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_4FA{iNoiseVFA,iTR,indexStdB1};
end
Table1Opt5FA_TR4_stdB13Perc = zeros(length(NoiseVFARange),8);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt5FA_TR4_stdB13Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_5FA{iNoiseVFA,iTR,indexStdB1};
end
Table1Opt2FA_TR4_stdB13Perc = zeros(length(NoiseVFARange),5);
for iNoiseVFA = 1:length(NoiseVFARange)
Table1Opt2FA_TR4_stdB13Perc(iNoiseVFA,:) = OptFA_stdT1_T1B1Range_2FA{iNoiseVFA,iTR,indexStdB1};
end
round(Table1Opt2FA_TR4_stdB13Perc(:,3).*100,1)
round(Table1Opt3FA_TR4_stdB13Perc(:,4).*100,1)
round(Table1Opt4FA_TR4_stdB13Perc(:,5).*100,1)
round(Table1Opt5FA_TR4_stdB13Perc(:,6).*100,1)

%% Deoni's optimal Flip Angle formula
TR = 4.1;
T1 = mean(T1Range);
f=0.71;
E1=exp(-TR/T1);
FA1Deoni = round(acosd((f^2*E1+(1-E1^2)*sqrt(1-f^2))/(1-E1^2*(1-f^2))))
FA2Deoni =  round((acosd((f^2*E1-(1-E1^2)*sqrt(1-f^2))/(1-E1^2*(1-f^2)))))

acosd(exp(-TR/T1))
indexStdB1 = size(StdB1Range,1);%worst case stdB1
iTR = 1;
[CoVT1_Opt2FAs_Deoni] = varianceSSSPGR_Fixed2FAComb([FA1Deoni,FA2Deoni],CoVT1_2Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_Opt4FAs_Deoni] = varianceSSSPGR_Fixed4FAComb([FA1Deoni,FA1Deoni,FA2Deoni,FA2Deoni],CoVT1_4Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);

%Performance Deoni's FAs compared to our algorithm
ImprovCoV_2FAs = CoVT1_Opt2FAs_Deoni.'-Table1Opt2FA_TR4_stdB13Perc(:,3);
ImprovCoV_4FAs = CoVT1_Opt4FAs_Deoni.'-Table1Opt4FA_TR4_stdB13Perc(:,5);

%% Research Questions: Simulations chapter/paper
clc
close all
%Q1: How does the CoV T1 change as a function of total #FAs and SNR?
%1.1) How many FAs and what are the optimal FAs for a target CoV T1 that allows to differentiate fibrosis stages?
%1.2) How much worse (CoV T1 increase) are the Standard FAs compared to Optimal FAs?
    %1.2.1) How dependent are the optimal FAs on the SPGR SNR? Can we use the
    %standard FAs 2,2,15,15 for all SNRs, i.e. how is the CoVT1 increase of using
    %2,2,15,15 for the 3 different SNRs compared to optimal FAs at each SNR minimal?

desiredCoV = 38/1025; 
indexStdB1 = size(StdB1Range,1);%worst case stdB1

%******Standard FAs***********
%2FAs: [2 15]
[CoVT1_215] = varianceSSSPGR_Fixed2FAComb([2,15],CoVT1_2Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
round(CoVT1_215.*100,1)
%3FAs: [2 2 15]
[CoVT1_2215] = varianceSSSPGR_Fixed3FAComb([2,2,15],CoVT1_3Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
round(CoVT1_2215.*100,1)
%4FAs: [2 2 15 15]
[CoVT1_221515] = varianceSSSPGR_Fixed4FAComb([2,2,15,15],CoVT1_4Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
[CoVT1_221315] = varianceSSSPGR_Fixed4FAComb([2,2,13,15],CoVT1_4Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
round(CoVT1_221515.*100,1)
round(CoVT1_221315.*100,1)
%5FAs: [2 2 2 15 15]
[CoVT1_2221515] = varianceSSSPGR_Fixed5FAComb([2,2,2,15,15],CoVT1_5Angles,iTR,indexStdB1,T1Range,B1TrueRange,SNREst);
round(CoVT1_2221515.*100,1)
%******5 uniformly spaced FAs **********
[CoVT1_3691215] = varianceSSSPGR_Fixed5FAComb([3,6,9,12,15],CoVT1_5Angles,iTR,istdB1,T1Range,B1TrueRange,SNREst);


%2 Opt FAs
CoVT14OptFAs=zeros(1,length(NoiseVFARange));
for indexNoise = 1:length(NoiseVFARange)
CoVT14OptFAs(1,indexNoise)=OptFA_stdT1_T1B1Range_4FA{indexNoise,1,indexStdB1}(1,5);
end

CoVT1_221515-CoVT14OptFAs
CoVT1_221315-CoVT14OptFAs

%**********2,3,4,5 Optimal FAs**********
NumbFA = [2,3,4,5];
indexStdB1 = size(StdB1Range,1); % worst case noise in B1+ map
CovT1_NumbFA_StandardFAs = cell(length(NoiseVFARange),length(TRRange));
CovT1_NumbFA_OptimalFAs = cell(length(NoiseVFARange),length(TRRange));
for indexNoise = 1:length(NoiseVFARange)
    for iTR = 1:length(TRRange)
        CovT1_NumbFA_OptimalFAs{indexNoise,iTR} = [OptFA_stdT1_T1B1Range_2FA{indexNoise,iTR,indexStdB1}(1,3),OptFA_stdT1_T1B1Range_3FA{indexNoise,iTR,indexStdB1}(1,4),OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,indexStdB1}(1,5),OptFA_stdT1_T1B1Range_5FA{indexNoise,iTR,indexStdB1}(1,6)];%OptFA_stdT1_T1B1Range_6FA{indexNoise,iTR,indexStdB1}(1,7)];
        CovT1_NumbFA_StandardFAs{indexNoise,iTR} =  [CoVT1_215(1,indexNoise) CoVT1_2215(1,indexNoise) CoVT1_221515(1,indexNoise) CoVT1_2221515(1,indexNoise)];

    end
end


addpath('DrosteEffect-BrewerMap-ca40391')

ColorLinePlotSet1 = brewermap(9,'Set1');
ColorLinePlotSet = [ColorLinePlotSet1(9,:); ColorLinePlotSet1(2:3,:)];

%=========FIGURE 1 PAPER===========
fig=figure();
ax=axes;
hold on
grid on
cnt = 0;

legendArray = strings(1,length(NoiseVFARange));
%Standard FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot1(indexNoise)=plot(NumbFA,CovT1_NumbFA_StandardFAs{indexNoise,1}.*100,'--s','MarkerSize',10,'LineWidth',2,'Color',ColorLinePlotSet(indexNoise,:));
    cnt = cnt +1;
    legendArray(1,cnt) = string(['Standard FAs, SNR= ',num2str(SNREst(1,indexNoise),'%0.1f')]);%,'; TR=',num2str(TRRange(1,iTR))]);
end
%Optimal FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot2 = plot(NumbFA+0.05,CovT1_NumbFA_OptimalFAs{indexNoise,1}.*100,'d','MarkerSize',10,'LineWidth',2,'Color',ColorLinePlotSet1(4,:),'MarkerFaceColor',ColorLinePlotSet1(4,:));
    
end
%5 Uniform FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot3 = plot(5.5,CoVT1_3691215(1,indexNoise).*100,'o','MarkerSize',10,'LineWidth',2,'Color','k');
end
%Target CoV T1: took out for manuscript
%hPlot4=line([1.7 6.3],[desiredCoV*100 desiredCoV*100],'Color',ColorLinePlotSet1(5,:),'LineStyle','--','LineWidth',2);
%Line to connect 5 Standard FAs to 5 Uniform FAs
for indexNoise = 1:length(NoiseVFARange)
    
    line([5 5.5],[CoVT1_2221515(1,indexNoise).*100 CoVT1_3691215(1,indexNoise).*100],'Color','k','LineStyle',':','LineWidth',2);
end
hold off
ylim([0 30])
xlim([1.7 6])
ax.XTick = [2 3 4 5 5.5];
ax.XTickLabel = '';
myLabels={'2','3','4','5','5';'', '' ,'' ,'' ,'Uniformly';'', '' ,'' ,'' ,'Spaced'};
for i = 1:length(myLabels)
    text(ax.XTick(1,i), ax.YLim(1), sprintf('%s\n%s\n%s', myLabels{:,i}), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top','FontSize',16,'FontName','Arial');    
end
ax.XLabel.String = sprintf('\n\n%s', 'Number of FAs');
set(ax,'ytick',0:5:40)
set (ax,'FontSize',16,'FontName','Arial')
%hTitle = title('CoV T1 as a function of total number of FAs used in T1 Estimate');
ylabel('T_1 CoV (%)');
%xlabel('Number of FAs');
%legendArray = strcat('SNR= ',string(num2cell(SNREst)))%,'; TR=',string(num2cell(TRRange)),'ms');
%THESIS
%legend([hPlot1,hPlot2,hPlot3,hPlot4],[legendArray,'Optimal FAs','Uniformly Spaced FAs','Target CoV: 3.7%'],'Location','NorthEast','NumColumns',2,'FontSize',14,'FontName','Arial');
legend([hPlot1,hPlot2,hPlot3],[legendArray,'Optimal FAs','Uniformly Spaced FAs'],'Location','NorthEast','NumColumns',2,'FontSize',14,'FontName','Arial');

set(fig, 'Units', 'centimeters');
figpos = get(fig, 'Position');
width = figpos(3);
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
set (gca,'FontSize',16,'FontName','Arial')
xt = -0.01;
yt = 1.05;
str = {'(a)'};
%text(xt,yt,str,'Units','normalized','fontsize',16,'FontName','Arial','FontWeight','bold')
fig.Color = [1 1 1];
figName = 'Fig1a_SimCoVT1_FANumb_SNR';


%% How many pixels to average together, i.e. average VFA pixels into 1 super pixel, to achieve the desired CoV 3.3% to differentiate different fibrosis stages 

% Figure not used in paper but may be useful together with the code

desiredCoV = 38/1025; %34ms VFA/meanT1 fibrosis stage 5-6: Table 3 repeatability chapter

clc
%using the optimal FAs for each SNR how many pixels do we need to average together?
[achievedCoV_Opt5FAs,achievedCoV_Opt4FAs,achievedCoV_Opt3FAs,achievedCoV_Opt2FAs,NPixels2Avg_Opt5FAs,NPixels2Avg_Opt4FAs,NPixels2Avg_Opt3FAs,NPixels2Avg_Opt2FAs] = OptFAs_CoV_NumberPixels2Avg(OptFA_stdT1_T1B1Range_5FA,OptFA_stdT1_T1B1Range_4FA,OptFA_stdT1_T1B1Range_3FA,OptFA_stdT1_T1B1Range_2FA,SNREst,StdB1Range);

[achievedCoV_Standard5FAs,achievedCoV_Standard4FAs,achievedCoV_Standard3FAs,achievedCoV_Standard2FAs,NPixels2Avg_Standard5FAs,NPixels2Avg_Standard4FAs,NPixels2Avg_Standard3FAs,NPixels2Avg_Standard2FAs,achievedCoV_5UniformFAs,NPixels2Avg_5UniformFAs] = NonOptFAs_CoV_NumberPixels2Avg(StdB1Range,T1Range,B1TrueRange,SNREst,CoVT1_5Angles,CoVT1_4Angles,CoVT1_3Angles,CoVT1_2Angles);

NumbFA = [2,3,4,5];
NPixels2Avg_VarySNR_OptFAs = cell(length(NoiseVFARange),1);
NPixels2Avg_VarySNR_StandardFAs = cell(length(NoiseVFARange),1);
for indexNoise = 1:length(NoiseVFARange)
        NPixels2Avg_VarySNR_OptFAs{indexNoise,1} = [NPixels2Avg_Opt2FAs(1,indexNoise) NPixels2Avg_Opt3FAs(1,indexNoise) NPixels2Avg_Opt4FAs(1,indexNoise) NPixels2Avg_Opt5FAs(1,indexNoise)];
        NPixels2Avg_VarySNR_StandardFAs{indexNoise,1} = [NPixels2Avg_Standard2FAs(1,indexNoise) NPixels2Avg_Standard3FAs(1,indexNoise) NPixels2Avg_Standard4FAs(1,indexNoise) NPixels2Avg_Standard5FAs(1,indexNoise)];
end


addpath('DrosteEffect-BrewerMap-ca40391')

ColorLinePlotSet1 = brewermap(9,'Set1');
%red orange blue grey 
ColorBlindSet = [ColorLinePlotSet1(1,:); ColorLinePlotSet1(5,:);ColorLinePlotSet1(2,:); ColorLinePlotSet1(9,:)];


%Q1.1: How many pixels do you need to average to reach the target CoV of 3.3%?

fig=figure();
ax=axes;
hold on
grid on
cnt = 0;
legendArray = strings(1,length(NoiseVFARange));
%Standard FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot1(indexNoise)=plot(NumbFA,NPixels2Avg_VarySNR_StandardFAs{indexNoise,1},'--s','MarkerSize',10,'LineWidth',2,'Color',ColorLinePlotSet(indexNoise,:));
    cnt = cnt +1;
    legendArray(1,cnt) = string(['Standard FAs, SNR= ',num2str(SNREst(1,indexNoise),'%0.1f')]);%,'; TR=',num2str(TRRange(1,iTR))]);
end
%Optimal FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot2 = plot(NumbFA+0.05,NPixels2Avg_VarySNR_OptFAs{indexNoise,1},'d','MarkerSize',10,'LineWidth',2,'Color',ColorLinePlotSet1(4,:),'MarkerFaceColor',ColorLinePlotSet1(4,:));
    
end
%5 Uniform FAs
for indexNoise = 1:length(NoiseVFARange)
    hPlot3 = plot(5.5,NPixels2Avg_5UniformFAs(1,indexNoise),'o','MarkerSize',10,'LineWidth',1.5,'Color','k');
end
%Line to connect 5 Standard FAs to 5Uniform FAs
for indexNoise = 1:length(NoiseVFARange)
    
    line([5 5.5],[NPixels2Avg_Standard5FAs(1,indexNoise) NPixels2Avg_5UniformFAs(1,indexNoise)],'Color','k','LineStyle',':','LineWidth',2);
end
hold off
ylim([0 60])
% yticks(0:10:60)
xlim([1.7 5.7])
ax.XTick = [2 3 4 5 5.5];
ax.XTickLabel = '';
myLabels={'2','3','4','5','5';'', '' ,'' ,'' ,'Uniformly';'', '' ,'' ,'' ,'Spaced'};
for i = 1:length(myLabels)
    text(ax.XTick(1,i), ax.YLim(1), sprintf('%s\n%s\n%s', myLabels{:,i}), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top','FontSize',16,'FontName','Arial');    
end
ax.XLabel.String = sprintf('\n\n%s', 'Number of FAs');
set(ax,'ytick',0:5:60)
set (ax,'FontSize',16,'FontName','TimesNewRoman')
%hTitle = title('CoV T1 as a function of total number of FAs used in T1 Estimate');
ylabel('Number of Pixels Averaged');
%xlabel('Number of FAs');
%legendArray = strcat('SNR= ',string(num2cell(SNREst)))%,'; TR=',string(num2cell(TRRange)),'ms');
legend([hPlot1,hPlot2,hPlot3],[legendArray,'Optimal FAs','Uniformly Spaced FAs'],'Location','NorthEast');
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'Fig1b_SimCoVT1_Pixels2Avg';


%Add a), b) or c) to figure
xt = -0.01;
yt = 1.05;
str = {'(b)'};
text(xt,yt,str,'Units','normalized','fontsize',16,'FontName','Arial','FontWeight','bold')
set (gca,'FontSize',16,'FontName','Arial')
fig.Color = [1 1 1];


%CoV T1 increase between optimal FAs and standard FAs is at most 1.6% for
%the lowest SNR with 4 FAs. 
DiffCoV5FAs=achievedCoV_Opt5FAs-achievedCoV_Standard5FAs;
DiffCoV4FAs=achievedCoV_Opt4FAs-achievedCoV_Standard4FAs;
DiffCoV3FAs=achievedCoV_Opt3FAs-achievedCoV_Standard3FAs;
DiffCoV2FAs=achievedCoV_Opt2FAs-achievedCoV_Standard2FAs;

% Using an uniform spacing FA set of 3,6,9,12,15 degrees resulted in a
% higher CoV T1 than using only 4 optimal FAs, regardless of SNR
achievedCoV_Opt5FAs-achievedCoV_5UniformFAs
achievedCoV_Opt4FAs-achievedCoV_5UniformFAs
achievedCoV_Opt3FAs-achievedCoV_5UniformFAs

%% FIGURE 2 What is the T1NR (or CoVT1) across the T1Range and B1TrueRange for the Optimal FA Combination?
%Remember that the Optimal FA should give the maximum bound in imprecision,
%i.e. the lowest T1NR, if this was found for the largest T1 (1200ms) and the smallest B1+ factor,
%then all the other T1s in T1Range and B1+s in B1TrueRange should have higher T1NR. 


indexNoise = 1; %SNR = 25
iTR = 1; %TR = 4.1ms
indexStdB1 = size(StdB1Range,1); % worst case B1+ nose
FA1Opt_T1B1Range = OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,indexStdB1}(1,1);
FA2Opt_T1B1Range = OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,indexStdB1}(1,2);
FA3Opt_T1B1Range = OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,indexStdB1}(1,3);
FA4Opt_T1B1Range = OptFA_stdT1_T1B1Range_4FA{indexNoise,iTR,indexStdB1}(1,4);
% FA5Opt_T1B1Range = OptFA_stdT1_T1B1Range_5FA{indexNoise,iTR,indexStdB1}(1,5);

stdT1_functionB1_FixedT1 = zeros(length(B1TrueRange),length(T1Range));
close all
%=========FIGURE 2 PAPER===========
fig=figure();
ax=axes;
hold on
grid on
clear legendArray
CoVT1_functionB1 = zeros(length(B1TrueRange),length(T1Range));
for iT1 = 1:length(T1Range)
stdT1_functionB1_FixedT1(:,iT1) = squeeze(stdT14Angles(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,iT1,:,indexNoise,iTR,indexStdB1));
T1NR_functionB1 = T1Range(iT1)./stdT1_functionB1_FixedT1(:,iT1);
CoVT1_functionB1(:,iT1) = stdT1_functionB1_FixedT1(:,iT1)./T1Range(iT1);
%Note: the code below is equivalent to the above line of code
%CoVT1_functionB1_Test(:,iT1) =CoVT1_4Angles(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,iT1,:,indexNoise,iTR,indexStdB1);
end
plot(B1TrueRange,CoVT1_functionB1(:,1).*100,'--^','MarkerSize',10,'LineWidth',2);
plot(B1TrueRange,CoVT1_functionB1(:,2).*100,'--o','MarkerSize',10,'LineWidth',2);
plot(B1TrueRange,CoVT1_functionB1(:,3).*100,'--s','MarkerSize',10,'LineWidth',2);
plot(B1TrueRange,CoVT1_functionB1(:,4).*100,'--x','MarkerSize',10,'LineWidth',2);
plot(B1TrueRange,CoVT1_functionB1(:,5).*100,'--d','MarkerSize',10,'LineWidth',2);
plot(B1TrueRange,CoVT1_functionB1(:,6).*100,'--*','MarkerSize',10,'LineWidth',2);

hold off
ylim([0 30])
%ylim([8.5 13.5])
%set(gca,'ytick',4.5:0.2:6.6)
%set(gca,'xtick',0.75:0.05:1.25)
%xlim([0.73 1.27])
%hTitle = title('CoVT1 variation with T1 & B1+ for 4 Optimal FAs');%found by min(max(stdT1 across T1 Range & B1TrueRange))')
ylabel('T_1 CoV (%)');
xlabel('True B_1^+ Factor');
legendArray = strcat('T_1= ',string(num2cell(T1Range(1:end))),'ms');
%Thesis
%legend(legendArray,'Location','SouthWest');
legend(legendArray,'Location','NorthEast','NumColumns',3,'FontSize',14,'FontName','Arial');
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
figName = 'Fig2_SimCoVT1_T1B1';
set (gca,'FontSize',16,'FontName','Arial')
fig.Color = [1 1 1];


%What is the B1+ and the T1 for the Smallest CoVT1?
[minCoV,Index] = min(CoVT1_functionB1(:))
[B1Index,T1Index] = ind2sub([length(B1TrueRange),length(T1Range)],Index);
B1TrueRange(1,B1Index)
T1Range(1,T1Index)
%What is the B1+, T1 for the Largest CoVT1?
[minCoV,Index] = max(CoVT1_functionB1(:))
[B1Index,T1Index] = ind2sub([length(B1TrueRange),length(T1Range)],Index);
B1TrueRange(1,B1Index)
T1Range(1,T1Index)

%% PAPER Fig_SimCoVT1_InvestB1orSPGRFAs.eps
close all
NumbFA = [2,3,4,5];
CovT1_NumbSPGRFA_OptimalFAs = cell(size(StdB1Range,1),length(NoiseVFARange));
for indexNoiseSPGR = 1:length(NoiseVFARange)
    for indexStdB1 = 1:size(StdB1Range,1)
        CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,indexNoiseSPGR} = [OptFA_stdT1_T1B1Range_2FA{indexNoiseSPGR,iTR,indexStdB1}(1,3),OptFA_stdT1_T1B1Range_3FA{indexNoiseSPGR,iTR,indexStdB1}(1,4),OptFA_stdT1_T1B1Range_4FA{indexNoiseSPGR,iTR,indexStdB1}(1,5),OptFA_stdT1_T1B1Range_5FA{indexNoiseSPGR,iTR,indexStdB1}(1,6)];%OptFA_stdT1_T1B1Range_6FA{indexNoise,iTR,indexStdB1}(1,7)];
        
    end
end

%ALL on same plot
fig=figure();
addpath('DrosteEffect-BrewerMap-ca40391')
colors =  get(gca,'ColorOrder');
colorsPGB=colors(4:6,:);%Purple, Green, Blue
fig.Color = [1 1 1];
hold on
cnt=0;
for indexStdB1 = size(StdB1Range,1):-1:2
    cnt=cnt+1;
    hPlotSNR12pt5(1,cnt)=plot(NumbFA,CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,1}.*100,'--s','MarkerSize',10,'LineWidth',2,'color',colorsPGB(cnt,:));
end
cnt=0;
for indexStdB1 = size(StdB1Range,1):-1:2
    cnt=cnt+1;
    hPlotSNR50(1,cnt)=plot(NumbFA,CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,3}.*100,'--*','MarkerSize',10,'LineWidth',2,'color',colorsPGB(cnt,:));
end
hold off
ylim([4 26])
xlim([1.5 5.5])
xticklabels ({'2', '3', '4', '5'});
xticks(2:1:5)
%axis square
%yticks(12:3:24)
ylabel('T_1 CoV (%)');
legend('2 B_1^+ FAs, SNR=12.5','4 B_1^+ FAs, SNR=12.5', '6 B_1^+ FAs, SNR=12.5', '2 B_1^+ FAs, SNR=50','4 B_1^+ FAs, SNR=50', '6 B_1^+ FAs, SNR=50','Location','NorthEast','NumColumns',2,'FontSize',12,'FontName','Arial');
xlabel('Number of SPGR FAs');
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
set (gca,'FontSize',14,'FontName','Arial')

%% FIGURE 3 PAPER Fig_SimCoVT1_InvestB1orSPGRFAs.eps

%Three separate plots: for SNR=12.5, SNR=25 and SNR=50
%=========FIGURE 3 PAPER===========
fig=figure();
fig.Color = [1 1 1];
colors =  get(gca,'ColorOrder');

colorsPGB=colors(4:6,:);%Purple, Green, Blue
t=tiledlayout(1,3);
nexttile
hold on
cnt=0;
for indexStdB1 = size(StdB1Range,1):-1:2
        cnt=cnt+1;
plot(NumbFA,CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,1}.*100,'--s','MarkerSize',8,'LineWidth',2,'color',colorsPGB(cnt,:),'MarkerFaceColor',colorsPGB(cnt,:));
end
hold off
%xlabel('Number of FAs');
ylim([12 24])
xlim([1 6])
xticklabels ({'2', '3', '4', '5'});
xticks(2:1:5)
%axis square
%yticks(12:3:24)
ylabel('T_1 CoV (%)');
%xlabel('Number of SPGR FAs');
%legend('2 B_1^+ FAs','4 B_1^+ FAs', '6 B_1^+ FAs','Location','NorthEast','NumColumns',1,'FontSize',12,'FontName','Arial');
%legendArray = strcat('SNR= ',string(num2cell(SNREst)))%,'; TR=',string(num2cell(TRRange)),'ms');
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
set (gca,'FontSize',14,'FontName','Arial')

nexttile
hold on
cnt=0;
for indexStdB1 = size(StdB1Range,1):-1:2
    cnt=cnt+1;
plot(NumbFA,CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,2}.*100,'--s','MarkerSize',8,'LineWidth',2,'color',colorsPGB(cnt,:),'MarkerFaceColor',colorsPGB(cnt,:));
end
hold off
xlabel('Number of SPGR FAs');
ylim([6 16])
xlim([1 6])
xticklabels ({'2', '3', '4', '5'});
xticks(2:1:5)
%yticks(12:3:24)
ylabel('T_1 CoV (%)');
xlabel('Number of SPGR FAs');
%legendArray = strcat('SNR= ',string(num2cell(SNREst)))%,'; TR=',string(num2cell(TRRange)),'ms');
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
set (gca,'FontSize',14,'FontName','Arial')

nexttile
hold on
cnt=0;
for indexStdB1 = size(StdB1Range,1):-1:2
    cnt=cnt+1;
plot(NumbFA,CovT1_NumbSPGRFA_OptimalFAs{indexStdB1,3}.*100,'--s','MarkerSize',8,'LineWidth',2,'color',colorsPGB(cnt,:),'MarkerFaceColor',colorsPGB(cnt,:));
end
hold off
%xlabel('Number of SPGR FAs');
ylim([4 14])
xlim([1 6])
xticklabels ({'2', '3', '4', '5'});
xticks(2:1:5)
%axis square
%yticks(12:3:24)
ylabel('T_1 CoV (%)');
set (gca,'FontSize',14,'FontName','Arial')
legend('2 B_1^+ FAs','4 B_1^+ FAs', '6 B_1^+ FAs','Location','NorthEast','NumColumns',1,'FontSize',12,'FontName','Arial');
%legendArray = strcat('SNR= ',string(num2cell(SNREst)))%,'; TR=',string(num2cell(TRRange)),'ms');
%lg.Layout.Tile = 'East'; % <-- place legend east of tiles
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');
set(fig, 'Units', 'centimeters');
figpos = get(fig, 'Position');
width = figpos(3);
fig.Position(3)         = 24;
fig.Position(4)         = 15;

figName = 'Fig_SimCoVT1_InvestB1orSPGRFAs';

%% Plot dS/dFA: how does it change as a function of FA? What is the steepest part of the SPGR curve? 

%Range FAs to plot
    %smallest FA with B1+ inhomogeneity: 0.59*1
    %largest FA with B1+ inhomogeneity: 1.15*15
alphaNom = deg2rad(0.1:0.1:20);
T1=700;
TR=4.1;
Tau =exp(-TR/T1);
%dSignal/dFA = dfdx in varianceT1.m
dfdx = @(M0,B1,Tau,alphaNom) (M0.*cos(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1) + (M0.*Tau.*sin(B1.*alphaNom).^2.*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2;
figure()
plot(rad2deg(alphaNom),dfdx(5000,1,Tau,alphaNom))
ylabel('dS/dFA')
xlabel('FA 1')


%Helms formula: check how dT1/dFA changes as a function of FA1, for fixed
%FA2 of 15 degrees. See if the min dT1/dFA occurs for low FA1 as the
%maximum precision is when dT1 is minimum
Alpha_1 = deg2rad(0.1:0.1:10);
Alpha_2 = deg2rad(15);
B1 =1;
M0 = 1;
% Signal Helms approximation
Sig = @(M0,TR,T1,Alpha)  M0.*Alpha.*(TR/T1)./(Alpha.^2/2 + TR/T1);

S1 = Sig(M0,TR,T1,B1.*Alpha_1);
S2 = Sig(M0,TR,T1,B1.*Alpha_2);

dT1dS1 = (2*TR/B1^2).*S2.*(Alpha_2.^2-Alpha_1.^2)./(Alpha_1.*Alpha_2.*((S2.*Alpha_2-S1.*Alpha_1).^2));
dT1dS2 = -(2*TR/B1^2).*S1.*(Alpha_2.^2-Alpha_1.^2)./(Alpha_1.*Alpha_2.*((S2.*Alpha_2-S1.*Alpha_1).^2));
%close all
figure()
plot(rad2deg(Alpha_1),dT1dS1)
hold on
plot(rad2deg(Alpha_1),dT1dS2)
ylabel('dT1/dS')
xlabel('FA 1')
legend('dT1dS1','dT1dS2')

[value,index]=min(dT1dS1);
FA_mindT1dS1=rad2deg(Alpha_1(index))


%% SPGR Curve with noise in the Signal and noise in the FA

clear;close all;clc
alpha = 1:0.1:20;
alphaAcquired = [2 3 6 12 15];

T1 = 800; %ms
TR = 4.1;
M0 = 5000;
funSPGR = @(M0,TR,alpha,T1)(M0.*(sind(alpha).*(1-exp(-TR/T1)))./(1-cosd(alpha).*exp(-TR/T1)));
funErnstAngle = @(TR,T1)(acosd(exp(-TR/T1)));
SPGRSignal = funSPGR(M0,TR,alpha,T1);
noiseSTD = 6;
SPGRSignal2deg=funSPGR(M0,TR,2,T1);
SNR=SPGRSignal2deg/noiseSTD;

addpath('DrosteEffect-BrewerMap-ca40391')
ColorLinePlot = brewermap(12,'Set1');

noiseB1=0.05;
fig=figure();
hold on
plot(alpha,SPGRSignal,'-','LineWidth',2,'Color',	ColorLinePlot(2,:))
for iAlphaAcq = 1:length(alphaAcquired)
    indexAlphaAcq=find(alpha==alphaAcquired(1,iAlphaAcq));
    x=alpha(indexAlphaAcq);
    NoiselessSignal=SPGRSignal(indexAlphaAcq);
    %noiseinSignal=SignalNoise(indexAlphaAcq)-NoiselessSignal;
    noiseFA=x*noiseB1;
    plot(x,SPGRSignal(indexAlphaAcq),'x','LineWidth',2,'Color','k','Markersize',10)
    errorbar(x,NoiselessSignal,noiseSTD,'LineWidth',2,'Color',	ColorLinePlot(1,:))
    errorbar(x,NoiselessSignal,noiseFA,'horizontal','LineWidth',2,'Color',ColorLinePlot(7,:))
end
xticks(1:3:20)
xlabel('Flip Angle (Degrees)')
ylabel('VFA SPGR Signal (a.u.)')
legend('SPGR Signal (noiseless)','Data Acquired','SPGR Signal Noise','B1+ map Noise')
set (gca,'FontSize',12,'FontName','Arial')
fig.Color = [1 1 1];
set(gca, ...
    'Box','off', ...
    'YColor','k',...
    'XColor','k',...
    'LineWidth',1,...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XGrid', 'on');



