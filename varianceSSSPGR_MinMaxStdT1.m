function [OptFA_stdT1_T1B1Range] = varianceSSSPGR_MinMaxStdT1(T1Range,B1TrueRange,maxCoVT1NFA_T1B1Range,index_maxCoVT1NFA_T1B1Range)

%search for minimum stdT1 and associated FA combination from the max stdT1
%previously found for each FA combination

%     Inputs:
%       T1Range - Range of T1 values
%       B1TrueRange - Range of B1+ values
%       maxCoVT1NFA_T1B1Range - maximum CoV T1 across the range of T1s and B1+ factors simulated
%       index_maxCoVT1NFA_T1B1Range - index to find out what T1, B1+ combination corresponds to the maximum CoV T1 

%     Outputs:
%       OptFA_stdT1_T1B1Range - cell with optimal FAs that minimize the
%       maximum CoV T1 found. Also included in the cell is the CoV T1
%       value, T1 value and B1+ factor corresponding to the maximum CoV T1.


% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022



FARange =  size(maxCoVT1NFA_T1B1Range);
NFA = length(FARange);
maxCoVT1NFA_T1B1Range(maxCoVT1NFA_T1B1Range == 0) = NaN;


TotalFA = FARange.*ones(1,NFA);
[MinCoVT1_T1B1Range,IndexMinCoVT1_T1B1Range] = min(maxCoVT1NFA_T1B1Range(:),[],'all','omitnan','linear');

if NFA == 6
    %get the Optimal FAs corresponding to the index found from min(maxstdT1NFA_T1B1Range(:)
    [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range,FA6Opt_T1B1Range] = ind2sub(TotalFA,IndexMinCoVT1_T1B1Range);
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxCoVT1NFA_T1B1Range(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range,FA6Opt_T1B1Range));
    T1_OptFA = T1Range(indexT1);
    B1_OptFA = B1TrueRange(indexB1);
    OptFA_stdT1_T1B1Range = [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range,FA6Opt_T1B1Range,MinCoVT1_T1B1Range,T1_OptFA,B1_OptFA];
end

if NFA == 5
    %get the Optimal FAs corresponding to the index found from min(maxstdT1NFA_T1B1Range(:)
    [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range] = ind2sub(TotalFA,IndexMinCoVT1_T1B1Range);
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxCoVT1NFA_T1B1Range(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range));
    T1_OptFA = T1Range(indexT1);
    B1_OptFA = B1TrueRange(indexB1);
    OptFA_stdT1_T1B1Range = [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,FA5Opt_T1B1Range,MinCoVT1_T1B1Range,T1_OptFA,B1_OptFA];
end


if NFA == 4
    [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range] = ind2sub(TotalFA,IndexMinCoVT1_T1B1Range);
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxCoVT1NFA_T1B1Range(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range));
    T1_OptFA = T1Range(indexT1);
    B1_OptFA = B1TrueRange(indexB1);
    OptFA_stdT1_T1B1Range = [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,FA4Opt_T1B1Range,MinCoVT1_T1B1Range,T1_OptFA,B1_OptFA];
end

if NFA == 3
    [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range] = ind2sub(TotalFA,IndexMinCoVT1_T1B1Range);
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxCoVT1NFA_T1B1Range(FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range));
    T1_OptFA = T1Range(indexT1);
    B1_OptFA = B1TrueRange(indexB1);
    OptFA_stdT1_T1B1Range = [FA1Opt_T1B1Range,FA2Opt_T1B1Range,FA3Opt_T1B1Range,MinCoVT1_T1B1Range,T1_OptFA,B1_OptFA];
end

if NFA == 2
    [FA1Opt_T1B1Range,FA2Opt_T1B1Range] = ind2sub(TotalFA,IndexMinCoVT1_T1B1Range);
    [indexT1,indexB1] = ind2sub([length(T1Range) length(B1TrueRange)],index_maxCoVT1NFA_T1B1Range(FA1Opt_T1B1Range,FA2Opt_T1B1Range));
    T1_OptFA = T1Range(indexT1);
    B1_OptFA = B1TrueRange(indexB1);
    OptFA_stdT1_T1B1Range = [FA1Opt_T1B1Range,FA2Opt_T1B1Range,MinCoVT1_T1B1Range,T1_OptFA,B1_OptFA];
end
% OptFA_stdT1_T1B1Range Structure
%Column1: OptFA1
%Column2: OptFA2
%Column3: OptFA3
%Column4: OptFA4
%Column5: OptFA5
%Column6: CoVT1
%Column7: stdT1/T1_OptFA eliminated: only makes sense when we doing max min search over stdT1 instead of CoV
%Column8: T1_OptFA
%Column9: B1_OptFA

end

