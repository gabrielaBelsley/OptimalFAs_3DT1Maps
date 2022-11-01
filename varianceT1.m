function [varT1] = varianceT1(M0,B1,Tau,alphaNom,noiseVFA,nparam,stdB1,TR)

%Implements equation 7 in paper to calculate the T1 variance 

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2022

%     Inputs:
%       M0 - M0 value
%       B1 - B1+ factor 
%       Tau - exp(-TR/T1)
%       alphaNom - vector with nominal FAs
%       nparam - number of parameters to estimate: T1 and M0
%       stdB1 - noise value from the B1+ map
%       TR - TR value

%     Outputs:
%       varT1 - T1 variance

% Calculate the derivatives with the symbolic toolbox
% Advantage of symbolic toolbox: gives you the derivative as if you derived on paper
% non-symbolic toolbox : gives the derivative evaluated at certain points
%f = @(M0,B1,Tau,alphaNom)  M0.*((sin(alphaNom.*B1).*(1-Tau))./(1-Tau.*cos(alphaNom.*B1)));
%Tau = exp(-TR/T1)

% syms M0 B1 Tau alphaNom real
% f = M0.*((sin(alphaNom.*B1).*(1-Tau))./(1-Tau.*cos(alphaNom.*B1)));
% dfdM0 = diff(f,M0)
% dfdM0 = (sin(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1);
% dfdB1 = diff(f,B1)
% dfdB1 =(M0.*alphaNom.*cos(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1) + (M0.*Tau.*alphaNom.*sin(B1.*alphaNom).^2.*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2
% dfdTau = diff(f,Tau)
% dfdTau = (M0.*sin(B1.*alphaNom))./(Tau.*cos(B1.*alphaNom) - 1) - (M0.*cos(B1.*alphaNom).*sin(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2;

% syms M0 x Tau real
% f = M0.*((sin(x).*(1-Tau))./(1-Tau.*cos(x)));
% dfdx = diff(f,x)
% dfdx = (M0.*cos(x).*(Tau - 1))./(Tau.*cos(x) - 1) + (M0.*Tau.*sin(x).^2.*(Tau - 1))./(Tau.*cos(x) - 1).^2;
%x = B1.*alphaNom

dfdM0 = @(M0,B1,Tau,alphaNom) (sin(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1);
dfdTau = @(M0,B1,Tau,alphaNom) (M0.*sin(B1.*alphaNom))./(Tau.*cos(B1.*alphaNom) - 1) - (M0.*cos(B1.*alphaNom).*sin(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2;
dfdB1 = @(M0,B1,Tau,alphaNom) (M0.*alphaNom.*cos(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1) + (M0.*Tau.*alphaNom.*sin(B1.*alphaNom).^2.*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2;
dfdx = @(M0,B1,Tau,alphaNom) (M0.*cos(B1.*alphaNom).*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1) + (M0.*Tau.*sin(B1.*alphaNom).^2.*(Tau - 1))./(Tau.*cos(B1.*alphaNom) - 1).^2;

% compute the weights: uncertainties in x and y
% x = alphaNom*B1 which has an uncertainty of stdB1
% y uncertainty from noise VFA
% sigmaAlphaTrue = alphaNom.*stdB1
w = @(alphaNom) 1./(noiseVFA^2 + (dfdx(M0,B1,Tau,alphaNom).*alphaNom.*stdB1).^2);

% Compute Matrix C:
% [dfdM0dfdM0      dfdM0dfdTau]
% [dfdTaudfdM0     dfdTaudfdTau]

C11 = sum(w(alphaNom).*(dfdM0(M0,B1,Tau,alphaNom)).^2,2); %sum along the columns
C12 = sum(w(alphaNom).*dfdM0(M0,B1,Tau,alphaNom).*dfdTau(M0,B1,Tau,alphaNom),2);
C21 = C12;
C22 = sum(w(alphaNom).*(dfdTau(M0,B1,Tau,alphaNom)).^2,2);


C = [C11 C12;C21 C22];

%Cinv = C\eye(2);
%Cinv_Analytical = C11/(C11*C22-C12*C21);
%in practice C\eye(2) and Cinv_Analytical give the same answer to within e-20

%Rationale for if (abs(C11*C22-C12*C21)/C11*C22)>1e-10 : check notebook 27082020 yellow box
%if the determinant C = the denominator of Cinv_Analytical can be determined
%and is larger than 1e-10 then the denominator can accuratly be calculated
%and consequently Cinv. Normalized the determinant by C11C22 as this shows how
%big the difference is relative to C11C22 because the larger C11C22 the
%more digits it will need to be saved so the less digits are available to calculate the difference
%e.g.
%if C11C22=1e5
%C12*C21 = 1e5 + eps*10^2:
%the difference C11*C22-C12*C21 can't be calculated
%1e5*(1-(1+2.2204e-17)) 2.2204e-17<eps (yellow box in notebook)
%1e5*(1-(1+2.2204e-17))/1e5 > 1e-10 (FALSE =>stdT1 = NaN)


if (abs(C11*C22-C12*C21)/C11*C22)>1e-10 %(1e6*larger than eps= 2.220446049250313e-16)
    
    %Other if statement used in the past
    
    %---1.if norm(Cinv*C - eye(2))< 1e-6
    %Worry: when the condition number of C is very large, the matrix is
    %ill-conditioned and Cinv calculation might be innacurate which would
    %lead to innacurate stdT1
    %To check if Cinv is accurate we do: Cinv*C = Identity
    
    %---2. if cond(C) > 1e14
    %     varT1 = NaN;
    % else
    
    Cinv_Analytical = 1/(C11*C22-C12*C21) .* [C22 -C12;-C21 C11];
    
    S = length(alphaNom);%sum(w(alphaNom).*(noiseVFA.^2),2);
    
    %variance in Tau due to SNR VFA
    varTau = (Cinv_Analytical(2,2)); % CRLB
    stdTau = sqrt(varTau);
    stdT1 = (TR/((log(Tau)).^2.*Tau)).*stdTau;
    varT1 = stdT1.^2;
    
else
    varT1 = NaN;
    
end

end












