function [E1,E2,Ld1,Ld2,wage,rate,XXB1,XXX1,xhat2,UVar1,UVar2,XPar1]=lmeq2sc(bsp1,bpr1,wealth2,Par1,Par2,FC)
% Labor Market EQuilibrium for 2 Sector Credit Economy.
%
% Functions used:
% Labdem.m, labdemc.m : computes labor demand for the sector 
%Input: support and distribution of bequests:(bsp,bpr)  
%       Parameters (Par) = [m alpha beta xi rho sigma gamma nu] for each sector. 
%Output:Entrepreneurs (E), Workers (Ld), wage, Bstar (XXB), Xstar (XXX),
%       Unconstrained Variables=Profits,Capital,Labor:  UVar=[XK4,ku,lu]
%       X Parameters:XPar=[XK1,XK2,XK3] in each sector.
%
%Revised 3/24/99
%Status: OK

eqeps   = 0.00001; %          Labor market clearing tolerance
% 
rwage   = Par1(7)+Par1(8);
% 
lowage  = rwage;
hiwage  = rwage + 5; 
wage    = rwage;
%
brkflag = 0;
while ~(brkflag),
    %
    [Ld1,E1,XXB1,XXX1,UVar1,XPar1] = labdem(wage,bsp1,bpr1,Par1);
    [Ld2,E2,rate,xhat2,UVar2] = labdemc(wage,wealth2,Par2);
    %
    Ld = (1-FC)*Ld1+FC*Ld2;
    E  = (1-FC)*E1+FC*E2;
    convergence = ( abs(Ld + E - 1) < eqeps);
    if convergence, 
       brkflag = 1; 
        
    end,
    %
    if ~(convergence),
         if (Ld + E < 1), 
            hiwage = wage; 
            HLd = Ld; 
            HE = E; 
         end,
         if (Ld + E > 1), 
            lowage = wage; 
         end,
         wage = (lowage + hiwage)/2;
    end,
    %                      Now, check if equilibrium is underemployment.
    if (abs(hiwage - lowage) < 1e-10),
         if (Ld + E < 1), 
            brkflag = 1; 
         else, 
            wage = hiwage; 
            Ld = HLd; 
            E = HE; 
            brkflag = 1;
         end,
    end,      
    %
    % 
end, 

%
Ld1 = Ld1*(1-FC);
Ld2 = Ld2*FC;
%
E1 = E1*(1-FC);
E2 = E2*FC;
%

