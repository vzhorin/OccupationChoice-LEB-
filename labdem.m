%% Library function to compute labor demand and number of entrepreneurs 
% Inputs: wage, support and distribution of bequests:(bsp,bpr), quadratic
% production function  parameters Par
%% labdem(wage,bsp,bpr,Par)

function [Ld,E,XXB,XXX,UVar,XPar]=labdem(wage,bsp,bpr,Par)

bpts    = length(bsp); 

m        = Par(1);
alpha   = Par(2);
beta     = Par(3);
xi         = Par(4);
rho       = Par(5);
sigma   = Par(6);

%                             Miscellaneous calculations
xi_wage = xi - wage;
alpha_1 = alpha - 1;
den       = beta*rho - sigma*sigma;
ku         = max((rho*alpha_1 + sigma*xi_wage)/den,0);
lu          = max((sigma * ku + xi_wage)/rho,0);
%
XK1     = alpha + sigma*xi_wage/rho;
XK2     = sigma*sigma/rho - beta;
XK3     = xi_wage*xi_wage/rho;
XK4     = alpha*ku - beta/2*ku*ku + xi*lu - rho/2*lu*lu + sigma*lu*ku ...
                - wage*lu - ku;
XXB     = (((ku*XK2 + XK1)^2 - XK1^2)/XK2 + XK3)/2 - wage ;
XXX     = XK4 - wage;
%
E       = 0;                  % Mass of Entrepreneurs
Ld      = 0;                  % Mass of Workers
%
for i = 1:bpts,
         b  = bsp(i);

         if b <= XXB,
            xx = b + (XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(wage+b))))/XK2;
         else
            xx  = XXX;
         end

         x1 = max(0,min(1,min(b - ku,xx)));     
         H1 = x1*(m*x1 + (1-m));
         x2 = max(0,min(1,min(b,xx))); 
         H2 = x2*(m*x2 + (1-m));
         %
         E  = E + H2*bpr(i);
         Ld = Ld + lu*H1*bpr(i);
         Ld = Ld + ((sigma*b + xi_wage)/rho*(H2 - H1) - sigma/rho* ...
            ( x2*x2*(2/3*m*x2 + (1 - m)/2) ...
            - x1*x1*(2/3*m*x1 + (1 - m)/2) ))*bpr(i);
end;

XPar=[XK1,XK2,XK3];
UVar=[XK4,ku,lu];