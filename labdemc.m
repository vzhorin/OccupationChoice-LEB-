%% Library function to compute labor demand and number of entrepreneurs with credit markets in place
% Inputs: wage, wealth, quadratic
% production function  parameters Par
%% labdemc(wage,wealth,Par)
function [Ld,E,rate,xhat,UVar,loanagg]=labdemc(wage,wealth,Par)

m       = Par(1);
alpha   = Par(2);
beta    = Par(3);
xi        = Par(4);
rho       = Par(5);
sigma   = Par(6);

eqeps   = 0.00001; %          Labor market clearing tolerance

rate    = 1;
lorate  = 1;
hirate  = 100;
brkflg1 = 0;
while ~(brkflg1),
         %
         %                             Miscellaneous calculations
         %
         den    = beta*rho - sigma*sigma;
         ku      = max((rho*(alpha-rate) + sigma*(xi-wage))/den,0);
         lu       = max((sigma*ku + xi-wage)/rho,0);
         %
         fu      = alpha*ku - beta/2*ku*ku + xi*lu - rho/2*lu*lu ...
                      + sigma*lu*ku ;
         xhat    = max(0,min(1,(fu - wage*lu - rate*ku - wage)/rate));
         pw      = 1 - xhat*(m*xhat + (1-m));
         pent    = 1 - pw;
         %
         suc     = xhat*xhat*(m*xhat*2/3 + (1-m)/2);
         loanagg = wealth - ku*pent - suc;
          
         %
         convergenc1 = ((loanagg >= 0) & (abs(rate-1)<eqeps)) ...
                         | (abs(loanagg) < eqeps);
         if convergenc1, brkflg1 = 1; end,
         %
         if ~(convergenc1),
            if (loanagg > 0), hirate = rate; end,
            if (loanagg <= 0), lorate = rate; end,
            rate = (lorate + hirate) /2;
         end,
         %
         if abs(lorate - hirate) < 1e-10,
            brkflg1 = 1;
            if loanagg > 0, rate = hirate; brkflg1 = 0; end,
         end,
end,
%
E       = pent ;
Ld      = pent * lu;
UVar  = [fu,ku,lu];



