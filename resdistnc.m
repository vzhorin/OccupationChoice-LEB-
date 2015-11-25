function [bsp,bpr,ysp,ypr,Res,yprw,ypre,bprw,bpre,gridco,occo,wkco]=resdistnc(bsp0,bpr0,wage,XXB,XXX,UVar,XPar,Par)
% Function that computes the new distributions of bequests and income
% and reports some statistics.
% Inputs: Old bequest dist and support, wage, Bstar, Xstar, Unc. Var, XPar, Par
% Output: Bequest and Income distributions and supports: bsp, ysp, bpr, ypr
%         Statistics = Eu, Ku, Kc, suc, kc.  
% Revised 3/24/99
% Status: OK

maxbpts   = Par(1);
m         = Par(2);    %the parameter from the talent distr.
omega     = Par(3);    %saving rate
alpha     = Par(4);
beta      = Par(5);
xi        = Par(6);
rho       = Par(7);
sigma     = Par(8);
gamma     = Par(9);
nu        = Par(10);
%
XK4       = UVar(1); %           Unconstrained Profits
ku        = UVar(2); %           Unconstrained Capital Stock 
lu        = UVar(3); %           Unconstrained Labor Demand 
%
XK1       = XPar(1); %           Miscellaneous Calculations
XK2       = XPar(2);
XK3       = XPar(3);
%
clear UVar XPar Par;
%
Ku      = 0; %          Initialization of statistics
Kc      = 0;
Eu      = 0;
suc     = 0;
%
bpts    = maxbpts;% 100
ypts    = maxbpts;
bpts0   = length(bsp0);

%
ylo     = wage; %0.0774

b       = bsp0(bpts0);%4.3393

k       = min(b - 0,ku);
l       = max((sigma*k + xi - wage)/rho,0);

yhi     = alpha*k - beta/2*k*k + xi*l - rho/2*l*l + sigma*l*k ...
              - wage*l - k - 0;%0.5334
          
if ylo > yhi
   ylo = yhi;
   yhi = wage;
end;
blo     = omega*( bsp0(1) + ylo - nu);%0.0132
bhi     = omega*( yhi - nu + b );%1.9291


bstep   = (bhi - blo)/(bpts-1);%0.0194
ystep   = (yhi - ylo)/(ypts-1);%0.0046

bsp     = linspace(blo,bhi,bpts); 
bpr     = zeros(1,bpts);        %initialize new wealth distr. 
bprw=bpr;
bpre=bpr;

gridco=zeros(bpts,bpts);   %fractions matrix, i is columns
occo=gridco; %occupations they held
wkco=occo;


%
ysp     = linspace(ylo,yhi,ypts); 
ypr     = zeros(1,ypts);
yprw=ypr;
ypre=ypr;
% 
for i = 1:bpts0;
      b_1     = bsp0(i);
      if b_1 <= XXB,
         xx      = b_1 + (XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(wage+b_1))))/XK2;
      else,
         xx      = XXX;
      end,
       x1      = max(0,min(1,min(xx,b_1 - ku))); 
      H1      = x1*(m*x1 + 1 - m);          %fraction of unconstr. entr.
       x2      = max(0,min(1,min(b_1,xx)));
      H2      = x2*(m*x2 + 1 - m);          %fraction of all entr.?
      %
      % Computation of some statistics
      %
      Eu      = Eu+H1*bpr0(i); 
      Ku      = Ku+ku*H1*bpr0(i);
      Kc      = Kc+(b_1*(H2 - H1) - ...
                     ( x2*x2*(2/3*m*x2 + (1 - m)/2) ...
                     - x1*x1*(2/3*m*x1 + (1 - m)/2) ))*bpr0(i);
      suc     = suc + x2*x2*(2/3*m*x2 + (1 - m)/2)*bpr0(i);

      
      
% Workers
      
      % NOTE THE INCOME VARIABLE FOR WORKERS INCLUDES COST OF LIVING NU
      
      
      b       = omega*(b_1 - nu + wage);%0.0132=blo    %next period wealth of worker
      y       = wage;
      
      j      = round((b-bsp(1))/bstep + 1);   %number of nearest bin to place the guy into
      mass   = (1 - H2)*bpr0(i); %0.2306      %workers coming from point i 
      bpr(j) = bpr(j) + mass;
      bprw(j)=bprw(j) + mass;
     
      
      gridco(j,i)=gridco(j,i)+1-H2;     %transition
      occo(j,i)=1;      %occupation is worker
      wkco(j,i)=1-H2;
      
      j      = round((y-ysp(1))/ystep + 1);
      mass   = (1 - H2)*bpr0(i); 
      ypr(j) = ypr(j) + mass;
      yprw(j)=yprw(j)+mass;
      
      %
      %
      % Constrained Entrepreneurs
      
      
      % NOTE THE INCOME VARIABLE FOR WORKERS INCLUDES COST OF LIVING NU
      
      if x1<1 & ~(x1 == x2) & x2>0,
      %
         xlobar  = x1;
         xupbar  = x2;
         %
         x       = xlobar;
         k       = b_1 - x;
         l       = max((sigma*k + xi - wage)/rho,0);
         bupbar  = omega*( alpha*k - beta/2*k*k + xi*l - rho/2*l*l + ...
                        sigma*l*k - wage*l - k - x - nu + b_1);% note the cost of living is substracted also for entrepreneurs
        
         yupbar  = alpha*k - beta/2*k*k + xi*l - rho/2*l*l + ...
                        sigma*l*k - wage*l - k - x;
         x       = xupbar;
         k       = b_1 - x;
         l       = max((sigma*k + xi - wage)/rho,0);
         blobar  = omega*( alpha*k - beta/2*k*k + xi*l - rho/2*l*l + ...
                        sigma*l*k - wage*l - k - x - nu + b_1) ;
         ylobar  = alpha*k - beta/2*k*k + xi*l - rho/2*l*l + ...
                        sigma*l*k - wage*l - k - x;
                    
         %vz boundary condition to prevent index from going into negative region
         if ylobar < ysp(1)
             ylobar = ysp(1);
         end;

         if blobar < bsp(1)
             blobar = bsp(1);
         end;
         
  	     if yupbar < ysp(1)
             yupbar = ysp(1);
         end;
         
         if bupbar < bsp(1)
             bupbar = bsp(1);
         end;
         
         %
         ilobarb = round((blobar-bsp(1))/bstep + 1);
         iupbarb = round((bupbar-bsp(1))/bstep + 1);
         ilobary = round((ylobar-ysp(1))/ystep + 1);
         iupbary = round((yupbar-ysp(1))/ystep + 1);
         % 
         px      = xupbar;
         bb      = bsp(ilobarb) + bstep/2;
         nx     = b_1 + (XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(nu+bb/omega))) )/XK2;
         hpx    = px*(m*px + 1 - m);
         hnx    = nx*(m*nx + 1 - m);
         mass   = (hpx - hnx)*bpr0(i);      %failed ents coming from i
         bpr(ilobarb)  = bpr(ilobarb) + mass;
         bpre(ilobarb) = bpre(ilobarb) + mass;
         
        gridco(j,ilobarb)=gridco(j,ilobarb)+hpx-hnx;     %transition
        occo(j,ilobarb)=occo(j,ilobarb)+2;      %occupation is constr. entr.
         
         % 
         for j=(ilobarb+1):(iupbarb-1),
            px      = nx;
            bb      = bsp(j) + bstep/2;
            nx      = b_1 + (XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(nu+bb/omega))))/XK2;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            mass    = (hpx-hnx)*bpr0(i);        
            bpr(j)  = bpr(j) + mass;
            bpre(j)=bpre(j) + mass;
             
            gridco(j,i)=gridco(j,i)+hpx-hnx;     %transition
            occo(j,i)=occo(j,i)+2;      %occupation is constr. entr. 
           
         end,
         % 
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         bpr(iupbarb) = bpr(iupbarb) + mass;
         bpre(iupbarb) = bpre(iupbarb) + mass;
         
          gridco(j,iupbarb)=gridco(j,iupbarb)+hpx-hnx;     %transition
          occo(j,iupbarb)=occo(j,iupbarb)+2;      %occupation is constr. entr.
         %
         % Income distribution
	 % 
         if(ilobary < 0) 
             disp(ilobary);
         end;
    	 px     = xupbar;
         yy     = ysp(ilobary) + ystep/2;
         nx     = b_1 + ( XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(yy+b_1))) )/XK2;
         hpx    = px*(m*px + 1 - m); 
         hnx    = nx*(m*nx + 1 - m);
         mass   = (hpx - hnx)*bpr0(i);
         ypr(ilobary)  = ypr(ilobary) + mass;
         ypre(ilobary)  = ypre(ilobary) + mass;
         %         %
         for j=(ilobary+1):(iupbary-1),
            px      = nx;
            yy      = ysp(j) + ystep/2;
            nx      = b_1 + (XK1 - sqrt(XK1*XK1 - XK2*(XK3 - 2*(yy+b_1))))/XK2;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            mass    = (hpx-hnx)*bpr0(i);
            ypr(j)  = ypr(j) + mass;
            ypre(j)=ypre(j)+mass;
         end,
         %
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         ypr(iupbary) = ypr(iupbary) + mass;
         ypre(iupbary) = ypre(iupbary) + mass;
      end,
      %
      % Unconstrained Entrepreneurs
      
      
      % NOTE THE INCOME VARIABLE FOR WORKERS INCLUDES COST OF LIVING NU
      
      
      if x1>0,
         %
         xlobar  = 0;
         xupbar  = x1;
         %
         x       = xlobar;
         bupbar  = omega*( XK4 - nu - x + b_1);% note the cost of living is substracted also for entrepreneurs
         %(- x + b_1)=1.1437 note this is positive: they are unconstrained
         
         yupbar  = XK4 - x ;
         x       = xupbar;
         blobar  = omega*( XK4 - nu - x + b_1);
         ylobar  = XK4 - x ;
         %
         ilobarb = round((blobar-bsp(1))/bstep + 1);
         iupbarb = round((bupbar-bsp(1))/bstep + 1);
         ilobary = round((ylobar-ysp(1))/ystep + 1);
         iupbary = round((yupbar-ysp(1))/ystep + 1);

         %
         px      = xupbar;
         bb      = bsp(ilobarb) + bstep/2;
         nx      = XK4 - nu + b_1 - bb/omega;
         hpx     = px*(m*px + 1 - m);
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         bpr(ilobarb)  = bpr(ilobarb) + mass;
         bpre(ilobarb)=bpre(ilobarb) + mass;
         
          gridco(j,ilobarb)=gridco(j,ilobarb)+hpx-hnx;     %transition
          occo(j,ilobarb)=occo(j,ilobarb)+4;      %occupation is unconstr.
         %
         for j=(ilobarb+1):(iupbarb-1),
            px      = nx;
            bb      = bsp(j) + bstep/2;
            nx      = XK4 - nu + b_1 - bb/omega;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            mass    = (hpx-hnx)*bpr0(i);
            bpr(j)  = bpr(j) + mass;
            bpre(j)=bpre(j)+mass;
            
            
            gridco(j,i)=gridco(j,i)+hpx-hnx;     %transition
            occo(j,i)=occo(j,i)+4;      %occupation is unconstr.
         end,
         %
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         bpr(iupbarb) = bpr(iupbarb) + mass;
         bpre(iupbarb)=bpre(iupbarb)+ mass;
         
         gridco(j,iupbarb)=gridco(j,iupbarb)+hpx-hnx;     %transition
         occo(j,iupbarb)=occo(j,iupbarb)+4;      %occupation is unconstr.
         
         %
         % Income Distribution
	 %
	 px      = xupbar;
         yy      = ysp(max(ilobary,1)) + ystep/2;
         nx      = XK4 - yy;
         hpx     = px*(m*px + 1 - m);
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         ypr(max(ilobary,1))  = ypr(max(ilobary,1)) + mass;
         ypre(max(ilobary,1))  = ypre(max(ilobary,1)) + mass;
         %
         for j=(ilobary+1):(iupbary-1),
            px      = nx;
            yy      = ysp(j) + ystep/2;
            nx      = XK4 - yy;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            mass    = (hpx-hnx)*bpr0(i);
            ypr(j)  = ypr(j) + mass;
            ypre(j)  = ypre(j) + mass;
         end,
         %
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         mass    = (hpx-hnx)*bpr0(i);
         ypr(iupbary) = ypr(iupbary) + mass;
         ypre(iupbary) = ypre(iupbary) + mass;
      end,
      %
      %
end,

ke     = (-(XK1-1) + sqrt((XK1-1)*(XK1-1) - XK2*(XK3 - 2*wage)))/XK2;
kc     = max(0,ke);
Res=[Eu,Ku,Kc,suc,kc];


%Eu:fraction of unconstrained entrepreneurs
%Ku:total capital used in unconstrained firms
%Kc:total capital used in constrained firms
%suc:total setup costs
%kc:total capital used in setup costs





