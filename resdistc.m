function [bsp,bpr,bprw,bpre,ysp,ypr,ynsp,yprw,ypre,ynprw,ynpre,gridco,occo,wkco]=resdistc(bsp0,bpr0,rate,wage,xhat,UVar,Par)
% Function that computes the distribution of bequests and income 
% in an LEB economy with credit markets.
% NOTE IN resdistnc we compute some statistics Res. Here we don't compute them
% Input  = bsp0,bpr0,rate,wage,xhat,UVar
% Output = bsp,bpr,ysp,ypr

maxbpts = Par(1); 
m       = Par(2);
omega   = Par(3);
nu      = Par(4);
%   
fu      = UVar(1);
ku      = UVar(2);
lu      = UVar(3);
% 
bpts0   = length(bsp0);%100
pw      = 1 - xhat*(m*xhat + (1-m));%0.5189
pe      = 1 - pw;
%
bpts    = maxbpts;%100
%
% note the cost of living is substracted

blo     = omega*( bsp0(1)*rate + wage - nu );%0.0136

b       = bsp0(bpts0);%4.3393
k       = ku;
l       = lu;

% note the cost of living is substracted

bhi     = omega*( fu - wage*l + rate * ( b - k - 0 ) - nu );%2.2039 wealth of richest guy if he chooses to be an entrepreneur
%
%omega*( bsp0(bpts0)*rate + wage - nu)=2.0417 % wealth of richest guy if he chooses to be a worker

bhi     = max(bhi, omega*( bsp0(bpts0)*rate + wage - nu ));
% so the richest guy chooses to be entrepreneur

bsp     = linspace(blo,bhi,bpts); 
bstep   = (bhi - blo)/(bpts-1);
bpr     = zeros(1,bpts);% vector to locate the pdf of wealth 
bprw=bpr;
bpre=bpr;
%

%Initialize
gridco=zeros(bpts,bpts);   %fractions matrix, i is columns
occo=gridco;     %occupations they held
wkco=occo;


for i = 1:bpts0;
      b_1     = bsp0(i);
      
      %
      % Workers
      %
      b       = omega*( rate*b_1 + wage - nu );%0.0136
      
      j       = round((b-blo)/bstep + 1);
      
      
      
      meas    = pw * bpr0(i);      
      bpr(j)  = bpr(j) + meas;
      bprw(j)=bprw(j)+meas;
      %
      gridco(j,i)=gridco(j,i)+pw;     %transition
      occo(j,i)=1;      %occupation is worker
      wkco(j,i)=pw;
            %
      % Entrepreneurs
      %
      if xhat>0,
         %
         xlobar  = 0;
         xupbar  = xhat;
         %
         x       = xlobar;
         k       = ku;
         l       = lu;
         bupbar  = omega*( fu - wage*l - nu + rate*(b_1 - x - k) );%0.1758
         
         %(b_1 - x - k)=0.0055 so note this is positive, they are not constrained
         
         x       = xupbar;
         blobar  = omega*( fu - wage*l - nu + rate*(b_1 - x - k) );
         ilobar  = round((blobar-blo)/bstep + 1);
         iupbar  = round((bupbar-blo)/bstep + 1);
         %
         px      = xupbar;
         bb      = bsp(ilobar) + bstep/2;
         nx      = (fu - wage*l - nu - bb/omega)/rate + b_1 - k;
         hpx     = px*(m*px + 1 - m);
         hnx     = nx*(m*nx + 1 - m);
         meas    = (hpx-hnx)*bpr0(i);% mass of entrepreneurs
         bpr(ilobar)  = bpr(ilobar) + meas;
         bpre(ilobar)=bpre(ilobar)+meas;
         
         gridco(j,ilobar)=gridco(j,ilobar)+hpx-hnx;     %transition
         occo(j,ilobar)=occo(j,ilobar)+3;      %occupation is worker
         
         %
         for j=(ilobar+1):(iupbar-1),
            px      = nx;
            bb      = bsp(j) + bstep/2;
            nx      = (fu - wage*l - nu - bb/omega)/rate + b_1 - k;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            meas    = (hpx-hnx)*bpr0(i);
            bpr(j)  = bpr(j) + meas;
            bpre(j)=bpre(j)+meas;
            
            gridco(j,i)=gridco(j,i)+hpx-hnx;     %transition
            occo(j,i)=occo(j,i)+3;      %occupation is entr.
         end,
         %
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         meas    = (hpx-hnx)*bpr0(i);
         bpr(iupbar) = bpr(iupbar) + meas;
         bpre(iupbar)=bpre(iupbar)+meas;
         
         gridco(j,iupbar)=gridco(j,iupbar)+hpx-hnx;     %transition
         occo(j,iupbar)=occo(j,iupbar)+3;      %occupation is worker
         %
      end
      %
      %
end
%
%                                Find Distribution of Incomes
%
ypts    = maxbpts;
ynpts    = maxbpts;



%
ylo     = wage + (rate -1)*bsp0(1);%0.0784
ynlo    = wage;

k       = ku;
l       = lu;

yhi     = fu - wage*l - rate*(k + 0) + (rate-1)*bsp0(bpts0);%1.2205
ynhi= fu - rate*(k + 0)- wage*l;


%

yhi = max(yhi,wage + (rate -1)*bsp0(bpts0)); %1.2205
ynhi=max(yhi,wage);
% so the richest guy chooses to be an entrepreneur

if yhi == ylo
   yhi = yhi + 0.05;
end;
if ynhi == ynlo
   ynhi = ynhi + 0.05;
end;

ystep   = (yhi - ylo)/(ypts-1);
ynstep   = (ynhi - ynlo)/(ynpts-1);
 
%
ysp     = linspace(ylo,yhi,ypts);
ypr     = zeros(1,ypts);

ynsp     = linspace(ynlo,ynhi,ynpts);
ynprw     = zeros(1,ynpts);

ynpre=ynprw;
yprw=ynprw;
ypre=ynpre;

%
for i = 1:bpts0;
      b_1     = bsp0(i);
      %
      %
      % Workers
      %
      y       = wage + (rate-1)*b_1;
      yn       = wage;
  
      
    
      %
      j       = round((y-ylo)/ystep + 1);
      meas    = pw*bpr0(i); 
      ypr(j)  = ypr(j) + meas;
      yprw(j)  = yprw(j) + meas;
      
      jn       = round((yn-ynlo)/ynstep + 1);
     
      meas    = pw*bpr0(i); 
      ynprw(jn)  = ynprw(jn) + meas;
      
      %
      %
      %
      %
      % Entrepreneurs
      %
      if xhat>0,
         %
         xlobar  = 0;
         xupbar  = xhat;
         %
         x       = xlobar;
         yupbar  = fu - wage*l + rate*(b_1 - x - k) - b_1;
         ynupbar  = fu - wage*l - rate*(x+k);
         x       = xupbar;
         ylobar  = fu - wage*l + rate*(b_1 - x - k) - b_1;
         ynlobar  = fu - wage*l - rate*(x+k);
         ilobar  = round((ylobar-ylo)/ystep + 1);
         iupbar  = round((yupbar-ylo)/ystep + 1);
         inlobar  = round((ynlobar-ynlo)/ynstep + 1);
         inupbar  = round((ynupbar-ynlo)/ynstep + 1);
         %
         px      = xupbar;
         yy      = ysp(ilobar) + ystep/2;
         yyn      = ynsp(inlobar) + ynstep/2;
         nx      = (fu - wage*l - b_1 - yy)/rate + b_1 - k ;
         nxn     = (fu - wage*l - yyn)/rate - k ;
         
         hpx     = px*(m*px + 1 - m);
         hnx     = nx*(m*nx + 1 - m);
         hnxn     = nxn*(m*nxn + 1 - m);
         meas    = (hpx-hnx)*bpr0(i);
         measn    = (hpx-hnxn)*bpr0(i);
         
      
         ypr(ilobar)  = ypr(ilobar) + meas;
         ypre(ilobar)  = ypre(ilobar) + meas;
         ynpre(inlobar)  = ynpre(inlobar) + measn;
         %
         for j=(ilobar+1):(iupbar-1),
            px      = nx;
            yy      = ysp(j) + ystep/2;
            nx      = (fu - wage*l - b_1 - yy)/rate + b_1 - k ;
            hpx     = hnx;
            hnx     = nx*(m*nx + 1 - m);
            meas    = (hpx-hnx)*bpr0(i);
           
            ypr(j)  = ypr(j) + meas; 
            ypre(j)  = ypre(j) + meas;
            

        end,
         
         for jn=(inlobar+1):(inupbar-1),
            px      = nxn;
            yyn      = ynsp(jn) + ynstep/2;
            nxn      = (fu - wage*l - yyn)/rate - k ;
            hpxn     = hnxn;
            hnxn     = nxn*(m*nxn + 1 - m);
            measn    = (hpxn-hnxn)*bpr0(i);
            ynpre(jn)  = ynpre(jn) + measn;
            

        end
        %
         px      = nx;
         nx      = xlobar;
         hpx     = hnx;
         hnx     = nx*(m*nx + 1 - m);
         meas    = (hpx-hnx)*bpr0(i);
         ypr(iupbar) = ypr(iupbar) + meas;
         ypre(iupbar) = ypre(iupbar) + meas;
         pxn      = nxn;
         nxn      = xlobar;
         hpxn     = hnxn;
         hnxn     = nxn*(m*nxn + 1 - m);
         measn    = (hpxn-hnxn)*bpr0(i);
         ynpre(inupbar) = ynpre(inupbar) + measn;
         %
      end
      %
      %
end



