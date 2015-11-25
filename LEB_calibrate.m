function [wsp, wpr, ysp, ypr] = LEB_calibrate(scale, pts, in1)

%Village data
if (nargin == 1)
     pts = 25;
end;
%vz leb76a is used to calibrate wealth distribution based on 1976 SES data
load leb76a;
wealth1 = wealth;
wm = median(wealth1);
clear wealth

load(in1);
wealth = sort(wealth);
% 
resc = wm/median(wealth);
wealth = resc*wealth;   %rescale in SES units

bpts = length(income);

ymodel  = income*scale;
wmodel  = wealth*scale;         %model units village wealths

[ymodel,iinc] = sort(ymodel);
%wmodel = wmodel(iinc);
wtp = wtp(iinc);   %weights   

ylo = ymodel(1);
yhi = ymodel(bpts);
ysp = linspace(ylo,yhi,pts);
yhist = hist(ymodel,ysp);  %create histogram
ycsum = cumsum(yhist);     %compute the cumulative sum of elements

ypr(1) = sum(wtp(1:ycsum(1)));  %sum of weights for initial bin
for i=2:pts
    ypr(i) = sum(wtp(ycsum(i-1)+1:ycsum(i))); %sum of weights for bin i
end;
ty  = sum(ypr);
ypr = ypr/ty;      %income density at each bin

%[wmodel,iwlth] = sort(wmodel);
%wtp = wtp(iwlth);

wlo = wmodel(1);
whi = wmodel(length(wmodel));
wsp = linspace(wlo,whi,pts);
whist = hist(wmodel,wsp);

wpr=whist/length(wmodel);
