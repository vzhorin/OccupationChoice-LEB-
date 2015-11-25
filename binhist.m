function[bin,his]=binhist(vnr, wealth, F, maxbpts)
%find bin and history

vilw=wealth;
wsp=linspace(min(vilw),max(vilw),maxbpts);

[bin_number, bin_center] = hist(vilw,wsp);
%nc=cumsum(N);

%nc1=nc-vnr;
%ind=[1:200];
%nc2=(nc1>=0).*ind;
%bin=min(nc2(nc2>0));
bin = find(bin_center>wealth(vnr), 1);
if isempty(bin) 
    bin = maxbpts; 
end
his = F(vnr,:);
