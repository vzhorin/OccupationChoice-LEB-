clear all;
addpath(genpath('..'));

%load business data 1-tambon_id, 2- year, 3-tambon_population, 4-retail business, 5 -cottage business
%tb_bus = load('tb_bus.txt');
%load wealth data 1-tambon_id, 2- year, 3- tambon_population, 4- flush_toilets, 5- tv_sets, 6- motorcycles, 7 - pickup_trucks
%tb_wealth = load('tb_wealth.txt');
%load business data 1-tambon_id, 2- year,  3 - baac_credit, 4- commercial_credit
%tb_fin = load('tb_credit.txt');


%1-tambon_id 2-year 3-business_proxy 4-credit_proxy 5-pca1 6-pca2 7-pca3 8-pca4 9 - population
tb_data = load('tb_data.txt');
road_data = load('dtp6_data.txt');
%tb_data = sort(tb_data, 1);
year_new = sort(unique(tb_data(:,2)))';

bus_new = [];F_new = []; wealth_new =[];
for i =1:size(year_new,2); 
	temp = []; tt=[];
	ind=find(tb_data(:,2) == year_new(i));
	bus_new = [bus_new tb_data(ind,3)];
	F_new = [F_new tb_data(ind,4)];
%find pricipal components for wealth	
	temp = [tb_data(ind,5)./tb_data(ind,9) tb_data(ind,6)./tb_data(ind,9) tb_data(ind,7)./tb_data(ind,9) tb_data(ind,8)./tb_data(ind,9)];	
 	tt = princomp(temp); 
    %tt(:,1) = (-1)^(prod(tt(:,1)<0))*tt(:,1);
    if(sign(tt(1,1)) < 0)
        tt = -tt;
    end;
 	wealth_new = [wealth_new temp*tt(:,1)];
	vilid_new = tb_data(ind,1);
end;


%get old data
load sample.mat
%clear wealth;
%clear dgr;
%clear F;
%clear bus2;
%clear vilid;

F = F_new; bus2 = bus_new; vilid = vilid_new; wealth = wealth_new(:,1); 
F2 = F2(1:5);
year = year_new;
dgr = zeros(size(vilid,1),1);

qq = quantile(road_data(:,2), [0 0.33 0.66 1]);

% bin separation stub
for i =1:size(vilid,1)
    %dgr = 3*ones(size(vilid,1),1);
    indr = find(road_data == vilid(i));
    if (indr < size(road_data,1))
        rr =  road_data(indr,2);
        if(rr >=  qq(1)) dgr(i)=1;
        end;
        if(rr >=  qq(2)) dgr(i)=2;
        end;
        if(rr >=  qq(3)) dgr(i)=3;
        end;    
    end;
end;
data2 =[vilid dgr]; 
%end of stub

save '-mat' sample_dtp6.mat year F F2 dgr bus2 vilid wealth
save '-mat' newbins_dtp6.mat data2



