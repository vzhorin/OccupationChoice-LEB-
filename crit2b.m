%% Criterion function
% used as a stub inside optimization/prediction procedure to freeze
% certain parameters from being varied before they get passed to actual
% model computations
%%
function[cr] = crit2b(Para, Par0, itec, scale, in1)

%% Freezing Parameters
%here we can limit optimization space by fixing parameters of input set Para

%Par = Para;

Par=[Para Par0];
%disp(Par);disp(in1);
%Par = [Par(1) 0.5791 0.0536 1.0519 0.0056 0.0001 0.0346 0.0663 0.0035 0.0009];
%Par(5) = max(0,Para(5));

%% Calling Main Modeling Module

cr = LEB_model(Par,itec,scale, in1);

%% Getting Results 
load criterion.mat
%if new criterion value is better, save as new best
if cr <= crb
    crb = cr;
    Para0 = Par;  
    disp(['Current Best is: ',num2str(crb)])
    disp('Achieved For: ')
    disp(Para)
    save -mat criterion.mat crb Para0
end
