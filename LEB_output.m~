function LEB_output(in1, out1, out2)

load(in1);
load ent.mat;
load criterion.mat;
load results.mat;

%do your own formatting here
%out = fopen(out1, 'w+');
%fprintf(out, '%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s%1s%12s\n', 'village_id', ',', 'model_ent_1988',...
%    ',','data_ent_1988',  ',', 'model_ent_1990', ',','data_ent_1990', ',', 'model_ent_1992',...
%    ',','data_ent_1992',',', 'model_ent_1994', ',','data_ent_1994',',', 'model_ent_1996', ',','data_ent_1996');
%for i = 1:size(vilid);
%	fprintf(out, '%8i%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f%1s%8.6f\n', vilid(i), ',', en0(i,2), ',',bus2(i,2),',', en0(i,3), ',',bus2(i,3), ',', en0(i,4), ',',bus2(i,4),',', en0(i,5), ',',bus2(i,5), ',', en0(i,6), ',',bus2(i,6));
%end;
%fclose(out);

out = fopen(out2, 'w+');
fprintf(out, '%12s%8.6f\n','criterion=',crb);
fprintf(out, '%30s\n','Parameters of the model');
for i = 1:size(Para0');
	fprintf(out, '%12.3e\n',Para0(i));
end;
fclose(out);

