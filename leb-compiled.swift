type file{};

(file modelOut, file paramOut) compiledLEB (file sampleIn, file paramIn, file dependencies[]){
	app{
		compiledLEB @filename(sampleIn) @filename(paramIn) @filename(modelOut) @filename(paramOut);
	}
}


(file modelOut, file paramOut) paramCompiledLEB (file sampleIn,  float f1, float f2, float f3,float f4, float f5, float f6, float f7, float f8, float f9, float f10, file dependencies[] ){
	app{
		paramCompiledLEB @filename(sampleIn)  @filename(modelOut) @filename(paramOut) f1 f2 f3 f4 f5 f6 f7 f8 f9 f10;
	}
}

file modelOut<"model_output_grid.mat">;
file paramOut<"param_output_grid.mat">;


paramEstimation(){

	float m=0.0;				// -1:1: 0.01
	float omega=0.2;		//0.1: 0.6 0.05
	float beta=0.2;			
	float alpha=1.5;
	float rho = 0.005;
	float sigma = 0.002;
	float gamma = 0.03;
	float xi = 0.6;
	float g = 0.003;			// 0.003 0.05  0.005 
	float nyu = 0.001;
	 
	
	file sampleIn<"inputs-compiled/sample.mat">;
	String dependencyFileNames="inputs-compiled/criterion.mat inputs-compiled/datec.mat  inputs-compiled/leb76a.mat inputs-compiled/newbins.mat inputs-compiled/Parcalibsesa.mat "; 
	file dependencies[] <fixed_array_mapper; files=dependencyFileNames>;
	
	float mStep=0.01;
	float omegaStep=0.05;
	float gStep=0.05;
	
	int mRange=[0:0];
	int omegaRange=[0:0];
	int gRange=[0:0];

	foreach mIndex in mRange {	
	foreach	omegaIndex in omegaRange{
	foreach gIndex in gRange{
		file outParam<single_file_mapper; file=@strcat("param-",mIndex,"-",omegaIndex,"-",gIndex,".out")>;
		file outModel<single_file_mapper; file=@strcat("model-",mIndex,"-",omegaIndex,"-",gIndex,".out")>;

		(outModel, outParam) = paramCompiledLEB (sampleIn, m+mStep*mIndex, omega+omegaIndex*omegaStep, beta, alpha, rho, sigma, gamma, xi, g+gIndex*gStep, nyu, dependencies);

	}
	}
	}
}


paramEstimation();

