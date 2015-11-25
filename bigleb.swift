type file{};

(file gisData) queryGIS (file queryCommand){
	app {
		queryGIS @filename(queryCommand) @filename(gisData);
	}
}

(file dataFile) xmlXtract (file xmlData){
	app {
		xmlXtract @filename(xmlData) @filename(dataFile);
	}
}

(file lebResult) biglebExecute(file mainScript, file supportFiles[], file dataFiles[], file gisData, int rangeStart, int rangeSize){
	app{
		biglebInvoke @filename(mainScript) rangeStart rangeSize @filename(lebResult);
	}
}

(file solutions[]) batchLebExecute (file gisTxtData){

	file lebMain<"scripts/rescalSES18c.m">;
	file lebScripts[] <filesys_mapper; location="scripts-lebbig", suffix=".m">;
	file lebData[] <filesys_mapper; location="inputs", suffix=".mat">;

    int rangeSize=2;
    int rangeStart=[0:1];

	solutions[0]=biglebExecute(lebMain, lebScripts, lebData, gisTxtData, 1, rangeSize);


    //foreach i in rangeStart {
    //    int position=i*rangeSize;
    //     solutions[i]=biglebExecute(lebMain, lebScripts, lebData, gisTxtData, position, rangeSize);
    //}
}


file queryScript<"scripts/gis.query.txt">;

file gisXmlData<"gis.xml">;
file gisTxtData<"gis.txt">;

file villageLEBSolutions[] <fixed_array_mapper; files="villageLEB.0.txt">;
//file villageLEBSolutions[] <fixed_array_mapper; files="villageLEB.0.txt, villageLEB.1.txt">;
string outputFileName="villageLEB";

//gisXmlData=queryGIS(queryScript);
//gisTxtData=xmlXtract(gisXmlData);
villageLEBSolutions=batchLebExecute(gisTxtData);



