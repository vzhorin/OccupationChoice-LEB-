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

(file lebResult) lebExecute(file mainScript, file supportFiles[], file dataFiles[], file gisData){
	app{
		lebInvoke @filename(mainScript) stdout=@filename(lebResult);
	}
}

file queryScript<"scripts/gis.query.txt">;
file lebMain<"scripts/rescalSES18c.m">;
file lebScripts[] <filesys_mapper; location="scripts", suffix=".m">;
file lebData[] <filesys_mapper; location="inputs", suffix=".mat">;

file gisXmlData<"gis.xml">;
file gisTxtData<"gis.txt">;
file lebOut<"lebOut.txt">;

gisXmlData=queryGIS(queryScript);
gisTxtData=xmlXtract(gisXmlData);
lebOut = lebExecute(lebMain, lebScripts, lebData, gisTxtData);
