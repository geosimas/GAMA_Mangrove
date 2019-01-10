/***
* Name: MangroveGrowth
* Author: Khristoffer Quinton
* Description: 
* Tags: Tag1, Tag2, TagN
***/

model MangroveGrowth

global {
	//PARAMETERS
	int yearInit <- 2018;
	int numMangrove <- 0;
	
	//GIS
	file gis_plot <- file("../includes/gisdata/plot_100m.shp");
	file gis_dem <- file("../includes/gisdata/kii_dem.asc");
	file gis_salinity <- file("../includes/gisdata/kii_salinity.asc");
	file gis_mangrove <- file("../includes/gisdata/mangrove100.shp");
	
	//initialization
	int numMangrove1 <- 20;
	int numMangrove2 <- 20;
	int numMangrove3 <- 20;

	//GIS DATA
	geometry plot <- rectangle(700,700);
	geometry shape <- envelope(plot);


	//Environment
	float pDeath_seed <- 0.1;
	float pDeath_sap <- 0.01;
	float pDeath_tree <- 0.001;
	float env_salinity <-1.0;
	float env_inundation <-1.0;
	float env_light <-1.0;

	//Environment Param
	float salt <- 35.0;
	
	init{
		create mangrove1 number: numMangrove1 with:[dbh::rnd(1.27,5)];//1.27];//
		create mangrove2 number: numMangrove2 with:[dbh::rnd(1.27,5)];//1.27];//		
		create mangrove3 number: numMangrove3 with:[dbh::rnd(1.27,5)];//1.27];//	
		}
	

	}


species mangrove control: fsm{
	
	//CONSTANT
	float G <- 162.0;
	float DBHmax <- 140.0;
	float Hmax <- 3500.0;
	float B2 <- 48.04;
	float B3 <- 0.172;
	float d <- 0.006; //should be 0.006
	float saltU <- 72.0;
	float saltD <- -0.18;
	
	float POINT_Y;
	float POINT_X;
	float elev;
	float dbh;
	float height;
	float crownRadius <-12.97;
	geometry trunkGeometry;
	geometry crownGeometry;
	geometry dispersalZone;
	int annualOffspring;
	int offspring;
	point position;
	string state <-"Seedling";
	float redFactor;
	float salt_response;

	init{
		height <- (137 + B2 * dbh - B3 * dbh^2);
		crownRadius <- 0.5*22.2*(dbh^0.654);
		crownGeometry <- sphere(crownRadius);
		trunkGeometry <- cylinder(dbh/2, height);
		position <- crownGeometry.location;
		salt_response <- 1.0/ 1.0 + exp(saltD*(saltU - salt));
		redFactor <- salt_response * env_inundation * env_light;
		annualOffspring <- int(crownGeometry.area/1000*d*redFactor);
		offspring <- offspring + annualOffspring;
	}
	
	state Seedling initial:true{
		float probdeath <- pDeath_seed;
		do grow();
		do death(probdeath);
		transition to: Sapling when:(dbh >= 0.5);}

	state Sapling{
		float probdeath <- pDeath_sap;
		do grow();
		do death(probdeath);
		transition to: Tree when:(dbh >= 2.5);
	}

	state Tree{
		float probdeath <- pDeath_tree;
		do grow();
		do death(probdeath);
		do recruit();
		transition to: OldTree when:(dbh >= 30);
	}

	state OldTree{
		float probdeath <- 0.1;
		do grow();
		do recruit();
	}

	action grow{
		//float dbhgrowth <- (G * dbh * ((1 - dbh * height)/(DBHmax * Hmax))) /(274 + 3 * B2 * dbh - (4*B3)* dbh^2);
		float dbhgrowth <- redFactor * (G * dbh * (1 - ((dbh*height) / (DBHmax*Hmax)))) /(274 + (3*B2*dbh) - (4*B3*dbh^2));
		dbh <- dbh + dbhgrowth;
		height <- (137 + B2 * dbh - B3 * dbh^2);
		crownRadius <- 0.5*22.2*(dbh^0.654);
		//write state + " " + "DBH:" + dbh + "+" +  string(dbhgrowth) + " height(m) : " + height/100 ;
	}
	
	action death(float probdeath){
		if flip(probdeath){
			do die;
		}
	}
	
	action recruit{
		write "dbh :" + string(dbh with_precision 2) + " offspring:" + string(offspring with_precision 2) + " geom:" + string(crownGeometry.area with_precision 2*0.9/10000) + " height:" + height;
		loop times: offspring{
			float range <- circle(10.0);
			geometry dispersalZone <- circle(crownRadius + range);
			point offLocation <- any_point_in(dispersalZone);
			
				if (offLocation overlaps trunkGeometry){
					point offLocation;
					loop while: (offLocation overlaps trunkGeometry){
						offLocation <- any_point_in(dispersalZone);
					}
				}
				
			float probRecruit <- rnd(1000)/1000;
			if (probRecruit>0.20){
				create species: self with:[location::offLocation, dbh::1.27, height::197];
			}
			
		}
	}
	
	
	aspect default2D{
		draw circle(crownRadius/10) color: #black;
	}
	
	aspect default3D{
	draw cylinder(dbh,height)at:{self.location.x,self.location.y, 0} color:°black;
		draw sphere(crownRadius*2) at:{self.location.x,self.location.y,height+(crownRadius*2)+0} color:°green;
	}
}
species mangrove1 parent: mangrove{
	//CONSTANT
	float G <- 162.0;
	float DBHmax <- 140.0;
	float Hmax <- 3500.0;
	float B2 <- 48.04;
	float B3 <- 0.172;
	float d <- 0.9; //should be 0.006
	float saltU <- 72.0;
	float saltD <- -0.18;	
	
	aspect default2D{
		draw circle(crownRadius/10) color: #forestgreen;
	}
	
	aspect default3D{
	draw cylinder(dbh,height)at:{self.location.x,self.location.y, 0} color:°black;
	draw cylinder(crownRadius,height) at:{self.location.x,self.location.y,height-(crownRadius*2)+0} color:°forestgreen;
		//draw cylinder(crownRadius,crownRadius) at:{self.location.x,self.location.y,height+(crownRadius*2)+0} color:°forestgreen;
	}
}

species mangrove2 parent: mangrove{
	//CONSTANT
	float G <- 243.0;
	float DBHmax <- 80.0;
	float Hmax <- 3500.0;
	float B2 <- 71.58;
	float B3 <- 0.447;
	float d <- 0.9; //should be 0.006
	float saltU <- 65.0;
	float saltD <- -0.20;
	
	aspect default2D{
		draw circle(crownRadius/10) color: #darkseagreen;
	}

	aspect default3D{
	draw cylinder(dbh,height)at:{self.location.x,self.location.y, 0} color:°black;
		draw sphere(crownRadius*2) at:{self.location.x,self.location.y,height+(crownRadius*2)+0} color:°darkseagreen;
	}
}

species mangrove3 parent: mangrove{
	//CONSTANT
	float G <- 267.0;
	float DBHmax <- 100.0;
	float Hmax <- 4000.0;
	float B2 <- 77.26;
	float B3 <- 0.396;
	float d <- 0.9; //should be 0.006
	float saltU <- 58.0;
	float saltD <- -0.25;
	
	aspect default2D{
		draw circle(crownRadius/10) color: #darkgreen;
	}

	aspect default3D{
	draw cylinder(dbh,height)at:{self.location.x,self.location.y, 0} color:°black;
	draw sphere(crownRadius) at:{self.location.x,self.location.y,height} color:°darkgreen;
	}
}
experiment main2D type: gui{
	float minimum_cycle_duration <- 0.5;

	//Param Initialise
	parameter "Initial number of Mangrove 1" var: numMangrove1 category: "Initialization" min:0 max: 50;
	parameter "Initial number of Mangrove 2" var: numMangrove2 category: "Initialization" min:0 max: 50;
	parameter "Initial number of Mangrove 3" var: numMangrove3 category: "Initialization" min:0 max: 50;	
	
	//Parameter
	parameter "Field Salinity" var: salt category: "Salinity" min: 30.0 max: 70.0;
	
	//GIS
	parameter "Shapefile for the Plot:" var: gis_plot category: "GIS" ;
	parameter "Grid File for the DEM:" var: gis_dem category: "GIS" ;
	parameter "Grid File for the Salinity:" var: gis_salinity category: "GIS" ;
	parameter "Shapefile for the Mangrove:" var: gis_mangrove category: "GIS" ;
	
	//Param Environment
	parameter "Effect of Salinity" var: env_salinity category: "Environment" min:0.5 max: 1.0;
	parameter "Effect of Inundation" var: env_inundation category: "Environment" min:0.5 max: 1.0;
	parameter "Effect of Light" var: env_light category: "Environment" min:0.5 max: 1.0;	
	
	//Param Death
	parameter "Seedling" var: pDeath_seed category: "Death Probability";
	parameter "Sapling" var: pDeath_sap category: "Death Probability";
	parameter "Tree" var: pDeath_tree category: "Death Probability";	
	
	output{
		layout #split;
				display map{
			species mangrove aspect:default2D;
			species mangrove1 aspect:default2D;
			species mangrove2 aspect:default2D;
			species mangrove3 aspect:default2D;	
		
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}	
		} // end map
		
		display speciesDistribution{
			
			chart "Species Distribution"  size: {1,1} position: {0, 0} type:pie
			{
				data "Mangrove 1" value: (mangrove1 count(each.dbh >0)) color:°green;
				data "Mangrove2" value: (mangrove2 count(each.dbh >0)) color:°forestgreen;
				data "Mangrove3" value: (mangrove3 count(each.dbh >0)) color:°darkseagreen;
				
			}
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}
					
		} //end display
		
		display graph{
			
			chart "Species Distribution"  size: {1,1} position: {0, 0} type:series
			{
				data "Mangrove 1" value: (mangrove1 count(each.dbh >0)) color:°green;
				data "Mangrove2" value: (mangrove2 count(each.dbh >0)) color:°forestgreen;
				data "Mangrove3" value: (mangrove3 count(each.dbh >0)) color:°darkseagreen;
				
			}
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}
					
		} //end display

		monitor "current year" value: yearInit + cycle;
		monitor "species1" value: (mangrove1 count(each.dbh >0));
		monitor "species2" value: (mangrove2 count(each.dbh >0));
		monitor "species3" value: (mangrove3 count(each.dbh >0));
		//inspect name:"mangrove1" value:mangrove1 attributes:["dbh","state"];
		}// end output
} // end experiment

experiment main3D type: gui{
	float minimum_cycle_duration <- 0.5;
	parameter "No of Mangrove at init:" var: numMangrove category: "General";
	
	//Param Environment
	parameter "Mangrove 1" var: numMangrove1 category: "Initialization" min:0 max: 50;
	parameter "Mangrove 2" var: numMangrove2 category: "Initialization" min:0 max: 50;
	parameter "Mangrove 3" var: numMangrove3 category: "Initialization" min:0 max: 50;	

	//Param Environment
	parameter "Effect of Salinity" var: env_salinity category: "Environment" min:0.5 max: 1.0;
	parameter "Effect of Inundation" var: env_inundation category: "Environment" min:0.5 max: 1.0;
	parameter "Effect of Light" var: env_light category: "Environment" min:0.5 max: 1.0;	
	
	
	output{
		monitor "current year" value: yearInit + cycle;
		monitor "species1" value: (mangrove1 count(each.dbh >0));
		monitor "species2" value: (mangrove2 count(each.dbh >0));
		monitor "species3" value: (mangrove3 count(each.dbh >0));
		display Plot type: opengl {
			species mangrove;
			species mangrove1 aspect: default3D;
			species mangrove2 aspect: default3D;
			species mangrove3 aspect: default3D;	
		
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}	
		}

		}
}

experiment mainGISdata type: gui{
	float minimum_cycle_duration <- 0.5;
	
	//GIS
	parameter "Shapefile for the Plot:" var: gis_plot category: "GIS" ;
	parameter "Grid File for the DEM:" var: gis_dem category: "GIS" ;
	parameter "Grid File for the Salinity:" var: gis_salinity category: "GIS" ;
	parameter "Shapefile for the Mangrove:" var: gis_mangrove category: "GIS" ;
	
	output{
		monitor "current year" value: yearInit + cycle;
		monitor "species1" value: (mangrove1 count(each.dbh >0));
		monitor "species2" value: (mangrove2 count(each.dbh >0));
		monitor "species3" value: (mangrove3 count(each.dbh >0));
		//inspect name:"mangrove1" value:mangrove1 attributes:["dbh","state"];
		display map{
			species mangrove aspect:default2D;
			species mangrove1 aspect:default2D;
			species mangrove2 aspect:default2D;
			species mangrove3 aspect:default2D;	
		
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}	
		} // end map
		
		display speciesDistribution{
			
			chart "Species Distribution"  size: {1,1} position: {0, 0} type:pie
			{
				data "Mangrove 1" value: (mangrove1 count(each.dbh >0)) color:°green;
				data "Mangrove2" value: (mangrove2 count(each.dbh >0)) color:°forestgreen;
				data "Mangrove3" value: (mangrove3 count(each.dbh >0)) color:°darkseagreen;
				
			}
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}
					
		} //end display
		
		display graph{
			
			chart "Species Distribution"  size: {1,1} position: {0, 0} type:series
			{
				data "Mangrove 1" value: (mangrove1 count(each.dbh >0)) color:°green;
				data "Mangrove2" value: (mangrove2 count(each.dbh >0)) color:°forestgreen;
				data "Mangrove3" value: (mangrove3 count(each.dbh >0)) color:°darkseagreen;
				
			}
			graphics "YEAR"{
				draw "YEAR : " +string(yearInit + cycle) color: #blue at: {50,50} font: font("Helvetica", 18 * #zoom, #bold) perspective:true;
			}
					
		} //end display

		}// end output
} // end experiment
