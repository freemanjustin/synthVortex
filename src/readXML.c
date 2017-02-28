#include "synth.h"


static int	done_it;
static int	track_count;

void get_params(e *E){

	xmlDocPtr doc;
	xmlNodePtr cur;

  //FILE  *gridrXML;


	// first pass is to see how much daa we have
	// 1. find out the number of ensembles present in the xml file (E->n_ens)
	// 2. find out the track length (E->track_length)
	done_it = FALSE;
	E->n_ens = 0;
	E->track_length = 0;
	doc = xmlParseFile(E->input_xml);

	if (doc == NULL ) {
		printf("Could not open input file %s\n", E->input_xml);
        exit(1);
	}

	cur = xmlDocGetRootElement(doc);

	if (cur == NULL) {
		xmlFreeDoc(doc);
		printf("empty input file %s\n", E->input_xml);
        exit(1);
	}

	if (xmlStrcmp(cur->name, (const xmlChar *) "cxml")) {
		xmlFreeDoc(doc);
		printf("%s is the wrong type or not a valid input xml file.\nroot node != <cxml>\n", E->input_xml);
        exit(1);
	}



	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "header"))){
				_parseInputFile_header( E, doc, cur);
		}
		else if ((!xmlStrcmp(cur->name, (const xmlChar *) "data")) && (done_it == FALSE)){
            //E->nEnsemble++;

            //fprintf(E->outFile,"#Ensemble number = %d\n", E->nEnsemble);
            //fprintf(E->outFile,"# validTime latitude longitude max_speed max_radius vortexParameters\n");
            //found_it = FALSE;
			_parseInputFile_data_first_pass (E, doc, cur);
            //fprintf(E->outFile,"\n");
						done_it = TRUE;
		}
		cur = cur->next;
	}

	xmlFreeDoc(doc);

	printf("done 1st pass\n");

	if (E->n_ens == 0){
		// temp hack to make this work with non ensemble data files
		E->n_ens = 1;
	}
	printf("number of ensembles in file %d\n", E->n_ens);
	printf("tracklength is  %d\n", E->track_length);

	// malloc arrays to store xml data
	E->ens_juldates = malloc2d_double(E->n_ens, E->track_length);
	E->ens_lat = malloc2d_float(E->n_ens, E->track_length);
	E->ens_lon = malloc2d_float(E->n_ens, E->track_length);
	E->ens_max_speed = malloc2d_float(E->n_ens, E->track_length);
	E->ens_max_radius = malloc2d_float(E->n_ens, E->track_length);
	E->ens_minimumPressure = malloc2d_float(E->n_ens, E->track_length);
    E->ens_lastClosedIsobar_Pressure = malloc2d_float(E->n_ens, E->track_length);
    E->ens_lastClosedIsobar_radius = malloc2d_float(E->n_ens, E->track_length);

	// second pass
	// now actually read the data into the E struct
	doc = xmlParseFile(E->input_xml);

	if (doc == NULL ) {
		printf("Could not open input file %s\n", E->input_xml);
        exit(1);
	}

	cur = xmlDocGetRootElement(doc);

	if (cur == NULL) {
		xmlFreeDoc(doc);
		printf("empty input file %s\n", E->input_xml);
        exit(1);
	}

	if (xmlStrcmp(cur->name, (const xmlChar *) "cxml")) {
		xmlFreeDoc(doc);
		printf("%s is the wrong type or not a valid input xml file.\nroot node != <cxml>\n", E->input_xml);
        exit(1);
	}

	cur = cur->xmlChildrenNode;
	E->nEnsemble = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "data"))){

			E->track_count = 0;

			//printf("# doing Ensemble number = %d\n", E->nEnsemble);
            //fprintf(E->outFile,"#Ensemble number = %d\n", E->nEnsemble);
            //fprintf(E->outFile,"# validTime latitude longitude max_speed max_radius vortexParameters\n");
            //found_it = FALSE;
			_parseInputFile_data (E, doc, cur);
			E->nEnsemble++;
            //fprintf(E->outFile,"\n");
		}
		cur = cur->next;
	}

	xmlFreeDoc(doc);

	/*
  printf("origin lat = %f\n", E->origin_lat);
  printf("origin lon = %f\n", E->origin_lon);

  printf("min lat crossing = %f\n", E->min_lat_xing);
  printf("max lat crossing = %f\n", E->max_lat_xing);
  printf("min lon crossing = %f\n", E->min_lon_xing);
  printf("max lon crossing = %f\n", E->max_lon_xing);
	*/

}


/**************************************************************************************
 *
 * 	private functions
 *
 **************************************************************************************/

void _parseInputFile_header (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "generatingApplication"))){// polygon
			_parseInputFile_generatingApplication(E, doc, cur);
		}
        cur = cur->next;
	}
}

void _parseInputFile_generatingApplication (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "ensemble"))){// polygon
			_parseInputFile_ensemble(E, doc, cur);
		}
        cur = cur->next;
	}
}


void _parseInputFile_ensemble (e *E, xmlDocPtr doc, xmlNodePtr cur) {
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "numMembers"))){// how many ensembles
			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			sscanf((char*)key,"%d\n", &E->n_ens);
			xmlFree(key);
		}
		cur = cur->next;
	}
}

void _parseInputFile_data_first_pass (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "disturbance"))){// polygon
			_parseInputFile_disturbance_first_pass(E, doc, cur);
		}
        cur = cur->next;
	}

}


void _parseInputFile_data (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "disturbance"))){// polygon
			_parseInputFile_disturbance(E, doc, cur);
		}
        cur = cur->next;
	}

}


void _parseInputFile_disturbance_first_pass (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "fix"))){// polygon

			_parseInputFile_fix_first_pass(E, doc, cur);
		}
    cur = cur->next;
	}

}


void _parseInputFile_disturbance (e *E, xmlDocPtr doc, xmlNodePtr cur) {

	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "fix"))){// polygon
			_parseInputFile_fix(E, doc, cur);
		}
        cur = cur->next;
	}

}


void _parseInputFile_fix_first_pass(e *E, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *key;
  float   value;
  static int  beenHere = FALSE;

  cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "validTime"))){// coordinates of polygon
			E->track_length++;
		}
		cur = cur->next;
	}
}



void _parseInputFile_fix(e *E, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *key;
  float   value;
	double	jday;
  //static int  beenHere = FALSE;

  cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "validTime"))){// coordinates of polygon

			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			//fprintf(E->outFile,"%s ", (char*)key);

			//printf("xml time string is: %s\n", (char*)key);
			//2011-02-03T22:00:00+1100
			sscanf( (char*)key, "%04d-%02d-%02dT%02d:%02d:%02lf\n", &E->year, &E->month, &E->day, &E->hour, &E->minute, &E->sec);
			//printf("%02d %02d %04d %02d:%02d:%02f\n",E->month, E->day, E->year, E->hour, E->minute, E->sec);
			julday(E->month, E->day, E->year, E->hour, E->minute, E->sec, &jday);
			E->ens_juldates[E->nEnsemble][E->track_count] = jday - E->juldate_ref;
			printf("days since = %f\n", E->ens_juldates[E->nEnsemble][E->track_count]);
			//E->track_count++;
			xmlFree(key);
		}
    else if((!xmlStrcmp(cur->name, (const xmlChar *) "latitude"))){// coordinates of polygon
			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			//fprintf(E->outFile,"%s ", (char*)key);
      sscanf((char*)key,"%f\n", &value);

			E->ens_lat[E->nEnsemble][E->track_count] = value;
			//printf("in _parse_fix: E->nEnsemble = %d\n",E->nEnsemble);
			//printf("in _parse_fix: E->track_count = %d", E->track_count);
			//printf("in _parse_fix: E->ens_lat = %f\n", E->ens_lat[E->nEnsemble][E->track_count]);

			xmlFree(key);
		}
    else if((!xmlStrcmp(cur->name, (const xmlChar *) "longitude"))){// coordinates of polygon
			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			//fprintf(E->outFile,"%s ", (char*)key);
      sscanf((char*)key,"%f\n", &value);

			E->ens_lon[E->nEnsemble][E->track_count] = value;


			xmlFree(key);
		}
        else if((!xmlStrcmp(cur->name, (const xmlChar *) "cycloneData"))){// polygon
			_parseInputFile_cycloneData(E, doc, cur);
		}
		cur = cur->next;
	}
	E->track_count++;
}

void _parseInputFile_cycloneData(e *E, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *key;
    cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "maximumWind"))){// polygon
			_parseInputFile_maximumWind(E, doc, cur);
		}
        else if((!xmlStrcmp(cur->name, (const xmlChar *) "minimumPressure"))){// polygon
			_parseInputFile_minimumPressure(E, doc, cur);
		}
		else if((!xmlStrcmp(cur->name, (const xmlChar *) "lastClosedIsobar"))){// polygon
			_parseInputFile_lastClosedIsobar(E, doc, cur);
		}
		cur = cur->next;
	}
}

void _parseInputFile_maximumWind(e *E, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *key;
  float value;
    cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if((!xmlStrcmp(cur->name, (const xmlChar *) "speed"))){// coordinates of polygon
			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			sscanf((char*)key,"%f\n", &value);

			E->ens_max_speed[E->nEnsemble][E->track_count] = value;

			//fprintf(E->outFile,"%s ", (char*)key);

			xmlFree(key);
		}
        else if((!xmlStrcmp(cur->name, (const xmlChar *) "radius"))){// coordinates of polygon
			key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			//fprintf(E->outFile,"%s ", (char*)key);
      		sscanf((char*)key,"%f\n", &value);

			E->ens_max_radius[E->nEnsemble][E->track_count] = value;

      //if(value > E->max_radius)  E->max_radius = value;
			xmlFree(key);
		}
		cur = cur->next;
	}
}

void _parseInputFile_minimumPressure(e *E, xmlDocPtr doc, xmlNodePtr cur){
    xmlChar *key;
	float value;
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if((!xmlStrcmp(cur->name, (const xmlChar *) "pressure"))){
            key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			sscanf((char*)key,"%f\n", &value);
            E->ens_minimumPressure[E->nEnsemble][E->track_count] = value;

            xmlFree(key);
        }
        cur = cur->next;
    }
}

void _parseInputFile_lastClosedIsobar(e *E, xmlDocPtr doc, xmlNodePtr cur){
    xmlChar *key;
	float	value;
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if((!xmlStrcmp(cur->name, (const xmlChar *) "pressure"))){
            key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			sscanf((char*)key,"%f\n", &value);
            //E->ens_lastClosedIsobar_Pressure[E->nEnsemble][E->track_count] = value;
			fprintf(stderr, "setting last closed isobar pressure to 1004.0\n");
			E->ens_lastClosedIsobar_Pressure[E->nEnsemble][E->track_count] = 1004.0;
            xmlFree(key);
        }
		else if((!xmlStrcmp(cur->name, (const xmlChar *) "radius"))){
            key = lr_pack(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
			sscanf((char*)key,"%f\n", &value);
            E->ens_lastClosedIsobar_radius[E->nEnsemble][E->track_count] = value;

            xmlFree(key);
        }
        cur = cur->next;
    }
}