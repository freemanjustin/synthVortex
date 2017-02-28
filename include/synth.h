// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <errno.h>


// libxml2 headers
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>


// netCDF header
#include "netcdf.h"

#include "jutil.h"

#include "vortex.h"

// macros
#define	TRUE 1
#define FALSE 0

typedef struct{

    char	 *input_xml;
    char     *force_out;
    double   *lon;
    double   *lat;


    // xmlIO params
    int nEnsemble;
	int track_count;
    // time parsing variables
	int month;
	int day;
	int year;
	int hour;
	int minute;
	double sec;
	double juldate;
	double juldate_ref;

    int		n_ens;
	int		track_length;

    double	**ens_juldates;
	float	**ens_lat;
	float	**ens_lon;
	float	**ens_max_speed;
	float	**ens_max_radius;
	float	**ens_minimumPressure;
    float	**ens_lastClosedIsobar_Pressure;
    float	**ens_lastClosedIsobar_radius;

}e;



// xmlIO.c
void get_params(e*);

// prototypes
void _parseInputFile_header (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_generatingApplication (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_ensemble (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_params (e*, xmlDocPtr, xmlNodePtr);
void _parseInputFile_data (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_disturbance (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_fix(e *E, xmlDocPtr doc, xmlNodePtr cur);

void _parseInputFile_vortexParameters(e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_cycloneData(e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_maximumWind(e *E, xmlDocPtr doc, xmlNodePtr cur);

void _parseInputFile_data_first_pass (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_disturbance_first_pass (e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_fix_first_pass(e *E, xmlDocPtr doc, xmlNodePtr cur);

void _parseInputFile_minimumPressure(e *E, xmlDocPtr doc, xmlNodePtr cur);
void _parseInputFile_lastClosedIsobar(e *E, xmlDocPtr doc, xmlNodePtr cur);

