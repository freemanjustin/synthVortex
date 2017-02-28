#include "synth.h"

int main(int argc, char **argv)
{
    e *E;
    int i, j;
    int ncid, varid, dimid;
    int lat_dimid, lon_dimid, time_dimid, dimIds[3];
    int lat_varid, lon_varid, time_varid;
    int taux_varid, tauy_varid, pressure_varid;

    int nLat, nLon;
    double *lon, *lat;

    double lon_start, lon_end;
    double lat_start, lat_end;

    int fix;

    E = malloc(sizeof(e));
    if (E == NULL)
    {
        fprintf(stderr, "malloc failed for E struct");
        exit(1);
    }

    // parse command line arguments
    if (argc < 3)
    {
        printf("usage:\n");
        printf("\tsynth [storm parameter txt file] [roms forcing netcdf file]\n");
        exit(1);
    }
    else
    {
        get_command_line_arg_as_string(&E->input_xml, argv[1]);
        get_command_line_arg_as_string(&E->force_out, argv[2]);
    }

    // set up reference julian date
    // ref date is 1970-01-01 00:00:0.0

    julday(1, 1, 1970, 0, 0, 0, &E->juldate_ref);
    // read the input xml file
    get_params(E);

    // dump the data to the term to check
    /*
    for(i=0;i<E->track_length;i++){

        printf("%f %f %f %f %f %f %f %f\n", E->ens_juldates[0][i], E->ens_lat[0][i], E->ens_lon[0][i],
                            E->ens_max_speed[0][i], E->ens_max_radius[0][i],
                            E->ens_minimumPressure[0][i], E->ens_lastClosedIsobar_Pressure[0][i],
                            E->ens_lastClosedIsobar_radius[0][i]);

    }
    */

    // set up target grid dimensions
    // access-r bounds:
    double slon = 95.0, elon = 170.0;
    double slat = -55.0, elat = 5.0;
    double dlat = 0.05, dlon = 0.05;
    int nx, ny;
    double *glon, *glat;  //grid lon/lat
    double *ang, *radius; // radians, km
    double *p, *v;        //pressure and Holland velocity
    double *vx, *vy;
    double *tau_x, *tau_y;
    double *rv;         //Rankine velocity
    double alpha = 0.4; //Rankine parameter
    double lat1, lon1;
    double Rm = 31.3;  //radius of vmax (km)
    double Vm = 56.0;  //vmax (m/s)
    double rho = 1.15; //density of air (kg/m^3)
    double pc, pn;     //central and outer pressures (mb)
    double cor;        //coriolis force;
    FILE *fp;

    slon = 95.0; elon = 170.0;
    slat = -55.0; elat = 5.0;
    dlat = 0.05; dlon = 0.05;

    nx = (int)((elon - slon) / dlon);
    ny = (int)((elat - slat) / dlat);

    dlon = (elon - slon) / nx; //recalculate now
    dlat = (elat - slat) / ny; //recalculate now

    // malloc arrays
    glon = malloc(nx * sizeof(double));
    glat = malloc(ny * sizeof(double));
    ang = malloc(nx * ny * sizeof(double));
    radius = malloc(nx * ny * sizeof(double));
    p = malloc(nx * ny * sizeof(double));
    v = malloc(nx * ny * sizeof(double));
    vx = malloc(nx * ny * sizeof(double));
    vy = malloc(nx * ny * sizeof(double));
    tau_x = malloc(nx * ny * sizeof(double));
    tau_y = malloc(nx * ny * sizeof(double));
    //rv = malloc(nx * ny * sizeof(double));
    //MCK(rv);

    //printf("ang mem =  %fmb\n", nx * ny * sizeof(double) / 1024.0 / 1024.0);

    for (i = 0; i < nx; i++){
        glon[i] = slon + i * dlon;
    }
    for (i = 0; i < ny; i++){
        glat[i] = slat + i * dlat;
    }

    pn = 1004.0; //outer pressure (hPa)

    // make the output file
    nc_create(E->force_out, NC_CLOBBER, &ncid);
    nc_def_dim(ncid, "lat", ny, &lat_dimid);
    nc_def_dim(ncid, "lon", nx, &lon_dimid);
    nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dimid);

    // def vars
    dimIds[0] = time_dimid;
    dimIds[1] = lat_dimid;
    dimIds[2] = lon_dimid;

    nc_def_var(ncid, "time", NC_DOUBLE, 1, &dimIds[0], &time_varid);
    nc_put_att_text(ncid, time_varid, "long_name", strlen("time, scalar, series"), "time, scalar, series");
    nc_put_att_text(ncid, time_varid, "field", strlen("time"), "time");
    nc_put_att_text(ncid, time_varid, "time", strlen("time"), "time");
    nc_put_att_text(ncid, time_varid, "calendar", strlen("gregorian"), "gregorian");
    nc_put_att_text(ncid, time_varid, "units", strlen("days since 1970-01-01 00:00:0.0"), "days since 1970-01-01 00:00:0.0");

    nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[1], &lat_varid);
    nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[2], &lon_varid);

    
    nc_def_var(ncid, "sustr", NC_DOUBLE, 3, dimIds, &taux_varid);
    nc_put_att_text(ncid, taux_varid, "long_name", strlen("windstress"), "windstress");
    nc_put_att_text(ncid, taux_varid, "units", strlen("N/m^2"), "N/m^2");
    nc_put_att_text(ncid, taux_varid, "field", strlen("sustr, scalar, series"), "sustr, scalar, series");
    nc_put_att_text(ncid, taux_varid, "time", strlen("time"), "time");
    nc_put_att_text(ncid, taux_varid, "coordinates", strlen("lon lat"), "lon lat");
   

    nc_def_var(ncid, "svstr", NC_DOUBLE, 3, dimIds, &tauy_varid);
    nc_put_att_text(ncid, tauy_varid, "long_name", strlen("windstress"), "windstress");
    nc_put_att_text(ncid, tauy_varid, "units", strlen("N/m^2"), "N/m^2");
    nc_put_att_text(ncid, tauy_varid, "field", strlen("svstr, scalar, series"), "svstr, scalar, series");
    nc_put_att_text(ncid, tauy_varid, "time", strlen("time"), "time");
    nc_put_att_text(ncid, tauy_varid, "coordinates", strlen("lon lat"), "lon lat");
    

    nc_def_var(ncid, "Pair", NC_DOUBLE, 3, dimIds, &pressure_varid);
    nc_put_att_text(ncid, pressure_varid, "long_name", strlen("surface pressure"), "surface presure");
    nc_put_att_text(ncid, pressure_varid, "units", strlen("hPa, mbar"), "hPa, mbar");
    nc_put_att_text(ncid, pressure_varid, "field", strlen("Pair, scalar, series"), "Pair, scalar, series");
    nc_put_att_text(ncid, pressure_varid, "time", strlen("time"), "time");
    nc_put_att_text(ncid, pressure_varid, "coordinates", strlen("lon lat"), "lon lat");

    nc_enddef(ncid);

    // write out the lat and lon variables
    nc_put_var_double(ncid, lat_varid, &glat[0]);
    nc_put_var_double(ncid, lon_varid, &glon[0]);

    size_t count[3];
    size_t start[3];
    count[0] = 1;
    count[1] = ny;
    count[2] = nx;
    start[1] = 0;
    start[2] = 0;
    for (fix = 0; fix < E->track_length; fix++){

        printf("------------ fix = %d ----------------------\n", fix);
        lat1 = E->ens_lat[0][fix]; //location of current storm
        lon1 = E->ens_lon[0][fix];

        pc = E->ens_minimumPressure[0][fix];      //central pressure (hPa)
        Rm = E->ens_max_radius[0][fix] * 1.852;   // convert nm to km
        Vm = E->ens_max_speed[0][fix] * 0.514444; // convert knots to m/s

        /*-----------------------------------------------------------------------------*/
        angle_radius(ny, nx, lat1, lon1, glat, glon, ang, radius); //angles and radii of every grid point cf to storm location

        for (i = 0; i < nx * ny; i++)
        {
            if (ang[i] > 180.0)
            {
                printf("angle > 180 %f\n", ang[i]);
            }
            if (ang[i] < -180.0)
            {
                printf("angle < 180  %f\n", ang[i]);
            }
            if (radius[i] > 100000.0)
            {
                printf("radius > 100000m %f\n", radius[i]);
            }
            if (radius[i] < 0.0)
            {
                printf("radius < 0m %f\n", radius[i]);
            }
        }
        /*-----------------------------------------------------------------------------*/

        holland_p(nx, ny, Vm, Rm, pc, pn, rho, radius, p);
        holland_v(nx, ny, Vm, Rm, pc, pn, rho, radius, v);

        cor = coriolis(lat1);
        uv(nx, ny, v, cor, ang, vx, vy);

        stress(nx, ny, p, vx, vy, tau_x, tau_y);

        //rankine_v(nx, ny, Vm, Rm, alpha, radius, rv);

        // write field to netcdf file
        start[0] = fix;
        nc_put_vara_double(ncid, pressure_varid, start, count, &p[0]);
        nc_put_vara_double(ncid, taux_varid, start, count, &tau_x[0]);
        nc_put_vara_double(ncid, tauy_varid, start, count, &tau_y[0]);

        // write time
        nc_put_vara_double(ncid, time_varid, start, count, &E->ens_juldates[0][fix]);

    }

    free(p);
    free(v);
    free(vx);
    free(vy);
    free(tau_x);
    free(tau_y);
    free(rv);
    free(glon);
    free(glat);
    free(ang);
    free(radius);

    nc_close(ncid);

    /*
    // write stress data to file
    // create the file
    nc_create(E->force_out, NC_CLOBBER, &ncid);
    // def dimensions
    nc_def_dim(ncid, "lat", nlat, &lat_dimid);
    nc_def_dim(ncid, "lon", nLon, &lon_dimid);
    nc_def_dim(ncid, "time", nLon, &time_dimid);
    // def vars
    dimIds[0] = time_dimid;
    dimIds[1] = lat_dimid;
    dimIds[2] = lon_dimid;

    nc_def_var(ncid, "time", NC_DOUBLE, 1, &dimIds[0], &lon_varid);
    nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[1], &lat_varid);
    nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[2], &lon_varid);
    
    nc_def_var(ncid, "taux", NC_DOUBLE, 3, dimIds, &taux_varid);
    nc_def_var(ncid, "tauy", NC_DOUBLE, 3, dimIds, &taux_varid);
    nc_def_var(ncid, "pres", NC_DOUBLE, 3, dimIds, &pressure_varid);

    nc_enddef(ncid);
    // write the data
    nc_put_var_double(ncid, time_dimid, &time[0]);
    nc_put_var_double(ncid, lat_varid, &lat[0]);
    nc_put_var_double(ncid, lon_varid, &lon[0]);
    nc_put_var_double(ncid, taux_varid, E->taux);
    nc_put_var_double(ncid, tauy_varid, E->tauy);
    nc_put_var_double(ncid, pressure_varid, E->pressure);

    // close the file
    nc_close(ncid);
    */
    return 0;
}
