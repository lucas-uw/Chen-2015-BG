#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXSTRING 500
#define MAXSTRING2 5000

void ease2latlon(double, double, double, double, double, double *, double *);
void latlon2ease(double, double, double, double, double, double *, double *);

int main (int argc, char *argv[])
/******************************************************************************
  Author: Sasha Richey 2009-Aug-26

  Description:
  This program finds for each grid cell, the fractional area of the wet and
  dry landcover types in the cell, and the fractional areas of the individual
  landcover types within the dry and wet fractions.

  Modifications:

******************************************************************************/
{
  FILE *flist,*ftar,*fpft,*fglwd,*friver,*flai,*fdem;
  int b,c,i,j,k,m,p,r,s,u,v,w;
  int ii_tar,jj_tar,ii_pft,jj_pft,ii_glwd,jj_glwd,ii_river,jj_river,ii_lai,jj_lai,ii_dem,jj_dem;
  double i_offset,j_offset,i_offset_tar,j_offset_tar,i_offset_pft,j_offset_pft,i_offset_glwd,j_offset_glwd,i_offset_river,j_offset_river,i_offset_lai,j_offset_lai,i_offset_dem,j_offset_dem;
  char argstr[10];
  int int_data;
  int coord_prec;
  int data_prec;
  int nrows_tar,ncols_tar,nrows_pft,ncols_pft,nrows_glwd,ncols_glwd,nrows_river,ncols_river,nrows_lai,ncols_lai,nrows_dem,ncols_dem;
  double xllcorner_tar,yllcorner_tar,xllcorner_pft,yllcorner_pft,xllcorner_cell,yllcorner_cell,xllcorner_glwd,yllcorner_glwd,xllcorner_river,yllcorner_river,xllcorner_lai,yllcorner_lai,xllcorner_dem,yllcorner_dem;
  double cellsize_tar,cellsize_pft,cellsize_glwd,cellsize_river,cellsize_lai,cellsize_dem;
  int nodata_tar,nodata_pft,nodata_glwd,nodata_river,nodata_lai,nodata_dem;
  float tmp_float;

 /* for slope calculation*/
  int slope[1000];
  float slopeavg[11];
  int sumslope,sumslope2,nloc;
  int total,surs[11],sure[11],ints[11],inte[11];
  float fenmu,fenzi;
  int flag,max_p,min_p;
  float max_slope,min_slope,slope_slope,b_slope;

  int **laicount,laitcountw,laitcountd;
  float ****laim, **laicum,laitw,laitd,laitwavg,laitdavg,laitavg;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  char filename[MAXSTRING];
  double row_min_cell,row_max_cell,col_min_cell,col_max_cell;
  double lat_min_cell,lat_max_cell,lon_min_cell,lon_max_cell;
  float lat_mid_cell,lon_mid_cell;
  int row_mid_cell,col_mid_cell;
  double pix_size;
  float res_ratio_tar,res_ratio_pft,res_ratio_glwd,res_ratio_river,res_ratio_lai,res_ratio_dem;
  double wet_frac,dry_frac,peat_frac;
  char str1[200];
  char str2[10];
  char str3[10];
  int n_pix_row,n_pix_col;
  double pix_row,pix_col;
  double pix_lat,pix_lon;
  int total_pix,land_pix,peat_tar;
  float **tar_data;
  int **pft_data;
  int **glwd_data;
  int **river_data;
  float **lai_data;
  float **dem_data;
  int **pft_count;
  float **pft_frac;
  int wet_count,dry_count,found_wet_river_alaska;
  int pft_wet_types[] = {0,11};
  int n_pft_wet_types = 2;
  int glwd_lake_type = 1; // Get lakes from GLWD1 and 2
  int found_wet_tar, found_wet_pft, found_wet_glwd, found_wet_river;
  int cell_id[5000];
  int n_cells;
  double lat_cell[5000], lon_cell[5000];
  int row_cell[5000], col_cell[5000];
  char tmpstr3[MAXSTRING2];
  char tmp_lat_cell[15];
  char tmp_lon_cell[15];
  char tmp_lat_cell_soil[15];
  char tmp_lon_cell_soil[15];

  char month[2];
  char laipath[50],lainame[20],laifile[80],lainame2[20];

  double row0,col0;
  double row,col;
  double ease_resolution;
  double lat,lon;
  double minminlon,minmaxlon,maxminlon,maxmaxlon;
  double minminlat,minmaxlat,maxminlat,maxmaxlat;

  /* Usage */
  if(argc!=9) {
    fprintf(stdout,"Usage: %s <cell_list> <tarnocai_file> <glwd1&2> <glwd3> <pft_file> <lai_path> <lai_file> <slope_file>\n",argv[0]);
    fprintf(stdout,"  <cell_list>   Text file containing cell_id, latitude, longitude, row, and column of cells\n");
    fprintf(stdout,"  <tarnocai_file>    Ascii file containing wetland distribution data (float, e.g. presence/absence or type)\n");
    fprintf(stdout,"  <glwd1&2>   Ascii file containing large and smaller lake mask from global lake and wetland database (integer type)\n");
    fprintf(stdout,"  <glwd3>   Ascii file containing river mask from global lake and wetland database (integer type)\n");
    fprintf(stdout,"  <pft_file>   Ascii file containing landcover data (integer, e.g. presence/absence or type)(assumed to be MOD pft product)\n");
    fprintf(stdout,"  <lai_path>   Ascii path of the lai files, before month directory\n");
    fprintf(stdout,"  <lai_file>   Ascii file containing LAI data (float)(assumed to be MOD lai product)\n");
    fprintf(stdout,"  <slope_file>   Ascii file containing slope information(assumed to be GTOPO30)\n");
    exit(0);
  }

  /* Open EASE grid coordinate list */
  if((flist=fopen(argv[1],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[1]);
    exit(1);
  }
  c=0;
  fgets(tmpstr3,MAXSTRING,flist);
  while(!feof(flist)) {
    sscanf(tmpstr3,"%d %d %d %lf %lf", &(cell_id[c]), &(row_cell[c]), &(col_cell[c]), &(lon_cell[c]), &(lat_cell[c]));
    fgets(tmpstr3,MAXSTRING,flist);
    c++;
  }
  fclose(flist);
  n_cells = c;
//  printf("number of cells is %d\n",n_cells);

  /* Open wetland distribution file, i.e. Tarnocai database */
  if((ftar=fopen(argv[2],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[2]);
    exit(1);
  }

  /* Read in the header - tarnocai file */
  fscanf(ftar,"%*s %d",&ncols_tar);
  fscanf(ftar,"%*s %d",&nrows_tar);
  fscanf(ftar,"%*s %s",tmpstr2);
  xllcorner_tar = atof(tmpstr2);
  fscanf(ftar,"%*s %s",tmpstr2);
  yllcorner_tar = atof(tmpstr2);
  fscanf(ftar,"%*s %s",tmpstr2);
  cellsize_tar = atof(tmpstr2);
  fscanf(ftar,"%*s %d",&nodata_tar);

  /* Allocate wet_data array */
  if ( (tar_data = (float**)calloc(nrows_tar,sizeof(float*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for tar_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_tar; i++) {
    if ( (tar_data[i] = (float*)calloc(ncols_tar,sizeof(float))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for tar_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read Tarnocai wetland data into data array */
  for(i=0;i<nrows_tar;i++) {
    for(j=0;j<ncols_tar;j++) {
      fscanf(ftar,"%f",&(tar_data[i][j]));
    }
  } 

  fclose(ftar);

//  printf("finished tarnocai loading\n");

  /* Open GLWD1&2 file */
  if((fglwd=fopen(argv[3],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[3]);
    exit(1);
  }

  /* Read in the header - GLWD lake mask */
  fscanf(fglwd,"%*s %d",&ncols_glwd);
  fscanf(fglwd,"%*s %d",&nrows_glwd);
  fscanf(fglwd,"%*s %s",tmpstr2);
  xllcorner_glwd = atof(tmpstr2);
  fscanf(fglwd,"%*s %s",tmpstr2);
  yllcorner_glwd = atof(tmpstr2);
  fscanf(fglwd,"%*s %s",tmpstr2);
  cellsize_glwd = atof(tmpstr2);
  fscanf(fglwd,"%*s %d",&nodata_glwd);

  /* Allocate glwd1&2 data array */
  if ( (glwd_data = (int**)calloc(nrows_glwd,sizeof(int*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for glwd_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_glwd; i++) {
    if ( (glwd_data[i] = (int*)calloc(ncols_glwd,sizeof(int))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for glwd_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read GLWD1&2 data into data array */
  for(i=0;i<nrows_glwd;i++) {
    for(j=0;j<ncols_glwd;j++) {
      fscanf(fglwd,"%d",&(glwd_data[i][j]));
    }
  }

  fclose(fglwd);
//  printf("finished glwd1&2 loading\n");

  /* Open GLWD3 river mask */
  if((friver=fopen(argv[4],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[4]);
    exit(1);
  }

  /* Read in the header - GLWD river mask */
  fscanf(friver,"%*s %d",&ncols_river);
  fscanf(friver,"%*s %d",&nrows_river);
  fscanf(friver,"%*s %s",tmpstr2);
  xllcorner_river = atof(tmpstr2);
  fscanf(friver,"%*s %s",tmpstr2);
  yllcorner_river = atof(tmpstr2);
  fscanf(friver,"%*s %s",tmpstr2);
  cellsize_river = atof(tmpstr2);
  fscanf(friver,"%*s %d",&nodata_river);

  /* Allocate river data array */
  if ( (river_data = (int**)calloc(nrows_river,sizeof(int*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for river_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_river; i++) {
    if ( (river_data[i] = (int*)calloc(ncols_river,sizeof(int))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for river_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read GLWD river data into data array */
  for(i=0;i<nrows_river;i++) {
    for(j=0;j<ncols_river;j++) {
      fscanf(friver,"%d",&(river_data[i][j]));
    }
  }

  fclose(friver);
//  printf("finished glwd3 loading\n");

  /* Open  MODIS pft file */
  if((fpft=fopen(argv[5],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[5]);
    exit(1);
  }

  /* Read in the header - MODIS pft file */
  fscanf(fpft,"%*s %d",&ncols_pft);
  fscanf(fpft,"%*s %d",&nrows_pft);
  fscanf(fpft,"%*s %s",tmpstr2);
  xllcorner_pft = atof(tmpstr2);
  fscanf(fpft,"%*s %s",tmpstr2);
  yllcorner_pft = atof(tmpstr2);
  fscanf(fpft,"%*s %s",tmpstr2);
  cellsize_pft = atof(tmpstr2);
  fscanf(fpft,"%*s %d",&nodata_pft);

  /* Allocate pft data array */
  if ( (pft_data = (int**)calloc(nrows_pft,sizeof(int*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_pft; i++) {
    if ( (pft_data[i] = (int*)calloc(ncols_pft,sizeof(int))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read MODIS pft data into data array */
  for(i=0;i<nrows_pft;i++) {
    for(j=0;j<ncols_pft;j++) {
      fscanf(fpft,"%d",&(pft_data[i][j]));
    }
  }

  fclose(fpft);
//  printf("finished modis pft loading\n");

  /* Open slope */
  if((fdem=fopen(argv[8],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[8]);
    exit(1);
  }

  /* Read in the header - DEM-slope */
  fscanf(fdem,"%*s %d",&ncols_dem);
  fscanf(fdem,"%*s %d",&nrows_dem);
  fscanf(fdem,"%*s %s",tmpstr2);
  xllcorner_dem = atof(tmpstr2);
  fscanf(fdem,"%*s %s",tmpstr2);
  yllcorner_dem = atof(tmpstr2);
  fscanf(fdem,"%*s %s",tmpstr2);
  cellsize_dem = atof(tmpstr2);
  fscanf(fdem,"%*s %d",&nodata_dem);

  /* Allocate dem-slope data array */
  if ( (dem_data = (float**)calloc(nrows_dem,sizeof(float*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for dem_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_dem; i++) {
    if ( (dem_data[i] = (float*)calloc(ncols_dem,sizeof(float))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for dem_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read slope data into data array */
  for(i=0;i<nrows_dem;i++) {
    for(j=0;j<ncols_dem;j++) {
      fscanf(fdem,"%f",&(dem_data[i][j]));
    }
  }

  fclose(fdem);
//  printf("finished gtopo30 loading\n");

  strcpy(laipath,argv[6]);
  strcpy(lainame,argv[7]);
  strcpy(lainame2,lainame);

  /* Allocate pft_count and pft_frac array */
  if ( (pft_count = (int**)calloc(2,sizeof(int*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_count array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<2; i++) {
    if ( (pft_count[i] = (int*)calloc(26,sizeof(int))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_count array\n",argv[0]);
      exit(1);
    }
  }

  if ( (pft_frac = (float**)calloc(2,sizeof(float*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_frac array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<2; i++) {
    if ( (pft_frac[i] = (float*)calloc(26,sizeof(float))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for pft_frac array\n",argv[0]);
      exit(1);
    }
  }

  /* Allocate laicum and laicount array */
  if ( (laicount = (int**)calloc(2,sizeof(int*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laicount array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<2; i++) {
    if ( (laicount[i] = (int*)calloc(500,sizeof(int))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laicount array\n",argv[0]);
      exit(1);
    }
  }

  if ( (laicum = (float**)calloc(2,sizeof(float*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laicum array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<2; i++) {
    if ( (laicum[i] = (float*)calloc(500,sizeof(float))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laicum array\n",argv[0]);
      exit(1);
    }
  }

  /* Allocate laim array */
  if ( (laim = (float****)calloc(2,sizeof(float***))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laim array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<2; i++) {
    if ( (laim[i] = (float***)calloc(n_cells,sizeof(float**))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laim array\n",argv[0]);
      exit(1);
    }
    for (j=0; j<n_cells; j++) {
      if ( (laim[i][j] = (float**)calloc(26,sizeof(float*))) == NULL ) {
        fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laim array\n",argv[0]);
        exit(1);
      }
      for (k=0; k<26; k++) {
        if ( (laim[i][j][k] = (float*)calloc(12,sizeof(float))) == NULL ) {
          fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for laim array\n",argv[0]);
          exit(1);
        }
      }
    }
  }

 /* Loop over lai files  */
 for(m=0; m<12; m++){

  /* Generate the LAI file name for each month */
  strcpy(lainame,lainame2);
  sprintf(month, "%02d", m+1);
  sprintf(laifile, "%s/%s/%s", laipath,month,lainame);

//  fprintf(stdout,"current laifile is %s\n",laifile);

  /* Open LAI file */
  if((flai=fopen(laifile,"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],laifile);
    exit(1);
  }

  /* Read in the header - MODIS LAI file */
  fscanf(flai,"%*s %d",&ncols_lai);
  fscanf(flai,"%*s %d",&nrows_lai);
  fscanf(flai,"%*s %s",tmpstr2);
  xllcorner_lai = atof(tmpstr2);
  fscanf(flai,"%*s %s",tmpstr2);
  yllcorner_lai = atof(tmpstr2);
  fscanf(flai,"%*s %s",tmpstr2);
  cellsize_lai = atof(tmpstr2);
  fscanf(flai,"%*s %d",&nodata_lai);

  /* Allocate lai_data array */
  if ( (lai_data = (float**)calloc(nrows_lai,sizeof(float*))) == NULL ) {
    fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for lai_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_lai; i++) {
    if ( (lai_data[i] = (float*)calloc(ncols_lai,sizeof(float))) == NULL ) {
      fprintf(stderr,"argv[0]: ERROR: cannot allocate sufficient memory for lai_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Read MODIS LAI data into data array */
  for(i=0;i<nrows_lai;i++) {
    for(j=0;j<ncols_lai;j++) {
      fscanf(flai,"%f",&(lai_data[i][j]));
    }
  } 

  fclose(flai);

//  printf("finished lai loading\n");


  /* Loop over grid cells */
  for (c=0; c<n_cells; c++) {
//    fprintf(stdout,"cell %d row %d col %d lat %lf lon %lf\n",cell_id[c],row_cell[c],col_cell[c],lat_cell[c],lon_cell[c]);

    row_mid_cell = row_cell[c];
    col_mid_cell = col_cell[c];
    lat_mid_cell = lat_cell[c];
    lon_mid_cell = lon_cell[c];

    /* Find EASE grid boundaries of cell */
    row_min_cell = row_mid_cell - 0.5;
    row_max_cell = row_mid_cell + 0.5;
    col_min_cell = col_mid_cell - 0.5;
    col_max_cell = col_mid_cell + 0.5;

    pix_size = 1.0/120.0;   // 30 arc sec pixel resolution in grid cell

    /* Determine lat-lon boundaries of EASE grid domain */
    ease_resolution = 100.2701;
    row0 = 89.625;
    col0 = 89.625;

    ease2latlon(ease_resolution,row0,col0,row_min_cell,col_min_cell,&minminlat,&minminlon);
    ease2latlon(ease_resolution,row0,col0,row_min_cell,col_max_cell,&minmaxlat,&minmaxlon);
    ease2latlon(ease_resolution,row0,col0,row_max_cell,col_min_cell,&maxminlat,&maxminlon);
    ease2latlon(ease_resolution,row0,col0,row_max_cell,col_max_cell,&maxmaxlat,&maxmaxlon);

    lat_min_cell = minminlat;
    if (minmaxlat < lat_min_cell) lat_min_cell = minmaxlat;
    if (maxminlat < lat_min_cell) lat_min_cell = maxminlat;
    if (maxmaxlat < lat_min_cell) lat_min_cell = maxmaxlat;
    lat_max_cell = minminlat;
    if (minmaxlat > lat_max_cell) lat_max_cell = minmaxlat;
    if (maxminlat > lat_max_cell) lat_max_cell = maxminlat;
    if (maxmaxlat > lat_max_cell) lat_max_cell = maxmaxlat;
    lon_min_cell = minminlon;
    if (minmaxlon < lon_min_cell) lon_min_cell = minmaxlon;
    if (maxminlon < lon_min_cell) lon_min_cell = maxminlon;
    if (maxmaxlon < lon_min_cell) lon_min_cell = maxmaxlon;
    lon_max_cell = minminlon;
    if (minmaxlon > lon_max_cell) lon_max_cell = minmaxlon;
    if (maxminlon > lon_max_cell) lon_max_cell = maxminlon;
    if (maxmaxlon > lon_max_cell) lon_max_cell = maxmaxlon;

    /* Set offset to compare 30 arc sec EASE grid pixels to wetland and Bartalev 30 arc sec */
    yllcorner_cell = lat_min_cell;
    xllcorner_cell = lon_min_cell;

    i_offset_tar = (1/cellsize_tar)*(yllcorner_cell-yllcorner_tar);
    j_offset_tar = (1/cellsize_tar)*(xllcorner_cell-xllcorner_tar);
    res_ratio_tar = pix_size/cellsize_tar;

    i_offset_glwd = (1/cellsize_glwd)*(yllcorner_cell-yllcorner_glwd);
    j_offset_glwd = (1/cellsize_glwd)*(xllcorner_cell-xllcorner_glwd);
    res_ratio_glwd = pix_size/cellsize_glwd;

    i_offset_river = (1/cellsize_river)*(yllcorner_cell-yllcorner_river);
    j_offset_river = (1/cellsize_river)*(xllcorner_cell-xllcorner_river);
    res_ratio_river = pix_size/cellsize_river;

    i_offset_pft = (1/cellsize_pft)*(yllcorner_cell-yllcorner_pft);
    j_offset_pft = (1/cellsize_pft)*(xllcorner_cell-xllcorner_pft);
    res_ratio_pft = pix_size/cellsize_pft;

    i_offset_lai = (1/cellsize_lai)*(yllcorner_cell-yllcorner_lai);
    j_offset_lai = (1/cellsize_lai)*(xllcorner_cell-xllcorner_lai);
    res_ratio_lai = pix_size/cellsize_lai;

    i_offset_dem = (1/cellsize_dem)*(yllcorner_cell-yllcorner_dem);
    j_offset_dem = (1/cellsize_dem)*(xllcorner_cell-xllcorner_dem);
    res_ratio_dem = pix_size/cellsize_dem;

    /* Determine number of pixels in the grid cell */
    n_pix_row = (int)( (lat_max_cell-lat_min_cell)/pix_size + 0.5 );
    n_pix_col = (int)( (lon_max_cell-lon_min_cell)/pix_size + 0.5 );

    /* Initialize counts */
    total_pix = 0;
    land_pix = 0;
    wet_count = 0;
    dry_count = 0;
    for(u=0;u<2;u++) {
      for(v=0;v<26;v++) { 
        pft_count[u][v] = 0;
      }
    }
    for(u=0;u<2;u++) {
      for(v=0;v<26;v++) { 
        pft_frac[u][v] = 0;
      }
    }
    peat_tar=0; 
    found_wet_river_alaska=0;
    laitcountw = 0;
    laitcountd = 0;
    laitw = 0;
    laitd = 0;
    laitwavg = 0;
    laitdavg = 0;

    /* Initialization of the lai count and cumulative  */
    for(s=0; s<2; s++) {
      for(p=0; p<26; p++){
        laicum[s][p]=0;
        laicount[s][p]=0;
      }
    }

    /* Initialization of slope counter*/
    for(s=0; s<1000; s++) {
      slope[s]=0;
    }


    /* Loop over pixels in cell lat/lon bounding box */
    for(i=0; i<n_pix_row; i++) {
      pix_lat = yllcorner_cell + (n_pix_row-1-i)*pix_size;
      for(j=0; j<n_pix_col; j++) {
        pix_lon = xllcorner_cell + j*pix_size;


        /* Check if pixel is within EASE grid cell */
        latlon2ease(pix_lat,pix_lon,ease_resolution,row0,col0,&pix_row,&pix_col);
        if((row_mid_cell-0.5) <= pix_row && pix_row < (row_mid_cell+0.5) 
          && (col_mid_cell-0.5) <= pix_col && pix_col < (col_mid_cell+0.5)) {
          total_pix++;

//fprintf(stdout,"pixelrowcol %f %f %d %d %f %f\n",pix_row,pix_col,row_mid_cell,col_mid_cell,pix_lat,pix_lon);

          /* Compute corresponding pixel row/cols in the Tarnocai, glwd, MODIS pft and lai masks */
          ii_tar = (int)( (nrows_tar-1)-i_offset_tar-(res_ratio_tar*(n_pix_row-1-i)) +1 );
          jj_tar = (int)( j_offset_tar + (res_ratio_tar*j) );
          ii_glwd = (int)( (nrows_glwd-1)-i_offset_glwd-(res_ratio_glwd*(n_pix_row-1-i)) +1 );
          jj_glwd = (int)( j_offset_glwd + (res_ratio_glwd*j) );
          ii_river = (int)( (nrows_river-1)-i_offset_river-(res_ratio_river*(n_pix_row-1-i)) +1 );
          jj_river = (int)( j_offset_river + (res_ratio_river*j) );
          ii_pft = (int)( (nrows_pft-1)-i_offset_pft-(res_ratio_pft*(n_pix_row-1-i)) +1 );
          jj_pft = (int)( j_offset_pft + (res_ratio_pft*j) );
          ii_lai = (int)( (nrows_lai-1)-i_offset_lai-(res_ratio_lai*(n_pix_row-1-i)) +1 );
          jj_lai = (int)( j_offset_lai + (res_ratio_lai*j) );
          ii_dem = (int)( (nrows_dem-1)-i_offset_dem-(res_ratio_dem*(n_pix_row-1-i)) +1 );
          jj_dem = (int)( j_offset_dem + (res_ratio_dem*j) );

//fprintf(stdout,"pixlatlons %f %f \n",pix_lat,pix_lon);

          /* Determine wet/dry fractions of cell, and fractions of each landcover type within wet/dry fractions */
          // Is it in the landcover map?
          if(ii_pft>=0 && jj_pft>=0 && ii_pft<nrows_pft && jj_pft<ncols_pft && ii_river>=0 && jj_river>=0 && ii_river<nrows_river && jj_river<ncols_river) {
              if(pft_data[ii_pft][jj_pft]>0 || glwd_data[ii_glwd][jj_glwd]!= nodata_glwd) {
                land_pix++;
  
                
                // Is it a water bodies, wetland in the pft map?
                found_wet_pft = 0;
                for (k=0; k<n_pft_wet_types; k++) {
                  if (pft_data[ii_pft][jj_pft] == pft_wet_types[k]) found_wet_pft = 1;
                }
  
                // Is it a wetland in the Tarnocai wetland map?
                found_wet_tar = 0;
                if(ii_tar>=0 && jj_tar>=0 && ii_tar<nrows_tar && jj_tar<ncols_tar) {
                  if(tar_data[ii_tar][jj_tar] >= 30) {
                    found_wet_tar = 1;
                    peat_tar++;
                  }
                }
  
                // Is it a lake in the GLWD1&2?
                found_wet_glwd = 0;
                if(ii_glwd>=0 && jj_glwd>=0 && ii_glwd<nrows_glwd && jj_glwd<ncols_glwd) {
                  if(glwd_data[ii_glwd][jj_glwd] != nodata_glwd) found_wet_glwd = 1;
                }
  
                // Is it wetland in the GLWD3?
                found_wet_river = 0;
                if(ii_river>=0 && jj_river>=0 && ii_river<nrows_river && jj_river<ncols_river) {
                  if (river_data[ii_river][jj_river] != nodata_river) {
                    if (river_data[ii_river][jj_river] != 10 && river_data[ii_river][jj_river] != 11) { // Get rid of 10(25-50%wetland) and 11(50-100%wetland).20130123
                      found_wet_river = 1;
                    }
                    if (river_data[ii_river][jj_river] == 10) {
                      found_wet_river_alaska++;
                    }
                  }
//                  if(river_data[ii_river][jj_river] != nodata_river ) found_wet_river = 1;
                }
  
                // Decide whether pixel is wet or dry
                if (found_wet_tar || found_wet_pft || found_wet_glwd || found_wet_river) {
                  if (lai_data[ii_lai][jj_lai]!= nodata_lai) {
                    laicum[0][pft_data[ii_pft][jj_pft]] += lai_data[ii_lai][jj_lai];
                    laitw += lai_data[ii_lai][jj_lai];
                    laicount[0][pft_data[ii_pft][jj_pft]]++;
                    laitcountw++;
                  }
                  pft_count[0][pft_data[ii_pft][jj_pft]]++;
                  wet_count++;
                }
                else {
                  if (lai_data[ii_lai][jj_lai]!= nodata_lai) {
                    laicum[1][pft_data[ii_pft][jj_pft]] += lai_data[ii_lai][jj_lai];
                    laitd += lai_data[ii_lai][jj_lai];
                    laicount[1][pft_data[ii_pft][jj_pft]]++;
                    laitcountd++;
                  }
                  pft_count[1][pft_data[ii_pft][jj_pft]]++; 
                  dry_count++;
                }
  
                if(dem_data[ii_dem][jj_dem]>=0 && dem_data[ii_dem][jj_dem]<100){
                  slope[(int)(dem_data[ii_dem][jj_dem]*10)]++;
                }
                if(dem_data[ii_dem][jj_dem]>=100){
                  slope[999]++;
                }
              }

          }
        } // End if pixel is in EASE grid cell
      } // End if pixel is in cell lat/lon bounding box

    } // End loop over pixels
//    fprintf(stdout,"end loop over pixels\n");



//     slope[400]=1;
//     for(p=0; p<400; p++) {
//       slope[p]=0;
//     }
//     for(p=401; p<1000; p++) {
//       slope[p]=0;
//     }
   /* Calculate deciles for slope*/
      nloc = 0;
      ints[1] = 0;
      surs[1] = 0;
      if(m==11){
        total = 0;
        for (s=0; s<1000; s++){
          total = total + slope[s];
        }
// fprintf(stdout,"total dem pixel is %d\n",total);
   //  ints and inte
        for (p=1; p<11; p++){
          sumslope = 0;
          for (s=0; s<1000; s++){
            sumslope2 = sumslope;
            sumslope= sumslope + slope[s];
            if (sumslope2 <= (int)((p-1)*0.1*total) && sumslope > (int)((p-1)*0.1*total)) { ints[p] =s;}
          }
          sumslope = 0;
          for (s=0; s<1000; s++) {
            sumslope2 = sumslope;
            sumslope = sumslope + slope[s];
            if (sumslope2 < (int)(p*0.1*total) && sumslope >= (int)(p*0.1*total)) { inte[p] =s;}
          }
        }
  
      /* loop over 1st decile, a little bit different from rest 9*/
   //  sure
        if (inte[1] == ints[1]) {sure[1] = (int)(0.1*total);}
        if (inte[1] == ints[1]+1) {
          sure[1] = (int)(0.1*total)-slope[0];
        }
        if (inte[1] >= ints[1]+2) {
          sure[1] = (int)(0.1*total);
          for (s=0; s<= inte[1]-1; s++) {
            sure[1] = sure[1]-slope[s];
          }
        }
   //  nloc update
        nloc = sure[1];

        for (p=2; p<11; p++){    /* loops over 9 deciles */
   //  surs
          if (ints[p] == inte[p]) {surs[p] = 0;}
          if (ints[p] != inte[p]) {surs[p] = slope[inte[p-1]]-nloc;}
   //  sure
          if (inte[p] == ints[p]) {sure[p] = (int)(0.1*total);}
          if (inte[p] == ints[p]+1) {sure[p] = (int)(0.1*total)-surs[p];}
          if (inte[p] >= ints[p]+2) {
            sure[p] = (int)(0.1*total);
            for (s= inte[p-1]+1; s<= inte[p]-1; s++){
              sure[p] = sure[p]-slope[s];
            }
            sure[p] = sure[p]-surs[p];
          }
   //  nloc update
          w=0;
          while(sure[p-w] == (int)(0.1*total) && inte[p-w] == inte[p-w-1]) {w++;}
          nloc = 0;
          for (s=0; s<=w; s++) {nloc = nloc+sure[p-s];}
        }

        for (p=1; p<11; p++){
   //  minor correction on surs and sure, it should be improved for perfection
          if (surs[p]<0) {surs[p]=0;}
          if (sure[p]<0) {sure[p]=0;}

//          if(inte[p]==ints[p]) {slopeavg[p]=(float)(ints[p]+1);}
          if(inte[p]==ints[p]) {slopeavg[p]=(float)(ints[p]+1)*0.0001;}
//          if(inte[p]==ints[p]+1) {slopeavg[p]=((float)(surs[p]*(ints[p]+1)+sure[p]*(inte[p]+1)))/(float)(surs[p]+sure[p]);}
          if(inte[p]==ints[p]+1) {slopeavg[p]=((float)(surs[p]*(ints[p]+1)+sure[p]*(inte[p]+1))*0.0001)/(float)(surs[p]+sure[p]);}
          if(inte[p]>=ints[p]+2) {
            fenmu=0;
            fenzi=0;
            for(w=ints[p]+1; w<=inte[p]-1; w++) {
//              fenzi = fenzi + (float)(slope[w]*(w+1));
              fenzi = fenzi + (float)(slope[w]*(w+1))*0.0001;
              fenmu = fenmu + (float)(slope[w]);
            }
//            fenzi = fenzi + (float)(surs[p]*(ints[p]+1)+sure[p]*(inte[p]+1));
            fenzi = fenzi + (float)(surs[p]*(ints[p]+1)+sure[p]*(inte[p]+1))*0.0001;
            fenmu = fenmu + (float)(surs[p]+sure[p]);
            if(fenmu!=0) {slopeavg[p]=fenzi/fenmu;}
            if(fenmu==0) {slopeavg[p]=-9999;}
          }
//          fprintf(stdout,"fenzi=%f,fenmu=%f,slopeavg[%d]=%f\n",fenzi,fenmu,p,slopeavg[p]);
        }
      }
//      fprintf(stdout,"slope calculation is done\n");

     flag = 0;
     max_slope = 0;
     min_slope = 0;
     max_p = 0;
     min_p = 0;
     slope_slope = 0;
     b_slope = 0;
     for(p=1; p<11; p++) {
       if(slopeavg[p]!=-9999) {
       flag++;
         max_slope = slopeavg[p];
         max_p = p;
       }
    }
   if(flag==0) {
       for(p=1; p<11; p++) {
         slopeavg[p] = 0.000001;
       }
   }
     if(flag==1) {
       for(p=1; p<11; p++) {
         slopeavg[p] = max_slope;
       }
     }
     if(flag>1 && flag<10) {
       for(p=10; p>0; p--) {
         if(slopeavg[p]!=-9999) {
           min_slope = slopeavg[p];
           min_p = p;
         }
       }
       slope_slope = (float)(max_slope-min_slope)/(float)(max_p-min_p);
       b_slope = min_slope - (slope_slope)*(float)(min_p);
       for(p=1; p<11; p++) {
         slopeavg[p] = b_slope + (slope_slope)*(float)(p);
       } 
     }

    /* Calculate the averaged LAI for the grid cell  */
    for(p=0; p<26; p++){
      if(laicount[0][p] != 0){
        laim[0][c][p][m] = laicum[0][p]/(float)laicount[0][p];
      }
      if(laicount[0][p] == 0) {
        laim[0][c][p][m] = -9999;
      }
      if(laicount[1][p] != 0){
        laim[1][c][p][m] = laicum[1][p]/(float)laicount[1][p];
      }
      if(laicount[1][p] == 0) {
        laim[1][c][p][m] = -9999;
      }
    }
    /* Calculate averaged LAI for wet and dry parts*/
    if(laitw!=0) laitwavg=laitw/(float)laitcountw;
    if(laitw==0) laitwavg= -9999;
    if(laitd!=0) laitdavg=laitd/(float)laitcountd;
    if(laitd==0) laitdavg= -9999;
    if(laitw!=0 || laitd!=0) laitavg = (laitw+laitd)/((float)laitcountw+(float)laitcountd);
    if(laitw==0 && laitd==0) laitavg =-9999;

//    fprintf(stdout,"lai calculation is finished\n");



    /* Calculate fractions */
    if(m==11){
      if(land_pix > 0) {
        wet_frac = ((float)wet_count+(float)found_wet_river_alaska/2)/(float)land_pix;
//        wet_frac = (float)wet_count/(float)land_pix;
        dry_frac = 1-wet_frac;
        dry_count = land_pix - wet_count;
        peat_frac = (float)peat_tar/(float)wet_count;
        for(v=0;v<26;v++) {
          if(wet_count > 0) {
            pft_frac[0][v] = (float)(pft_count[0][v])/(float)wet_count;
          }
          else {
            pft_frac[0][v] = 0;
          }
          if(dry_count > 0) {
            pft_frac[1][v] = (float)(pft_count[1][v])/(float)dry_count;
          }
          else {
            pft_frac[1][v] = 0;
          }
        }
      }
      else {
        wet_frac = 0;
        dry_frac = 0;
        for(v=0;v<26;v++) {
          pft_frac[0][v] = 0;
          pft_frac[1][v] = 0;
        }
      }
    }
//    fprintf(stdout,"fraction calculation is finished\n");


    /* Print fractions */
   if(m==11) {
     if(land_pix>0) {
       fprintf(stdout,"cell %d %.4f %.4f %d %d",cell_id[c],lat_mid_cell,lon_mid_cell,row_mid_cell,col_mid_cell);
       for(p=1; p<11; p++){
         fprintf(stdout," slopedecile%d %.6f",p,slopeavg[p]);
       }
       fprintf(stdout," totalavglai %f",laitavg);
       fprintf(stdout," wet %f",wet_frac);
       fprintf(stdout," Tar %f",peat_frac);
       fprintf(stdout," wetlai %f",laitwavg);
       for(v=0;v<17;v++) {
         fprintf(stdout," pft_%d %f",v,pft_frac[0][v]);
         for(r=0; r<12; r++) {
           fprintf(stdout," M%d %.08f",r,laim[0][c][v][r]);
         }
       }
       fprintf(stdout," dry %f",dry_frac);
       fprintf(stdout," drylai %f",laitdavg);
       for(v=0;v<17;v++) {
         fprintf(stdout," pft_%d %f",v,pft_frac[1][v]);
         for(r=0; r<12; r++) {
           fprintf(stdout," M%d %.08f",r,laim[1][c][v][r]);
         }
       }
       fprintf(stdout,"  %.05f   ",(float)land_pix/(float)total_pix);
       fprintf(stdout,"\n");
     }
   }
  
  } // End loop over grid cells

  /* Free lai data array */
  for (i=0; i<nrows_lai; i++) {
    free(lai_data[i]);
  }
  free(lai_data);

 }

  /* Free Tarnocai data array */
  for (i=0; i<nrows_tar; i++) {
    free(tar_data[i]);
  }
  free(tar_data);

  /* Free pft data array */
  for (i=0; i<nrows_pft; i++) {
    free(pft_data[i]);
  }
  free(pft_data);

  /* Free glwd1&2 data array */
  for (i=0; i<nrows_glwd; i++) {
    free(glwd_data[i]);
  }
  free(glwd_data);

  /* Free glwd river mask array */
  for (i=0; i<nrows_river; i++) {
    free(river_data[i]);
  }
  free(river_data);

}



void ease2latlon(double res, double row0, double col0, double row, double col, double *lat, double *lon) {

/******************************************************************************
* This function converts latitude and longitude to the corresponding
* row and column in the northern hemisphere EASE-grid mapping system
*
* Inputs
* row  : EASE grid row coordinate (can be a fraction)
* col  : EASE grid col coordinate (can be a fraction)
* res  : EASE grid cell edge length (km) (common values are 25.067525 and 100.2701)
* row0 : Row index of North Pole
* col0 : Column index of North Pole
* North pole is typically located at (360.0,360.0) (25km cell length)
*                                 or (89.625,89.625) (100km cell length)
* Fractional row/col values indicate fraction of distance from center of this EASE row/col
* to the next row/col.  For example, (row,col) = (1.0, 2.0) is at the center of the
* grid cell on row 1 and col 2.  (row,col) = (1.5, 2.0) is on the border between
* row 1 and 2.
*
* Outputs
* lat  : latitude (deg) corresponding to (fractional) row/col position
* lon  : longitude (deg) corresponding to (fractional) row/col position
******************************************************************************/

  // Constants
  double R = 6371.228; // Radius of Earth (km)
  double PI = 3.1415927;

  // Other variables
  double tmp;
  double phi; // latitude in radians
  double lambda; // longitude in radians

  // Compute lambda and phi
  tmp = 0.5*res/R*sqrt((col-col0)*(col-col0) + (row-row0)*(row-row0));
  phi = 2*(PI/4 - atan2(tmp, sqrt(1-tmp*tmp)));
  lambda = atan2((col-col0),(row-row0));

  // Convert rad to deg
  *lat = phi*180/PI;
  *lon = lambda*180/PI;

}




void latlon2ease(double lat, double lon, double res, double row0, double col0, double *row, double *col) {

/******************************************************************************
* This function converts latitude and longitude to the corresponding
* row and column in the northern hemisphere EASE-grid mapping system
*
* Inputs
* lat  : latitude (deg)
* lon  : longitude (deg)
* res  : EASE grid cell edge length (km) (common values are 25.067525 and 100.2701)
* row0 : Row index of North Pole
* col0 : Column index of North Pole
* North pole is typically located at (360.0,360.0) (25km cell length)
*                                 or (89.625,89.625) (100km cell length)
*
* Outputs
* row  : EASE grid row coordinate (can be a fraction)
* col  : EASE grid col coordinate (can be a fraction)
* Fractional values indicate fraction of distance from center of this EASE row/col
* to the next row/col.  For example, (row,col) = (1.0, 2.0) is at the center of the
* grid cell on row 1 and col 2.  (row,col) = (1.5, 2.0) is on the border between
* row 1 and 2.
******************************************************************************/

  // Constants
  double R = 6371.228; // Radius of Earth (km)
  double PI = 3.1415927;

  // Other variables
  double phi; // latitude in radians
  double lambda; // longitude in radians

  // Convert deg to rad
  phi = lat*PI/180;
  lambda = lon*PI/180;

  // Convert lat and lon to EASE-grid row and col
  *col = 2*R/res * sin(lambda) * sin(PI/4 - phi/2) + col0;
  *row = 2*R/res * cos(lambda) * sin(PI/4 - phi/2) + row0;

}

   
