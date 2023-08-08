#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
  
  DOM_LOOP(k,j,i){
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

    SetOutputVar("rho", PNG_OUTPUT, YES);
    image = GetImage ("rho");
    image->slice_plane = X13_PLANE;
    image->slice_coord = 0.0;
    image->max = image->min = 0.0;
    image->logscale = 1;
    image->colormap = "red";

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





