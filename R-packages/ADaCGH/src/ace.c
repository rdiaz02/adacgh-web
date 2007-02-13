/* Copyright (C) 2005-2006  Oscar Rueda Palacio <omrueda@cnio.es> */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, */
/* USA. */





#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*************************************************************
Function to perform ACE analysis
Based in ace.java from CGH-Explorer (Ole Lingjaerde)
*************************************************************/

void aceAnalysis(double *x, double *sdev, int *Ngenes, int *Nlevels, 
		 double *alpha1, double *alpha2, double *beta1, 
		 double *beta2, double *ACEPgene, double *FDR,
		 int *calledGenes, int *called, int *first, int *last)
{
  /*counters*/
  int i, j, k, p, q, r;
  double yhat[*Ngenes];
  int nknots;
  int knots[*Ngenes];
  int yhatSign[*Ngenes]; /*boolean*/

  /*step 2 variables*/
  int *z1, Nclusters;
  int *size;
  double *height, *z2;
  /*step 3 variables*/
  double a1, a2, b1, b2,v1, v2;
  /*step 6 variables*/
  double err_opt, e1, e2, e3, v, err;
  int fi, la, p_opt, q_opt;
  
  /******************************************
  Step 1
  Segmentation. Compute running mean
  *******************************************/
  yhat[0] = x[0];
  yhat[1] = x[1];
  for (i=2; i< *Ngenes - 2; i++) {
    yhat[i] = 0;
    for (j= i-2; j <= i +2; j++) {
      yhat[i] += x[j];
    }
    yhat[i] = yhat[i] / 5;

  }
  yhat[*Ngenes -2] = x[*Ngenes -2];
  yhat[*Ngenes -1] = x[*Ngenes -1];

  /******************************************
  Find change points
  *******************************************/
  nknots = 0;

  for (i=0; i< *Ngenes; i++) {
    if (yhat[i] >= 0) {
      yhatSign[i] = 1;
    }
    else {
      yhatSign[i] = 0;
      }
  }
  for (i = 2; i < *Ngenes - 2; i++) {
    if ((yhat[i -2] >0) && (yhat[i - 1] >=0) && (yhat[i + 1] >=0) && (yhat[i + 2] >= 0)) {
      yhatSign[i] = 1;
    }
    else if ((yhat[i -2] < 0) && (yhat[i - 1] < 0) && (yhat[i +1] <0) && (yhat[i +2] <0)) {
      yhatSign[i] = 0;
    }  
  }
  knots[nknots++] = 0;
  for (i = 1; i < *Ngenes; i++) {
    if (yhatSign[i -1] != yhatSign[i]) {
      knots[nknots++] = i;
    }
  }
  knots[nknots++] = *Ngenes;
  
  /**********************************************************
  Step 2
  Feature extraction
  Compute size, height, first, last tuples
  **********************************************************/
  Nclusters = nknots -1;
  int calledC[*Nlevels][Nclusters]; /*boolean. Matrix version of called*/

  size = (int *)malloc(sizeof(int)*Nclusters);
  height = (double *)malloc(sizeof(double)*Nclusters);
  z1 = (int *)malloc(sizeof(int)*Nclusters);
  z2 = (double *)malloc(sizeof(double)*Nclusters);

  for (i=0; i < Nclusters; i++) {
    first[i] = knots[i];
    last[i] = knots[i + 1] - 1;
    size[i] = knots[i + 1] - knots[i];
    height[i] = 0;
    for (j = first[i]; j <= last[i]; j++) {
      height[i] += x[j];
    }
    height[i] = height[i] / size[i];
  }
  for (j=0; j < Nclusters; j++) {
    z1[j] = size[j];
    z2[j] = fabs(height[j]);
  }
  
  /************************************************
  Step 3
  Obtain the null distribution of the (L,H) - pairs
  ************************************************/

  for (i=0; i < *Nlevels; i++) {
    a1 = alpha1[i] * *sdev;
    b1 = beta1[i] * *sdev;
    a2 = alpha2[i] * *sdev;
    b2 = beta2[i] * *sdev;
    for (j=0; j < Nclusters; j++) {
      v1 = z2[j] - (a1 + b1 *z1[j]);
      v2 = z2[j] - (a2 + b2 *z1[j]);
      calledC[i][j] = ! (v1 <= 0 && v2 >=0);
    }
    /*******************************************
    Step 6
    Report genes
    Adjust start/end positions of clusters
    *******************************************/

    if (i == 0) {
      for (j = 0; j < Nclusters; j++) {
	fi = first[j];
	la = last[j];
	err_opt = 9999999;
	p_opt = fi;
	q_opt = la;
        r = (la - fi) / 2;
	if (r > 16) {
	  r = 16;
	}
	r = r - 1;
	if (r <0) {
	  r = 0;
	}
	for (p = fi; p <= fi + r; p++) {
	  for (q = la; q >= la - r; q--) {
	    e1 = 0;
	    e2 = 0;
	    e3 = 0;
            for (k = fi; k < p; k++) {
	      v = x[k];
	      e1 += v * v;
	    }
	    for (k = p; k <= q; k++) {
	      v = x[k] - height[j];
	      e2 += v * v;
	    }
	    for (k = q + 1; k <= la; k++) {
	      v = x[k];
	      e3 += v * v;
	    }
	    err = (e1 + e2 + e3) / (la - fi);
	    if (err < err_opt) {
	      err_opt = err;
	      p_opt = p;
	      q_opt = q;
              }
            } 
          } 
	first[j] = p_opt;
	last[j] = q_opt;
	size[j] = q_opt - p_opt + 1;
      } 
    }
    calledGenes[i] = 0;
    for (j = 0; j < Nclusters; j++) {
      if (calledC[i][j]) {
	calledGenes[i] += last[j] - first[j] + 1;
      }
    }
    FDR[i] = *Ngenes * (1 - ACEPgene[i]) / calledGenes[i];
    }
  /*Convert calledC to matrix from through called*/
  for (i=0; i < Nclusters; i++) {
    for (j=0; j < *Nlevels; j++) {
      called[(i* *Nlevels) + j] =calledC[j][i];
    }
  }
  free(size);
  free(height);
  free(z1);
  free(z2);
 
}

