//
//  GrowthInteraction.cpp 
//
//
//  This is code related to interaction and growth integrating both models, GI and SN.
//  Synthesised by Henrike Haebel on 24/04/2018 based on 
//  Arne Pommerening on 11/03/2017. Last modified by him on 21/07/21.
//  Copyright 2018 Philodendron International. All rights reserved.
//

#include <Rcpp.h>
#include <cmath>
#include <string>
using namespace Rcpp;
using namespace std;



/********************************************************************************************/
/* Auxiliaries                                                                              */
/********************************************************************************************/

/*
 * Calculates distance between trees with periodic boundary conditions (Illian et al., 2008, p. 184).
 */
double getEuclideanDistance(double xmax, double ymax, double x1, double y1, double x2, double y2) {
  double dx = fabs(x1 - x2);
  double dy = fabs(y1 - y2);
  dx = min(dx, xmax - dx);
  dy = min(dy, ymax - dy);
  double dz = sqrt(dx * dx + dy * dy);
  return dz;
}

/**
 * Auxiliary function to match species with species index.
 */
int findIndex(NumericVector codes, int speciesi) {
  int z = -1;
  do {
    z++;
  } while(codes[z] != speciesi); 
  // Rcout << "Species index is " << speciesi << " and z is "<< z << std::endl;
  return z;
}



/********************************************************************************************/
/* Interaction                                                                              */
/********************************************************************************************/

/*
 * Calculates interaction for SN regression (parameter estimation).
 */
// [[Rcpp::export]]
NumericVector estTreeInteraction(NumericVector abd, DataFrame xData, bool choice) {
  vector<double> x = as< vector<double> > (xData["x"]);
  vector<double> y = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["dbh"]);
  vector<double> year = as< vector<double> > (xData["year"]);
  // vector<int> plotno = as< vector<int> > (xData["plotno"]); // For BE, NS
  vector<int> plotno = as< vector<int> > (xData["Plot"]); // DF, for everything else
  // vector<string> plotno = as< vector<string> > (xData["Plot"]); // Pinus pinaster
  vector<double> xmax = as< vector<double> > (xData["xmax"]);
  vector<double> ymax = as< vector<double> > (xData["ymax"]);
  int n = mark.size();
  int cc = 0;
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0; 
    cc = 0;
    for (int j = n - 1; j > -1; j--) {
      if ((i != j) && (year[i] == year[j]) && (plotno[i] == plotno[j])) {
        double distance = getEuclideanDistance(xmax[i], ymax[i], x[i], y[i], x[j], y[j]);
        if(choice)
          ci[i] += pow(mark[j] / (1 + 2 * distance), abd[0]) / pow(mark[i], abd[0]);
        else
          ci[i] += (pow(mark[j], abd[0]) * exp(-abd[2] * pow(distance, 1) / pow(mark[j], abd[1]))) / pow(mark[i], abd[0]);        
        cc++;
      }
    }
    // ci[i] /= cc; // Leads to slightly worse results for relascope kernel and does not compute for exponential kernel
    if(choice)
      ci[i] = exp(-abd[1] / ci[i]); 
    else
      ci[i] = exp(-abd[3] / ci[i]); 
  }
  return 1 - ci;
} 

/*
 * Calculates SN field values for a grid.
 */
// [[Rcpp::export]]
NumericVector calcShotNoiseField(NumericVector abdn, NumericVector x, NumericVector y, DataFrame xData, double xmax, double ymax) {
  vector<double> xx = as< vector<double> > (xData["x"]);
  vector<double> yy = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["mark"]);
  int n = x.size();
  int m = xx.size();
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0;
    for (int j = 0; j < m; j++) {
      double distance = getEuclideanDistance(xmax, ymax, x[i], y[i], xx[j], yy[j]);
      ci[i] += pow(mark[j] / (1 + 2 * distance), abdn[0]);
    }
  }
  return ci;
}



/*********************************************************************************************/
/* Relative growth rate RGR                                                                 */
/********************************************************************************************/

/**
 * Mean annual RGR regression.
 */
// [[Rcpp::export]]
NumericVector estimateRelativeGrowthRate(DataFrame xData, NumericVector param, NumericVector abdn, NumericVector species, NumericVector codes, NumericVector pi, bool choice) {
  NumericVector y = as< NumericVector > (xData["dbh"]);
  int n = y.size();
  NumericVector g(n);
  NumericVector ypot(n);
  NumericVector p(n);
  // for (int i = 0; i < n; i++) {
  //   int z = findIndex(codes, species[i]);
  //   dpot = param(0, z) *  param(1, z) * param(2, z) * exp(-param(1, z) * y) * pow(1 - exp(-param(1, z) * y), param(2, z) - 1);
  // }
  // AGR potential
  // ypot = param[0] *  param[1] * param[2] * exp(-param[1] * y) * pow(1 - exp(-param[1] * y), param[2] - 1); // Chapman-Richards
  ypot = param[0] * pow(y, param[1]) * exp(-param[2] *  y); // Modified from Zeide (1993, p. 609, Eq. 7))
  g = ypot * estTreeInteraction(abdn, xData, choice);
  // p = g / y; // When using diameter RGR
  p = g / (pi * pow(y / 200, 2)); // When using basal-area RGR
  // p = -log(1 - p); // Wrong
  p = log(1 + p);
  return p;
}  

/**
 * Mean annual AGR regression.
 */
// [[Rcpp::export]]
NumericVector estimateAbsoluteGrowthRate(DataFrame xData, NumericVector param, NumericVector abdn, NumericVector species, NumericVector codes, NumericVector pi, bool choice) {
  NumericVector y = as< NumericVector > (xData["dbh"]);
  int n = y.size();
  NumericVector g(n);
  NumericVector ypot(n);
  // for (int i = 0; i < n; i++) {
  //   int z = findIndex(codes, species[i]);
  //   dpot = param(0, z) *  param(1, z) * param(2, z) * exp(-param(1, z) * y) * pow(1 - exp(-param(1, z) * y), param(2, z) - 1);
  // }
  // AGR potential
  // ypot = param[0] *  param[1] * param[2] * exp(-param[1] * y) * pow(1 - exp(-param[1] * y), param[2] - 1); // Chapman-Richards
  ypot = param[0] * pow(y, param[1]) * exp(-param[2] *  y); // Modified from Zeide (1993, p. 609, Eq. 7))
  g = ypot * estTreeInteraction(abdn, xData, choice);
  return g; // 200 * sqrt(g / pi);
}  

