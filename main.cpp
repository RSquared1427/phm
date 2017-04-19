#include <algorithm>
#include <assert.h>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>

using namespace std;

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl_cdf.h>
#include <gsl_randist.h>
#include <gsl_rng.h>
#include <gsl_sf_exp.h>
#include <gsl_sf_gamma.h>

#define max_length 10000
#define PI 3.1415926

int name_input_heatmap;
int name_input_covariate;
int NumPoint;        // size of heat map
int NumPiece = 2;    // number of helixes within the piecewise helical
int NumGibbs = 5000; // number of gibbs sample
int NumTune = 50;    // number of tuning in gibbs sample
double BetaEnzIni;
double BetaGccIni;
double BetaMapIni;

// define gsl random number generator
gsl_rng *rng;
long seed = 1;

#include "piecewisehelical.h"
#include "piecewisehelical_refine.cpp"

int main(int argc, char **argv) {
  time_t time_start, time_end;
  struct tm *timeinfo_start, *timeinfo_end;

  time(&time_start);
  timeinfo_start = localtime(&time_start);

  FILE *s_time;
  s_time = fopen("time.txt", "w");
  fprintf(s_time, "start time:\t%s\n", asctime(timeinfo_start));

  int i, j, k;
  char *pch = NULL;
  char str[max_length];

  for (i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      name_input_heatmap = i + 1;
      i++;
    } else if (strcmp(argv[i], "-v") == 0) {
      name_input_covariate = i + 1;
      i++;
    } else if (strcmp(argv[i], "-NP") == 0) {
      NumPiece = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-NG") == 0) {
      NumGibbs = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-NT") == 0) {
      NumTune = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-SEED") == 0) {
      seed = atoi(argv[i + 1]);
      i++;
    }
  }

  rng = gsl_rng_alloc(gsl_rng_rand48);
  gsl_rng_set(rng, seed);

  // count the number of points
  i = 0;
  fstream file_count_num(argv[name_input_heatmap], ios::in);
  if (!file_count_num.is_open()) {
    cout << "file " << argv[name_input_heatmap] << " dose not exist!" << endl;
    exit(1);
  }
  while (!file_count_num.eof()) {
    file_count_num.getline(str, max_length);
    i++;
  }
  file_count_num.close();

  NumPoint = i - 1;
  printf("NumPoint = %d\n", NumPoint);

  if (NumPoint <= 3) {
    printf("Require NumPoint >= 4 \n");
    exit(1);
  }

  // Initialize covariate matrix
  vector<vector<double>> covariate;
  for (i = 0; i < NumPoint; i++) {
    covariate.push_back(vector<double>());
    for (j = 0; j < 6; j++)
      covariate[i].push_back(0);
  }

  // begin read covariate
  i = 0;
  j = 0;
  cout << argv[name_input_covariate] << endl;
  fstream file_op_cov(argv[name_input_covariate], ios::in);
  if (!file_op_cov.is_open()) {
    cout << "file " << argv[name_input_covariate] << " does not exist!" << endl;
    exit(1);
  }
  while (!file_op_cov.eof() && i < NumPoint) {
    file_op_cov.getline(str, max_length);
    pch = strtok(str, "	");
    while (pch != NULL && j < 6) {
      covariate[i][j] = atof(pch);
      pch = strtok(NULL, "	");
      j++;
    }
    j = 0;
    i++;
  }
  file_op_cov.close();
  // end read covariate

  vector<double> ArcLen;
  for (i = 0; i < NumPoint; i++)
    ArcLen.push_back(0);
  for (i = 0; i < NumPoint; i++) {
    ArcLen[i] = (covariate[i][1] + covariate[i][2]) / 2 -
                (covariate[0][1] + covariate[0][2]) / 2;
  }
  for (i = 0; i < NumPoint; i++)
    cout << i << "    " << ArcLen[i] << endl;

  double arclenunit = ArcLen[1];
  for (i = 0; i < NumPoint; i++)
    ArcLen[i] = ArcLen[i] / arclenunit;

  for (i = 0; i < NumPoint; i++)
    cout << i << "    " << ArcLen[i] << endl;

  FILE *s;
  s = fopen("ArcLenInitial.txt", "w");
  fprintf(s, "%.6lf	", ArcLen[NumPoint - 1]);
  fclose(s);

  // file in 3 covariate matrix
  vector<vector<double>> xenz;
  vector<vector<double>> xgcc;
  vector<vector<double>> xmap;
  for (i = 0; i < NumPoint; i++) {
    xenz.push_back(vector<double>());
    xgcc.push_back(vector<double>());
    xmap.push_back(vector<double>());
    for (j = 0; j < NumPoint; j++) {
      xenz[i].push_back(0);
      xgcc[i].push_back(0);
      xmap[i].push_back(0);
    }
  }

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      xenz[i][j] = log(covariate[i][3] * covariate[j][3]);
      xgcc[i][j] = log(covariate[i][4] * covariate[j][4]);
      xmap[i][j] = log(covariate[i][5] * covariate[j][5]);
    }
  }

  // standardize xenz, xgcc, xmap
  StandardizeMatrix(xenz);
  StandardizeMatrix(xgcc);
  StandardizeMatrix(xmap);

  // fill in the lower triangle matrix
  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      xenz[j][i] = xenz[i][j];
      xgcc[j][i] = xgcc[i][j];
      xmap[j][i] = xmap[i][j];
    }
  }

  vector<vector<double>> HeatMap;
  for (i = 0; i < NumPoint; i++) {
    HeatMap.push_back(vector<double>());
    for (j = 0; j < NumPoint; j++)
      HeatMap[i].push_back(0);
  }

  // begin read heatmap
  i = 0;
  j = 0;
  fstream file_op(argv[name_input_heatmap], ios::in);
  if (!file_op.is_open()) {
    cout << "file " << argv[name_input_heatmap] << " does not exist!" << endl;
    exit(1);
  }
  while (!file_op.eof() && i < NumPoint) {
    file_op.getline(str, max_length);
    pch = strtok(str, "	");
    while (pch != NULL && j < NumPoint) {
      HeatMap[i][j] = atof(pch);
      pch = strtok(NULL, "	");
      j++;
    }
    j = 0;
    i++;
  }
  file_op.close();
  // end read heatmap

  // output HeatMap[i][j], xenz[i][j], xgcc[i][j], xmap[i][j] for R to run
  // Poisson GLM
  FILE *sglm;
  sglm = fopen("tmp_glm_data.txt", "w");
  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++)
      fprintf(sglm, "%.0lf	%.4lf	%.4lf	%.4lf\n", HeatMap[i][j],
              xenz[i][j], xgcc[i][j], xmap[i][j]);
  }
  fclose(sglm);

  // run R, do Poisson GLM, find initial value for betaenz, betagcc, betamap
  system("R CMD BATCH rcode_glm.txt");

  vector<double> IniPara;
  for (i = 0; i < 3; i++)
    IniPara.push_back(0);

  // begin read initial value
  i = 0;
  fstream file_op_ini("tmp_glm_para.txt", ios::in);
  while (!file_op_ini.eof() && i < 3) {
    file_op_ini.getline(str, max_length);
    IniPara[i] = atof(str);
    i++;
  }
  file_op_ini.close();
  // end read initial value

  BetaEnzIni = IniPara[0];
  BetaGccIni = IniPara[1];
  BetaMapIni = IniPara[2];

  double beta0ini = 5;
  double beta1ini = -1;
  double phioverini = 15;

  component Example(NumPoint, NumPiece, beta0ini, beta1ini, BetaEnzIni,
                    BetaGccIni, BetaMapIni, phioverini);

  for (i = 0; i < NumPoint; i++)
    Example.ArcLen[i] = ArcLen[i];

  for (i = 1; i < NumPiece + 1; i++)
    Example.ChangePoint[i] =
        Example.ChangePoint[i - 1] + (int)NumPoint / NumPiece - 1;
  Example.ChangePoint[NumPiece] = NumPoint - 1;

  for (i = 0; i < Example.ChangePoint.size(); i++)
    printf("%d	%d	%d\n", i, Example.ChangePoint[i],
           (int)NumPoint / NumPiece);

  for (int PieceID = 0; PieceID < NumPiece; PieceID++) {
    double tmp =
        (Example.ChangePoint[PieceID + 1] - Example.ChangePoint[PieceID]) / 4 *
        2 * PI /
        (Example.ArcLen[Example.ChangePoint[PieceID + 1]] -
         Example.ArcLen[Example.ChangePoint[PieceID]]);
    Example.kappa[PieceID] = sqrt(tmp * tmp * 0.5);
    Example.tau[PieceID] = Example.kappa[PieceID];
  }

  Example.FillPointPos_tnb();
  Example.Normalization();

  for (i = 0; i < NumPoint; i++)
    cout << i << "    " << Example.ArcLen[i] << endl;

  printf("loglike = %.8lf\n", Example.GetLikelihood(HeatMap, xenz, xgcc, xmap));

  component CompMode = Example;

  Example.GibbsRefine(HeatMap, xenz, xgcc, xmap, CompMode);

  time(&time_end);
  timeinfo_end = localtime(&time_end);

  fprintf(s_time, "end time:\t%s\n", asctime(timeinfo_end));
  fprintf(s_time, "run time:\t%d\tseconds\n", (int)(time_end - time_start));

  fclose(s_time);

  remove("tmp_glm_data.txt");
  remove("tmp_glm_para.txt");
  remove("rcode_glm.txt.Rout");

  printf("END!\n");
  return 0;
}
