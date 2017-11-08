// normalize weights

void normalize_prob(vector<double> &w) {
  int i;
  double wmax = w[0], sum = 0;

  for (i = 1; i < w.size(); i++) {
    if (w[i] > wmax)
      wmax = w[i];
  }

  for (i = 0; i < w.size(); i++)
    w[i] = exp(w[i] - wmax);
  for (i = 0; i < w.size(); i++)
    sum += w[i];
  for (i = 0; i < w.size(); i++)
    w[i] = w[i] / sum;

  return;
}

// sample from multinomial distribution w

int find_index(vector<double> &w) {
  double prob = gsl_rng_uniform(rng);
  int i = 0;

  while (i + 1 < w.size() && prob > w[i]) {
    prob -= w[i];
    i++;
  }

  return i;
}

// Standardize covariance matrix

void StandardizeMatrix(vector<vector<double>> &x) {
  int i, j;
  int len = x.size();
  double sumx = 0, sumx2 = 0, total = len * (len - 1) * 0.5;
  for (i = 0; i < len - 1; i++) {
    for (j = i + 1; j < len; j++) {
      sumx += x[i][j];
      sumx2 += x[i][j] * x[i][j];
    }
  }

  double mean = sumx / total;
  double var = sumx2 / (total - 1) - mean * mean * total / (total - 1);
  double sd = sqrt(var);

  for (i = 0; i < len - 1; i++) {
    for (j = i + 1; j < len; j++)
      x[i][j] = (x[i][j] - mean) / sd;
  }
}

// log(theta_ij) = beta0 + beta1 * log(d_ij) + betaenz*xenz_ij + betagcc*xgcc_ij
// + betamap*xmap_ij

double d_to_theta(double d, double beta0, double beta1, double betaenz,
                  double betagcc, double betamap, double xenz_ij,
                  double xgcc_ij, double xmap_ij) {
  return exp(beta0 + beta1 * log(d) + betaenz * xenz_ij + betagcc * xgcc_ij +
             betamap * xmap_ij);
}

// Euclidean distance between p and q

double Euclidean_distance(vector<double> &p, vector<double> &q) {
  return sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) +
              (p[2] - q[2]) * (p[2] - q[2]));
}

time_t component::Normalization(void) {
  int i, j;
  double sumlogdij = 0;
  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++)
      sumlogdij += log(Euclidean_distance(PointPos[i], PointPos[j]));
  }

  sumlogdij = sumlogdij / (NumPoint * (NumPoint - 1) / 2);

  beta0 = beta0 + beta1 * sumlogdij;

  for (i = 0; i < NumPiece; i++) {
    kappa[i] = kappa[i] * exp(sumlogdij);
    tau[i] = tau[i] * exp(sumlogdij);
  }

  for (i = 0; i < NumPoint; i++) {
    ArcLen[i] = ArcLen[i] / exp(sumlogdij);
  }

  FillPointPos();
}

void GenerateTNB(vector<vector<double>> &tvecall,
                 vector<vector<double>> &nvecall,
                 vector<vector<double>> &bvecall) {
  int i, j;
  // calculate t, n, b
  // save t, n, b vectors into tvecall, nvecall, bvecall row by row
  for (i = 0; i < NumPiece - 1; i++) {
    // random generate t, n, b vector, then normalize, need 5 parameters in total
    tvecall[i + 1][0] =
        gsl_ran_gaussian(rng, 1); 
    tvecall[i + 1][1] =
        gsl_ran_gaussian(rng, 1); 
    tvecall[i + 1][2] =
        gsl_ran_gaussian(rng, 1); 
    double norm = sqrt(pow(tvecall[i + 1][0], 2) + pow(tvecall[i + 1][1], 2) +
                       pow(tvecall[i + 1][2], 2));
    tvecall[i + 1][0] = tvecall[i + 1][0] / norm;
    tvecall[i + 1][1] = tvecall[i + 1][1] / norm;
    tvecall[i + 1][2] = tvecall[i + 1][2] / norm;

    // random generate the first two coordinates of n vector, then use relationship with t to solve for the last
    nvecall[i + 1][0] =
        gsl_ran_gaussian(rng, 1); 
    nvecall[i + 1][1] =
        gsl_ran_gaussian(rng, 1); 
    nvecall[i + 1][2] = -(tvecall[i + 1][0] * nvecall[i + 1][0] +
                          tvecall[i + 1][1] * nvecall[i + 1][1]) /
                        tvecall[i + 1][2];
    norm = sqrt(pow(nvecall[i + 1][0], 2) + pow(nvecall[i + 1][1], 2) +
                pow(nvecall[i + 1][2], 2));
    nvecall[i + 1][0] = nvecall[i + 1][0] / norm;
    nvecall[i + 1][1] = nvecall[i + 1][1] / norm;
    nvecall[i + 1][2] = nvecall[i + 1][2] / norm;

    // b = cross product of t and  n
    bvecall[i + 1][0] = tvecall[i + 1][1] * nvecall[i + 1][2] -
                        tvecall[i + 1][2] * nvecall[i + 1][1];
    bvecall[i + 1][1] = tvecall[i + 1][2] * nvecall[i + 1][0] -
                        tvecall[i + 1][0] * nvecall[i + 1][2];
    bvecall[i + 1][2] = tvecall[i + 1][0] * nvecall[i + 1][1] -
                        tvecall[i + 1][1] * nvecall[i + 1][0];
    norm = sqrt(pow(bvecall[i + 1][0], 2) + pow(bvecall[i + 1][1], 2) +
                pow(bvecall[i + 1][2], 2));
    bvecall[i + 1][0] = bvecall[i + 1][0] / norm;
    bvecall[i + 1][1] = bvecall[i + 1][1] / norm;
    bvecall[i + 1][2] = bvecall[i + 1][2] / norm;

    if (tvecall[i + 1][0] * nvecall[i + 1][0] +
            tvecall[i + 1][1] * nvecall[i + 1][1] +
            tvecall[i + 1][2] * nvecall[i + 1][2] <
        pow(10, -6))
      cout << "t & n are orthogonal" << endl;
    if (bvecall[i + 1][0] * nvecall[i + 1][0] +
            bvecall[i + 1][1] * nvecall[i + 1][1] +
            bvecall[i + 1][2] * nvecall[i + 1][2] <
        pow(10, -6))
      cout << "b & n are orthogonal" << endl;

    if (tvecall[i + 1][0] * bvecall[i + 1][0] +
            tvecall[i + 1][1] * bvecall[i + 1][1] +
            tvecall[i + 1][2] * bvecall[i + 1][2] <
        pow(10, -6))
      cout << "t & b are orthogonal" << endl;
  }
}

void CalculateR(double kappavalue, double tauvalue, double svalue,
                vector<double> &tvalueini, vector<double> &nvalueini,
                vector<double> &bvalueini, vector<double> &rvalueini,
                vector<double> &output) {
  int i, j;
  output.clear();
  for (i = 0; i < 3; i++)
    output.push_back(0);

  double alphavalue = sqrt(pow(kappavalue, 2) + pow(tauvalue, 2));
  vector<double> V;
  for (i = 0; i < 3; i++)
    V.push_back(0);

  V[0] = (alphavalue * svalue * pow(tauvalue, 2) +
          pow(kappavalue, 2) * sin(alphavalue * svalue)) /
         pow(alphavalue, 3);
  V[1] = kappavalue * (1 - cos(alphavalue * svalue)) / pow(alphavalue, 2);
  V[2] = kappavalue * tauvalue *
         (alphavalue * svalue - sin(alphavalue * svalue)) / pow(alphavalue, 3);

  output[0] = V[0] * tvalueini[0] + V[1] * nvalueini[0] + V[2] * bvalueini[0] +
              rvalueini[0];
  output[1] = V[0] * tvalueini[1] + V[1] * nvalueini[1] + V[2] * bvalueini[1] +
              rvalueini[1];
  output[2] = V[0] * tvalueini[2] + V[1] * nvalueini[2] + V[2] * bvalueini[2] +
              rvalueini[2];

  return;
}

time_t component::FillPointPos_tnb(void) {
  int i, j;

  // initialization

  tvec[0][0] = 0;
  tvec[0][1] = kappa[0] / sqrt(pow(kappa[0], 2) + pow(tau[0], 2));
  tvec[0][2] = tau[0] / sqrt(pow(kappa[0], 2) + pow(tau[0], 2));

  nvec[0][0] = -1;
  nvec[0][1] = 0;
  nvec[0][2] = 0;

  bvec[0][0] = 0;
  bvec[0][1] = -tau[0] / sqrt(pow(kappa[0], 2) + pow(tau[0], 2));
  bvec[0][2] = kappa[0] / sqrt(pow(kappa[0], 2) + pow(tau[0], 2));

  PointPos[0][0] = 0;
  PointPos[0][1] = 0;
  PointPos[0][2] = 0;

  GenerateTNB(tvec, nvec, bvec);

  // calculate r
  for (i = 0; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++)
      CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]], tvec[i],
                 nvec[i], bvec[i], PointPos[ChangePoint[i]], PointPos[j]);
  }
}

time_t component::FillPointPos(void) {
  int i, j;

  // calculate r
  for (i = 0; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++)
      CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]], tvec[i],
                 nvec[i], bvec[i], PointPos[ChangePoint[i]], PointPos[j]);
  }
}

// update kappaNew which is of the PieceID'th helix

time_t component::FillPointPosKAPPA(int PieceID, double kappaNew,
                                    vector<vector<double>> &tvecNew,
                                    vector<vector<double>> &nvecNew,
                                    vector<vector<double>> &bvecNew,
                                    vector<vector<double>> &PointPosNew) {
  int i, j;

  if (PieceID == 0) {
    tvecNew[0][0] = 0;
    tvecNew[0][1] = kappaNew / sqrt(pow(kappaNew, 2) + pow(tau[0], 2));
    tvecNew[0][2] = tau[0] / sqrt(pow(kappaNew, 2) + pow(tau[0], 2));

    nvecNew[0][0] = -1;
    nvecNew[0][1] = 0;
    nvecNew[0][2] = 0;

    bvecNew[0][0] = 0;
    bvecNew[0][1] = -tau[0] / sqrt(pow(kappaNew, 2) + pow(tau[0], 2));
    bvecNew[0][2] = kappaNew / sqrt(pow(kappaNew, 2) + pow(tau[0], 2));

    PointPosNew[0][0] = 0;
    PointPosNew[0][1] = 0;
    PointPosNew[0][2] = 0;
  }

  // calculate r
  for (i = PieceID; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++) {
      if (i == PieceID)
        CalculateR(kappaNew, tau[i], ArcLen[j] - ArcLen[ChangePoint[i]],
                   tvecNew[i], nvecNew[i], bvecNew[i],
                   PointPosNew[ChangePoint[i]], PointPosNew[j]);

      if (i > PieceID)
        CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]],
                   tvecNew[i], nvecNew[i], bvecNew[i],
                   PointPosNew[ChangePoint[i]], PointPosNew[j]);
    }
  }
}

time_t component::FillPointPosTAU(int PieceID, double tauNew,
                                  vector<vector<double>> &tvecNew,
                                  vector<vector<double>> &nvecNew,
                                  vector<vector<double>> &bvecNew,
                                  vector<vector<double>> &PointPosNew) {
  int i, j;

  if (PieceID == 0) {
    tvecNew[0][0] = 0;
    tvecNew[0][1] = kappa[0] / sqrt(pow(kappa[0], 2) + pow(tauNew, 2));
    tvecNew[0][2] = tauNew / sqrt(pow(kappa[0], 2) + pow(tauNew, 2));

    nvecNew[0][0] = -1;
    nvecNew[0][1] = 0;
    nvecNew[0][2] = 0;

    bvecNew[0][0] = 0;
    bvecNew[0][1] = -tauNew / sqrt(pow(kappa[0], 2) + pow(tauNew, 2));
    bvecNew[0][2] = kappa[0] / sqrt(pow(kappa[0], 2) + pow(tauNew, 2));

    PointPosNew[0][0] = 0;
    PointPosNew[0][1] = 0;
    PointPosNew[0][2] = 0;
  }

  // calculate r
  for (i = PieceID; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++) {
      if (i == PieceID)
        CalculateR(kappa[i], tauNew, ArcLen[j] - ArcLen[ChangePoint[i]],
                   tvecNew[i], nvecNew[i], bvecNew[i],
                   PointPosNew[ChangePoint[i]], PointPosNew[j]);

      if (i > PieceID)
        CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]],
                   tvecNew[i], nvecNew[i], bvecNew[i],
                   PointPosNew[ChangePoint[i]], PointPosNew[j]);
    }
  }
}

double component::GetLikelihood(const vector<vector<double>> &HeatMap,
                                const vector<vector<double>> &xenz,
                                const vector<vector<double>> &xgcc,
                                const vector<vector<double>> &xmap) {
  int i, j;

  double sum = 0, d, theta;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      sum += gsl_sf_lngamma(HeatMap[i][j] + phiover) - gsl_sf_lngamma(phiover) +
             phiover * log(phiover) + HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
    }
  }

  return sum;
}

double component::GetLikelihood(vector<vector<double>> &PointPosExample,
                                const vector<vector<double>> &HeatMap,
                                const vector<vector<double>> &xenz,
                                const vector<vector<double>> &xgcc,
                                const vector<vector<double>> &xmap) {
  int i, j;

  double sum = 0, d, theta;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPosExample[i], PointPosExample[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      sum += gsl_sf_lngamma(HeatMap[i][j] + phiover) - gsl_sf_lngamma(phiover) +
             phiover * log(phiover) + HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
    }
  }

  return sum;
}

int component::UpdateKappa(const vector<vector<double>> &HeatMap,
                           const vector<vector<double>> &xenz,
                           const vector<vector<double>> &xgcc,
                           const vector<vector<double>> &xmap, int PieceID,
                           double stepsize) {
  double kappanew = kappa[PieceID] + gsl_ran_gaussian(rng, stepsize);
  if (kappanew < 0)
    return 0;
  //   if (kappanew/tau[PieceID] > 10) return 0;
  if (sqrt(pow(kappanew, 2) + pow(tau[PieceID], 2)) *
          (ArcLen[ChangePoint[PieceID + 1]] - ArcLen[ChangePoint[PieceID]]) /
          (2 * PI) >
      ((ChangePoint[PieceID + 1] - ChangePoint[PieceID]) / 2))
    return 0;

  vector<vector<double>> PointPosNew = PointPos;
  vector<vector<double>> tvecNew = tvec;
  vector<vector<double>> nvecNew = nvec;
  vector<vector<double>> bvecNew = bvec;

  FillPointPosKAPPA(PieceID, kappanew, tvecNew, nvecNew, bvecNew, PointPosNew);

  int i, j;
  double sum = 0, sumnew = 0, dnew, d, thetanew, theta;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      dnew = Euclidean_distance(PointPosNew[i], PointPosNew[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(dnew, beta0, beta1, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);

      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    kappa[PieceID] = kappanew;
    PointPos = PointPosNew;
    tvec = tvecNew;
    nvec = nvecNew;
    bvec = bvecNew;
    TuneKappa[PieceID].count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      kappa[PieceID] = kappanew;
      PointPos = PointPosNew;
      tvec = tvecNew;
      nvec = nvecNew;
      bvec = bvecNew;
      TuneKappa[PieceID].count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateTau(const vector<vector<double>> &HeatMap,
                         const vector<vector<double>> &xenz,
                         const vector<vector<double>> &xgcc,
                         const vector<vector<double>> &xmap, int PieceID,
                         double stepsize) {
  double taunew = tau[PieceID] + gsl_ran_gaussian(rng, stepsize);
  if (taunew < 0)
    return 0;
  //   if (taunew/kappa[PieceID] > 10) return 0;
  if (sqrt(pow(kappa[PieceID], 2) + pow(taunew, 2)) *
          (ArcLen[ChangePoint[PieceID + 1]] - ArcLen[ChangePoint[PieceID]]) /
          (2 * PI) >
      ((ChangePoint[PieceID + 1] - ChangePoint[PieceID]) / 2))
    return 0;

  vector<vector<double>> PointPosNew = PointPos;
  vector<vector<double>> tvecNew = tvec;
  vector<vector<double>> nvecNew = nvec;
  vector<vector<double>> bvecNew = bvec;

  FillPointPosTAU(PieceID, taunew, tvecNew, nvecNew, bvecNew, PointPosNew);

  int i, j;
  double sum = 0, sumnew = 0, d, dnew, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      dnew = Euclidean_distance(PointPosNew[i], PointPosNew[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(dnew, beta0, beta1, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    tau[PieceID] = taunew;
    PointPos = PointPosNew;
    tvec = tvecNew;
    nvec = nvecNew;
    bvec = bvecNew;
    TuneTau[PieceID].count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      tau[PieceID] = taunew;
      PointPos = PointPosNew;
      tvec = tvecNew;
      nvec = nvecNew;
      bvec = bvecNew;
      TuneTau[PieceID].count++;
      return 1;
    } else
      return 0;
  }
}

// b = cross product of t and  n

void GenerateBvec(vector<double> &tvec, vector<double> &nvec,
                  vector<double> &bvec) {
  bvec[0] = tvec[1] * nvec[2] - tvec[2] * nvec[1];
  bvec[1] = tvec[2] * nvec[0] - tvec[0] * nvec[2];
  bvec[2] = tvec[0] * nvec[1] - tvec[1] * nvec[0];
  double norm = sqrt(pow(bvec[0], 2) + pow(bvec[1], 2) + pow(bvec[2], 2));
  bvec[0] = bvec[0] / norm;
  bvec[1] = bvec[1] / norm;
  bvec[2] = bvec[2] / norm;
}

int component::UpdateT(const vector<vector<double>> &HeatMap,
                       const vector<vector<double>> &xenz,
                       const vector<vector<double>> &xgcc,
                       const vector<vector<double>> &xmap, int PieceID,
                       double stepsize) {
  double t1new = tvec[PieceID][0] + gsl_ran_gaussian(rng, stepsize);
  double t2new = tvec[PieceID][1] + gsl_ran_gaussian(rng, stepsize);

  vector<vector<double>> PointPosNew = PointPos;
  vector<vector<double>> tvecNew = tvec;
  vector<vector<double>> nvecNew = nvec;
  vector<vector<double>> bvecNew = bvec;

  tvecNew[PieceID][0] = t1new;
  tvecNew[PieceID][1] = t2new;

  tvecNew[PieceID][2] = -(tvecNew[PieceID][0] * nvecNew[PieceID][0] +
                          tvecNew[PieceID][1] * nvecNew[PieceID][1]) /
                        nvecNew[PieceID][2];
  double norm = sqrt(pow(tvecNew[PieceID][0], 2) + pow(tvecNew[PieceID][1], 2) +
                     pow(tvecNew[PieceID][2], 2));
  tvecNew[PieceID][0] = tvecNew[PieceID][0] / norm;
  tvecNew[PieceID][1] = tvecNew[PieceID][1] / norm;
  tvecNew[PieceID][2] = tvecNew[PieceID][2] / norm;

  GenerateBvec(tvecNew[PieceID], nvecNew[PieceID], bvecNew[PieceID]);

  int i, j;

  for (i = PieceID; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++) {
      CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]],
                 tvecNew[i], nvecNew[i], bvecNew[i],
                 PointPosNew[ChangePoint[i]], PointPosNew[j]);
    }
  }

  double sum = 0, sumnew = 0, d, dnew, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      dnew = Euclidean_distance(PointPosNew[i], PointPosNew[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(dnew, beta0, beta1, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    PointPos = PointPosNew;
    tvec = tvecNew;
    nvec = nvecNew;
    bvec = bvecNew;
    TuneT[PieceID].count++;
    return 1;
  }

  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      PointPos = PointPosNew;
      tvec = tvecNew;
      nvec = nvecNew;
      bvec = bvecNew;
      TuneT[PieceID].count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateN(const vector<vector<double>> &HeatMap,
                       const vector<vector<double>> &xenz,
                       const vector<vector<double>> &xgcc,
                       const vector<vector<double>> &xmap, int PieceID,
                       double stepsize) {
  double n1new = nvec[PieceID][0] + gsl_ran_gaussian(rng, stepsize);
  double n2new = nvec[PieceID][1] + gsl_ran_gaussian(rng, stepsize);
  vector<vector<double>> PointPosNew = PointPos;
  vector<vector<double>> tvecNew = tvec;
  vector<vector<double>> nvecNew = nvec;
  vector<vector<double>> bvecNew = bvec;

  nvecNew[PieceID][0] = n1new;
  nvecNew[PieceID][1] = n2new;
  nvecNew[PieceID][2] = -(tvecNew[PieceID][0] * nvecNew[PieceID][0] +
                          tvecNew[PieceID][1] * nvecNew[PieceID][1]) /
                        tvecNew[PieceID][2];

  double norm = sqrt(pow(nvecNew[PieceID][0], 2) + pow(nvecNew[PieceID][1], 2) +
                     pow(nvecNew[PieceID][2], 2));

  nvecNew[PieceID][0] = nvecNew[PieceID][0] / norm;
  nvecNew[PieceID][1] = nvecNew[PieceID][1] / norm;
  nvecNew[PieceID][2] = nvecNew[PieceID][2] / norm;

  GenerateBvec(tvecNew[PieceID], nvecNew[PieceID], bvecNew[PieceID]);

  int i, j;

  for (i = PieceID; i < NumPiece; i++) {
    for (j = ChangePoint[i] + 1; j <= ChangePoint[i + 1]; j++) {
      CalculateR(kappa[i], tau[i], ArcLen[j] - ArcLen[ChangePoint[i]],
                 tvecNew[i], nvecNew[i], bvecNew[i],
                 PointPosNew[ChangePoint[i]], PointPosNew[j]);
    }
  }

  double sum = 0, sumnew = 0, d, dnew, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      dnew = Euclidean_distance(PointPosNew[i], PointPosNew[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(dnew, beta0, beta1, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }
  if (sumnew >= sum) {
    PointPos = PointPosNew;
    tvec = tvecNew;
    nvec = nvecNew;
    bvec = bvecNew;
    TuneN[PieceID].count++;
    return 1;
  }

  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      PointPos = PointPosNew;
      tvec = tvecNew;
      nvec = nvecNew;
      bvec = bvecNew;
      TuneN[PieceID].count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateBeta0(const vector<vector<double>> &HeatMap,
                           const vector<vector<double>> &xenz,
                           const vector<vector<double>> &xgcc,
                           const vector<vector<double>> &xmap,
                           double stepsize) {
  double beta0new = beta0 + gsl_ran_gaussian(rng, stepsize);
  int i, j;
  double sum = 0, sumnew = 0, d, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);

      thetanew = d_to_theta(d, beta0new, beta1, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    beta0 = beta0new;
    TuneBeta0.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      beta0 = beta0new;
      TuneBeta0.count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateBeta1(const vector<vector<double>> &HeatMap,
                           const vector<vector<double>> &xenz,
                           const vector<vector<double>> &xgcc,
                           const vector<vector<double>> &xmap,
                           double stepsize) {
  double beta1new = beta1 + gsl_ran_gaussian(rng, stepsize);

  int i, j;
  double sum = 0, sumnew = 0, d, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(d, beta0, beta1new, betaenz, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    beta1 = beta1new;
    TuneBeta1.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      beta1 = beta1new;
      TuneBeta1.count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateBetaEnz(const vector<vector<double>> &HeatMap,
                             const vector<vector<double>> &xenz,
                             const vector<vector<double>> &xgcc,
                             const vector<vector<double>> &xmap,
                             double stepsize) {
  double betaenznew = betaenz + gsl_ran_gaussian(rng, stepsize);

  int i, j;

  double sum = 0, sumnew = 0, d, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);

      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);

      thetanew = d_to_theta(d, beta0, beta1, betaenznew, betagcc, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);

      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    betaenz = betaenznew;
    TuneBetaEnz.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      betaenz = betaenznew;
      TuneBetaEnz.count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateBetaGcc(const vector<vector<double>> &HeatMap,
                             const vector<vector<double>> &xenz,
                             const vector<vector<double>> &xgcc,
                             const vector<vector<double>> &xmap,
                             double stepsize) {
  double betagccnew = betagcc + gsl_ran_gaussian(rng, stepsize);
  int i, j;
  double sum = 0, sumnew = 0, d, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);

      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);

      thetanew = d_to_theta(d, beta0, beta1, betaenz, betagccnew, betamap,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);

      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    betagcc = betagccnew;
    TuneBetaGcc.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      betagcc = betagccnew;
      TuneBetaGcc.count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdateBetaMap(const vector<vector<double>> &HeatMap,
                             const vector<vector<double>> &xenz,
                             const vector<vector<double>> &xgcc,
                             const vector<vector<double>> &xmap,
                             double stepsize) {
  double betamapnew = betamap + gsl_ran_gaussian(rng, stepsize);
  int i, j;
  double sum = 0, sumnew = 0, d, theta, thetanew;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);

      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      thetanew = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamapnew,
                            xenz[i][j], xgcc[i][j], xmap[i][j]);
      sum += HeatMap[i][j] * log(theta) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += HeatMap[i][j] * log(thetanew) -
                (HeatMap[i][j] + phiover) * log(phiover + thetanew);
    }
  }

  if (sumnew >= sum) {
    betamap = betamapnew;
    TuneBetaMap.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      betamap = betamapnew;
      TuneBetaMap.count++;
      return 1;
    } else
      return 0;
  }
}

int component::UpdatePhiOver(const vector<vector<double>> &HeatMap,
                             const vector<vector<double>> &xenz,
                             const vector<vector<double>> &xgcc,
                             const vector<vector<double>> &xmap,
                             double stepsize) {
  double phiovernew = phiover + gsl_ran_gaussian(rng, stepsize);
  if (phiovernew < 0)
    phiovernew = phiovernew * (-1); // reflection at 0, make sure phiover > 0
  int i, j;
  double sum = 0, sumnew = 0, d, theta;

  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      sum += gsl_sf_lngamma(HeatMap[i][j] + phiover) - gsl_sf_lngamma(phiover) +
             phiover * log(phiover) -
             (HeatMap[i][j] + phiover) * log(phiover + theta);
      sumnew += gsl_sf_lngamma(HeatMap[i][j] + phiovernew) -
                gsl_sf_lngamma(phiovernew) + phiovernew * log(phiovernew) -
                (HeatMap[i][j] + phiovernew) * log(phiovernew + theta);
    }
  }

  if (sumnew >= sum) {
    phiover = phiovernew;
    TunePhiOver.count++;
    return 1;
  }
  if (sumnew < sum) {
    if (gsl_rng_uniform(rng) < exp(sumnew - sum)) {
      phiover = phiovernew;
      TunePhiOver.count++;
      return 1;
    } else
      return 0;
  }
}

void SummaryStat(const vector<vector<double>> &x, vector<double> &value) {
  value.clear();
  int i, j;

  vector<double> y;
  for (i = 0; i < NumPoint - 1; i++) {
    for (j = i + 1; j < NumPoint; j++)
      y.push_back(x[i][j]);
  }

  vector<double> ycopy(y);
  sort(ycopy.begin(), ycopy.end());

  int size = NumPoint * (NumPoint - 1) / 2;
  value.push_back(ycopy[0]);                  // min
  value.push_back(ycopy[(int)(size * 0.1)]);  // 10% percentile
  value.push_back(ycopy[(int)(size * 0.25)]); // 25% percentile
  value.push_back(ycopy[(int)(size * 0.5)]);  // 50% percentile
  value.push_back(ycopy[(int)(size * 0.75)]); // 75% percentile
  value.push_back(ycopy[(int)(size * 0.9)]);  // 90% percentile
  value.push_back(ycopy[size - 1]);           // max

  double sum = 0, sum2 = 0;
  for (i = 0; i < y.size(); i++) {
    sum += y[i];
    sum2 += y[i] * y[i];
  }
  sum = sum / size; // mean

  value.push_back(sum);                     // mean
  value.push_back(sum2 / size - sum * sum); // var
}

time_t component::SimCount(vector<vector<double>> &x,
                           const vector<vector<double>> &xenz,
                           const vector<vector<double>> &xgcc,
                           const vector<vector<double>> &xmap) {
  x.clear();
  int i, j;
  double d, theta;

  for (i = 0; i < NumPoint; i++) {
    x.push_back(vector<double>());
    for (j = 0; j < NumPoint; j++)
      x[i].push_back(0);
  }

  for (i = 0; i < NumPoint; i++) {
    for (j = i + 1; j < NumPoint; j++) {
      d = Euclidean_distance(PointPos[i], PointPos[j]);
      theta = d_to_theta(d, beta0, beta1, betaenz, betagcc, betamap, xenz[i][j],
                         xgcc[i][j], xmap[i][j]);
      x[i][j] =
          gsl_ran_negative_binomial(rng, phiover / (phiover + theta), phiover);
      x[j][i] = x[i][j];
    }
  }
}

time_t component::GibbsRefine(const vector<vector<double>> &HeatMap,
                              const vector<vector<double>> &xenz,
                              const vector<vector<double>> &xgcc,
                              const vector<vector<double>> &xmap,
                              component &CompMode) {
  int i, j, k;

  vector<vector<int>> UpdateKappaAccept;
  vector<vector<int>> UpdateTauAccept;

  vector<vector<int>> UpdateTAccept;
  vector<vector<int>> UpdateNAccept;

  vector<int> UpdateBeta0Accept;
  vector<int> UpdateBeta1Accept;
  vector<int> UpdateBetaEnzAccept;
  vector<int> UpdateBetaGccAccept;
  vector<int> UpdateBetaMapAccept;
  vector<int> UpdatePhiOverAccept;

  for (i = 0; i < NumPiece; i++) {
    UpdateKappaAccept.push_back(vector<int>());
    UpdateTauAccept.push_back(vector<int>());
    UpdateTAccept.push_back(vector<int>());
    UpdateNAccept.push_back(vector<int>());

    for (j = 0; j < NumGibbs; j++) {
      UpdateKappaAccept[i].push_back(0);
      UpdateTauAccept[i].push_back(0);
      UpdateTAccept[i].push_back(0);
      UpdateNAccept[i].push_back(0);
    }
  }

  for (i = 0; i < NumGibbs; i++) {
    UpdateBeta0Accept.push_back(0);
    UpdateBeta1Accept.push_back(0);
    UpdateBetaEnzAccept.push_back(0);
    UpdateBetaGccAccept.push_back(0);
    UpdateBetaMapAccept.push_back(0);
    UpdatePhiOverAccept.push_back(0);
  }

  vector<double> StatHeatMap;
  SummaryStat(HeatMap, StatHeatMap);

  vector<vector<double>> HeatMapSim;
  vector<double> StatHeatMapSim;
  vector<vector<double>> StatRecord;
  StatRecord.push_back(StatHeatMap);

  FILE *s;
  FILE *sp;
  FILE *sll;
  s = fopen("record.txt", "w");
  sp = fopen("record_p.txt", "w");
  sll = fopen("record_loglike.txt", "w");

  for (i = 0; i < NumGibbs; i++) {
    for (j = 0; j < NumPiece; j++) {
      UpdateKappaAccept[j][i] =
          UpdateKappa(HeatMap, xenz, xgcc, xmap, j, TuneKappa[j].step);
      UpdateTauAccept[j][i] =
          UpdateTau(HeatMap, xenz, xgcc, xmap, j, TuneTau[j].step);
      if (j > 0) {
        UpdateTAccept[j][i] =
            UpdateT(HeatMap, xenz, xgcc, xmap, j, TuneT[j].step);
        UpdateNAccept[j][i] =
            UpdateN(HeatMap, xenz, xgcc, xmap, j, TuneN[j].step);
      }
    }

    UpdateBeta0Accept[i] =
        UpdateBeta0(HeatMap, xenz, xgcc, xmap, TuneBeta0.step);
    UpdateBeta1Accept[i] =
        UpdateBeta1(HeatMap, xenz, xgcc, xmap, TuneBeta1.step);

    UpdateBetaEnzAccept[i] =
        UpdateBetaEnz(HeatMap, xenz, xgcc, xmap, TuneBetaEnz.step);
    UpdateBetaGccAccept[i] =
        UpdateBetaGcc(HeatMap, xenz, xgcc, xmap, TuneBetaGcc.step);
    UpdateBetaMapAccept[i] =
        UpdateBetaMap(HeatMap, xenz, xgcc, xmap, TuneBetaMap.step);
    UpdatePhiOverAccept[i] =
        UpdatePhiOver(HeatMap, xenz, xgcc, xmap, TunePhiOver.step);

    if ((i + 1) % NumTune == 0) {
      for (j = 0; j < NumPiece; j++) {
        TuneKappa[j].UpdateTuningPara(UpdateKappaAccept[j], i, NumTune, 0.2,
                                      0.4);
        TuneTau[j].UpdateTuningPara(UpdateTauAccept[j], i, NumTune, 0.2, 0.4);
        TuneT[j].UpdateTuningPara(UpdateTAccept[j], i, NumTune, 0.2, 0.4);
        TuneN[j].UpdateTuningPara(UpdateNAccept[j], i, NumTune, 0.2, 0.4);
      }

      TuneBeta0.UpdateTuningPara(UpdateBeta0Accept, i, NumTune, 0.2, 0.4);
      TuneBeta1.UpdateTuningPara(UpdateBeta1Accept, i, NumTune, 0.2, 0.4);
      TuneBetaEnz.UpdateTuningPara(UpdateBetaEnzAccept, i, NumTune, 0.2, 0.4);
      TuneBetaGcc.UpdateTuningPara(UpdateBetaGccAccept, i, NumTune, 0.2, 0.4);
      TuneBetaMap.UpdateTuningPara(UpdateBetaMapAccept, i, NumTune, 0.2, 0.4);
      TunePhiOver.UpdateTuningPara(UpdatePhiOverAccept, i, NumTune, 0.2, 0.4);
    }

    Normalization();

    loglike = GetLikelihood(HeatMap, xenz, xgcc, xmap);
    if (loglike > CompMode.loglike)
      CompMode = *this;

    if ((i + 1) % NumTune == 0)
      printf("%d	%.8lf	%.8lf\n", i + 1, loglike, CompMode.loglike);

    if ((i + 1) % NumTune == 0) {
      cout << "kappa, tau" << endl;
      for (j = 0; j < NumPiece; j++)
        printf("%.4lf  %.4lf	  \n", kappa[j], tau[j]);
      printf("%.4lf   %.4lf   %.4lf   %.4lf\n", beta0, beta1, phiover, loglike);
    }

    fprintf(s, "%.4lf   %.4lf   %.4lf   %.4lf   %.4lf   %.4lf\n", beta0, beta1,
            betaenz, betagcc, betamap, phiover);

    fprintf(sll, "%.4lf\n", GetLikelihood(HeatMap, xenz, xgcc, xmap));

    for (j = 0; j < NumPoint; j++) {
      for (k = 0; k < 3; k++)
        fprintf(sp, "%.4lf	", PointPos[j][k]);
    }
    fprintf(sp, "\n");

    // do posterior predictive check
    if (i > NumGibbs * 0.2 && (i + 1) % NumTune == 0) {
      SimCount(HeatMapSim, xenz, xgcc, xmap);
      SummaryStat(HeatMapSim, StatHeatMapSim);
      StatRecord.push_back(StatHeatMapSim);
    }
  }
  fclose(s);
  fclose(sp);
  fclose(sll);

  s = fopen("mode_p.txt", "w");
  for (i = 0; i < NumPoint; i++) {
    for (j = 0; j < 3; j++)
      fprintf(s, "%.4lf	", CompMode.PointPos[i][j]);
    fprintf(s, "\n");
  }
  fclose(s);

  s = fopen("mode_loglike.txt", "w");
  fprintf(s, "mode_loglike = %.6lf\n", CompMode.loglike);
  fclose(s);

  s = fopen("mode_nui.txt", "w");
  fprintf(s, "mode_beta0 = %.6lf\n", CompMode.beta0);
  fprintf(s, "mode_beta1 = %.6lf\n", CompMode.beta1);
  fprintf(s, "mode_betaenz = %.4lf\n", CompMode.betaenz);
  fprintf(s, "mode_betagcc = %.4lf\n", CompMode.betagcc);
  fprintf(s, "mode_betamap = %.4lf\n", CompMode.betamap);
  fprintf(s, "mode_phiover = %.6lf\n", CompMode.phiover);
  fclose(s);

  s = fopen("mode_tnb.txt", "w");
  for (i = 0; i < NumPiece; i++) {
    fprintf(s, "helix%d ", i + 1);
    fprintf(s, "\n");
    fprintf(s, "    t vector");
    for (j = 0; j < 3; j++)
      fprintf(s, "    %.8lf", CompMode.tvec[i][j]);
    fprintf(s, "\n");
    fprintf(s, "    n vector");
    for (j = 0; j < 3; j++)
      fprintf(s, "    %.8lf", CompMode.nvec[i][j]);
    fprintf(s, "\n");
    fprintf(s, "    b vector");
    for (j = 0; j < 3; j++)
      fprintf(s, "    %.8lf", CompMode.bvec[i][j]);
    fprintf(s, "\n");
  }
  fclose(s);

  double totallength = ArcLen[NumPoint - 1];

  int NumPointG = 1000;
  for (i = 1; i < NumPiece + 1; i++)
    CompMode.ChangePoint[i] =
        CompMode.ChangePoint[i - 1] + (int)NumPointG / NumPiece - 1;

  CompMode.ChangePoint[NumPiece] = NumPointG - 1;

  s = fopen("mode_kappa_tau.txt", "w");
  for (i = 0; i < NumPiece; i++)
    fprintf(s, "helix%d	%.8lf	%.8lf\n", i + 1, CompMode.kappa[i],
            CompMode.tau[i]);

  fclose(s);
}
