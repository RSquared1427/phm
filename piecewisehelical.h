class TuningPara {
public:
  double step, large, small; // tuning parameters in HMC for position
  double count;              // count for acceptance

  TuningPara &operator=(const TuningPara &rhs);

  time_t UpdateTuningPara(const vector<int> &acceptance, int i, int size,
                          double lowrate, double highrate);

  void setstep(double ss);

  // constructor and destructor
  TuningPara(void) : step(10), large(10), small(0.000001), count(0) {}
  ~TuningPara() {}
};

time_t TuningPara::UpdateTuningPara(const vector<int> &acceptance,
                                    int current_i, int size, double lowrate,
                                    double highrate) {
  int j;
  double sum = 0;
  for (j = current_i - size + 1; j <= current_i; j++)
    sum += acceptance[j];

  if (sum < lowrate * size) {
    large = step;
    step = (step + small) / 2;
  }
  if (sum > highrate * size) {
    small = step;
    step = (step + large) / 2;
  }
  if (large - small < large * 0.1) {
    large *= 10;
    small *= 0.1;
  }
}

TuningPara &TuningPara::operator=(const TuningPara &rhs) {
  step = rhs.step;
  large = rhs.large;
  small = rhs.small;
  return *this;
}

class component {
public:
  vector<vector<double>> PointPos; // NumPoint*3 coordinates
  vector<double> kappa;            // curvature
  vector<double> tau;              // torsion
  vector<vector<double>> tvec;     // tangent vector, NumPiece*3
  vector<vector<double>> nvec;     // normal vector, NumPiece*3
  vector<vector<double>> bvec;     // binormal vector, NumPiece*3

  vector<double> ArcLen;   // arc length
  vector<int> ChangePoint; // change point for poly-helix, NumPiece+1 (use
                           // Gibbs sampler to update NumPiece-1 change point)

  double beta0;   // intercept
  double beta1;   // slope for 3D distance
  double betaenz; // enzyme effect
  double betagcc; // GC content effect
  double betamap; // mappability effect
  double phiover; // overdisperion
  double loglike; // log likelihood from single component model

  vector<TuningPara>
      TuneKappa; // tune RW MH parametres, acceptance ratio 20% ~ 40%
  vector<TuningPara>
      TuneTau; // tune RW MH parametres, acceptance ratio 20% ~ 40%
  vector<TuningPara> TuneT;
  vector<TuningPara> TuneN;

  TuningPara TuneBeta0;   // tune RW MH parametres, acceptance ratio 20% ~ 40%
  TuningPara TuneBeta1;   // tune RW MH parametres, acceptance ratio 20% ~ 40%
  TuningPara TuneBetaEnz; // tune RW MH parametres, acceptance ratio 20% ~ 40%
  TuningPara TuneBetaGcc; // tune RW MH parametres, acceptance ratio 20% ~ 40%
  TuningPara TuneBetaMap; // tune RW MH parametres, acceptance ratio 20% ~ 40%
  TuningPara TunePhiOver; // tune RW MH parametres, acceptance ratio 20% ~ 40%

  component &operator=(const component &rhs);

  time_t Normalization(void);

  time_t FillPointPos_tnb(void);
  time_t FillPointPos(void);

  time_t FillPointPosKAPPA(int PieceID, double KAPPA,
                           vector<vector<double>> &tvecNew,
                           vector<vector<double>> &nvecNew,
                           vector<vector<double>> &bvecNew,
                           vector<vector<double>> &PointPosNew);
  time_t FillPointPosTAU(int PieceID, double TAU,
                         vector<vector<double>> &tvecNew,
                         vector<vector<double>> &nvecNew,
                         vector<vector<double>> &bvecNew,
                         vector<vector<double>> &PointPosNew);

  // calculate likelihood of the HeatMap based on the single-component model
  double GetLikelihood(const vector<vector<double>> &HeatMap,
                       const vector<vector<double>> &xenz,
                       const vector<vector<double>> &xgcc,
                       const vector<vector<double>> &xmap);

  double GetLikelihood(vector<vector<double>> &PointPosExample,
                       const vector<vector<double>> &HeatMap,
                       const vector<vector<double>> &xenz,
                       const vector<vector<double>> &xgcc,
                       const vector<vector<double>> &xmap);

  // RW MH update kappa
  int UpdateKappa(const vector<vector<double>> &HeatMap,
                  const vector<vector<double>> &xenz,
                  const vector<vector<double>> &xgcc,
                  const vector<vector<double>> &xmap, int PieceID,
                  double stepsize);

  // RW MH update tau
  int UpdateTau(const vector<vector<double>> &HeatMap,
                const vector<vector<double>> &xenz,
                const vector<vector<double>> &xgcc,
                const vector<vector<double>> &xmap, int PieceID,
                double stepsize);

  int UpdateT(const vector<vector<double>> &HeatMap,
              const vector<vector<double>> &xenz,
              const vector<vector<double>> &xgcc,
              const vector<vector<double>> &xmap, int PieceID, double stepsize);

  int UpdateN(const vector<vector<double>> &HeatMap,
              const vector<vector<double>> &xenz,
              const vector<vector<double>> &xgcc,
              const vector<vector<double>> &xmap, int PieceID, double stepsize);

  // RW MH update beta0
  int UpdateBeta0(const vector<vector<double>> &HeatMap,
                  const vector<vector<double>> &xenz,
                  const vector<vector<double>> &xgcc,
                  const vector<vector<double>> &xmap, double stepsize);

  // RW MH update beta1
  int UpdateBeta1(const vector<vector<double>> &HeatMap,
                  const vector<vector<double>> &xenz,
                  const vector<vector<double>> &xgcc,
                  const vector<vector<double>> &xmap, double stepsize);

  // RW MH update betaenz
  int UpdateBetaEnz(const vector<vector<double>> &HeatMap,
                    const vector<vector<double>> &xenz,
                    const vector<vector<double>> &xgcc,
                    const vector<vector<double>> &xmap, double stepsize);

  // RW MH update betagcc
  int UpdateBetaGcc(const vector<vector<double>> &HeatMap,
                    const vector<vector<double>> &xenz,
                    const vector<vector<double>> &xgcc,
                    const vector<vector<double>> &xmap, double stepsize);

  // RW MH update betamap
  int UpdateBetaMap(const vector<vector<double>> &HeatMap,
                    const vector<vector<double>> &xenz,
                    const vector<vector<double>> &xgcc,
                    const vector<vector<double>> &xmap, double stepsize);

  // RW MH update phiover
  int UpdatePhiOver(const vector<vector<double>> &HeatMap,
                    const vector<vector<double>> &xenz,
                    const vector<vector<double>> &xgcc,
                    const vector<vector<double>> &xmap, double stepsize);

  time_t SimCount(vector<vector<double>> &x, const vector<vector<double>> &xenz,
                  const vector<vector<double>> &xgcc,
                  const vector<vector<double>> &xmap);

  time_t GibbsRefine(const vector<vector<double>> &HeatMap,
                     const vector<vector<double>> &xenz,
                     const vector<vector<double>> &xgcc,
                     const vector<vector<double>> &xmap, component &CompMode);

  // constructor and destructor

  component(int NumPoint, int NumPiece, double beta0ini, double beta1ini,
            double betaenzini, double betagccini, double betamapini,
            double phioverini)
      : PointPos(NumPoint, vector<double>(3, 0)),
        ArcLen(vector<double>(NumPoint, 0)), kappa(vector<double>(NumPiece, 0)),
        tau(vector<double>(NumPiece, 0)), tvec(NumPiece, vector<double>(3, 0)),
        nvec(NumPiece, vector<double>(3, 0)),
        bvec(NumPiece, vector<double>(3, 0)),
        ChangePoint(vector<int>(NumPiece + 1, 0)), beta0(beta0ini),
        beta1(beta1ini), TuneKappa(vector<TuningPara>(NumPiece)),
        TuneTau(vector<TuningPara>(NumPiece)),
        TuneT(vector<TuningPara>(NumPiece)),
        TuneN(vector<TuningPara>(NumPiece)), betaenz(betaenzini),
        betagcc(betagccini), betamap(betamapini), phiover(phioverini),
        loglike(0) {}

  ~component() {}
};

component &component::operator=(const component &rhs) {
  PointPos = rhs.PointPos;
  kappa = rhs.kappa;
  tau = rhs.tau;
  tvec = rhs.tvec;
  nvec = rhs.nvec;
  bvec = rhs.bvec;

  ArcLen = rhs.ArcLen;
  ChangePoint = rhs.ChangePoint;

  beta0 = rhs.beta0;
  beta1 = rhs.beta1;
  betaenz = rhs.betaenz;
  betagcc = rhs.betagcc;
  betamap = rhs.betamap;
  phiover = rhs.phiover;
  loglike = rhs.loglike;

  TuneKappa = rhs.TuneKappa;
  TuneTau = rhs.TuneTau;
  TuneT = rhs.TuneT;
  TuneN = rhs.TuneN;

  TuneBeta0 = rhs.TuneBeta0;
  TuneBeta1 = rhs.TuneBeta1;
  TuneBetaEnz = rhs.TuneBetaEnz;
  TuneBetaGcc = rhs.TuneBetaGcc;
  TuneBetaMap = rhs.TuneBetaMap;
  TunePhiOver = rhs.TunePhiOver;

  return *this;
}
