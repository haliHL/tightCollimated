//This program is used to simulate the tight collimated beam laser
//  --using the cNumber method 
//  --using the Adiabatic Elimination Approximation
#include "tightCollimated.hpp"
#include "config.hpp"

//Changes required subject to the definition of Param 
void getParam(const char* filename, Param *param) 
{
  std::ifstream configInput(filename);
  std::string dummy;

  while (!configInput.eof()) {
    configInput >> dummy;
    if (configInput.eof()) break;
    if (!configInput.good()) {
      std::cout << "Bad read in input file" << std::endl;
      exit(-1);
    }
    if (dummy.compare("dt") == 0)
      configInput >> param->dt;
    else if (dummy.compare("tMax") == 0)
      configInput >> param->tMax;
    else if (dummy.compare("nStore") == 0)
      configInput >> param->nStore;
    else if (dummy.compare("nTrajectory") == 0)
      configInput >> param->nTrajectory;
    else if (dummy.compare("nBin") == 0)
      configInput >> param->nBin;
    else if (dummy.compare("yWall") == 0)
      configInput >> param->yWall;
    else if (dummy.compare("lambda") == 0)
      configInput >> param->lambda;
    else if (dummy.compare("deltaZ") == 0)
      configInput >> param->deltaZ;
    else if (dummy.compare("deltaPz") == 0)
      configInput >> param->deltaPz;
    else if (dummy.compare("tau") == 0)
      configInput >> param->tau;
    else if (dummy.compare("density") == 0)
      configInput >> param->density;
    else if (dummy.compare("mAtom") == 0)
      configInput >> param->mAtom;
    else if (dummy.compare("rabi") == 0)
      configInput >> param->rabi;
    else if (dummy.compare("kappa") == 0)
      configInput >> param->kappa;
    else if (dummy.compare("detuning") == 0)
      configInput >> param->detuning;
    else if (dummy.compare("T2") == 0)
      configInput >> param->T2;
    else if (dummy.compare("controlType") == 0)
      configInput >> param->controlType;    
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else if (dummy.compare("pois") == 0)
      configInput >> param->pois;
    else if (dummy.compare("fast") == 0)
      configInput >> param->fast;
    else {
      std::cout << "Error: invalid label " << dummy << " in "
          << filename << std::endl;
      exit(-1);
    }
  }
}

// void generateInitialField(Ensemble& ensemble, const Param& param) 
// {
//   int nTimeStep = param.tMax/param.dt+0.5;
//   ensemble.cavity.Jx.setZero(param.nTrajectory, nTimeStep+1);
//   ensemble.cavity.Jy.setZero(param.nTrajectory, nTimeStep+1);
//   ensemble.cavity.Jz.setZero(param.nTrajectory, nTimeStep+1);
// }

// void generateInitialSpinFinal(Ensemble& ensemble, const Param& param) 
// {
//   ensemble.spinFinalVector.sxFinalVector.setZero(param.nTrajectory);
//   ensemble.spinFinalVector.syFinalVector.setZero(param.nTrajectory);
//   ensemble.spinFinalVector.szFinalVector.setZero(param.nTrajectory);
// }

void getGammaValues(double& gc, double& gd, double& g0, const Param& param)
{
  const double rabi = param.rabi;
  const double kappa = param.kappa;
  const double detuning = param.detuning;
  gc = (pow(rabi, 2) * kappa / 4) / (pow(kappa, 2) / 4 + pow(detuning, 2));
  gd = (pow(rabi, 2) * detuning / 2) / (pow(kappa, 2) / 4 + pow(detuning, 2));
  g0 = pow(rabi, 2) / kappa;
}

void getRabimult(VectorXd& rabiMult, const Ensemble& ensemble, const Param& param)
{
  //nAtom
  const int nAtom = rabiMult.size();
  //Define a wavenumber k
  double k = 2 * M_PI / param.lambda;
  //Define a beam waist;
  double waistRatio = 999999;
  double y0 = waistRatio * param.yWall;

  for (int i = 0; i < nAtom; i++) 
  {
    //z direction, standing wave
    rabiMult[i] = cos(k * ensemble.atoms[i].external.X[2]);
    //y direction, gaussian
    rabiMult[i] *= exp(- pow(ensemble.atoms[i].external.X[1], 2) / pow(y0, 2)); 
  }
}


void generateExternalState(Atom& newAtom, const Param& param)
{
  //meanP
  double meanP = param.yWall * 2/ param.tau; //vy = deltay/tau//Initial position

  //no distribution
  // Vector3d X (0, -param.yWall, 0);
  //with distribution
  double deltaZ = param.deltaZ;
  Vector3d X (0, -param.yWall, rng.get_uniform_rn(-deltaZ, deltaZ));


  //Initial velocity

  //no doppler
  //Vector3d P (0, meanP, 0);

  //doppler
  Vector3d P (0, meanP, rng.get_gaussian_rn(param.deltaPz));

  //linear
  // Vector3d P (0, meanP, param.deltaPz);

  //biased doppler?????
  // Vector3d P (0, meanP, param.deltaPz + rng.get_gaussian_rn(param.deltaPz));

  //Complete initiation
  External newExternal = {X,P};
  newAtom.external = newExternal;
}

void generateInternalState(Atom& newAtom, const Param& param)
{ 
  //For convenience
  const int nTrajectory = param.nTrajectory;

  //Set up empty internal state vectorsm
  VectorXd newSx, newSy, newSz;
  newSx.setZero(nTrajectory);
  newSy.setZero(nTrajectory);
  newSz.setOnes(nTrajectory);
  newSz *= param.mAtom;
  
  //Random initialization for sx and sy.
  //+1, -1 approach
  for (int j = 0; j < nTrajectory; j++)
  {
    for (int k = 0; k < param.mAtom; k++)
    {
      newSx(j) += double(rng.get_binomial_int(0.5, 1)) * 2 - 1;//50percent giving 1 or -1
      newSy(j) += double(rng.get_binomial_int(0.5, 1)) * 2 - 1;//50percent giving 1 or -1
    }
  } 
  // random phase approach
  // for (int j = 0; j < nTrajectory; j++) 
  // {
  //   double phi = rng.get_uniform_rn(0, 2 * M_PI);
  //   newSx[j] = sqrt(2) * cos(phi); //
  //   newSy[j] = sqrt(2) * sin(phi); //
  // }

  //Complete initiation
  Internal newInternal = {newSx, newSy, newSz};
  newAtom.internal = newInternal;
}

void addAtomsFromVanillaSource(Ensemble& ensemble, const Param& param)
{
  int nAtom;
  const double dN = param.density*param.dt;

  ///////////////////////////////////////////////////////////////////////////
  // Uniform atom generation
  ///////////////////////////////////////////////////////////////////////////
  // std::cout << "Uniform generation turned ON." << std::endl << std::endl;
  if (dN >= 1) {
    nAtom = dN;     

    for (unsigned long int n = 0; n < nAtom; n++) 
    {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
    }
  }
  else 
  {
    std::cout << "dN < 1 in input file." << std::endl;
    exit(-1);
  }
}
  
void addAtomsFromPoissonianSource(Ensemble& ensemble, const Param& param)
{
  int nAtom;
  const double dN = param.density * param.dt;
///////////////////////////////////////////////////////////////////////////
  //Poissonian atom generation
  ///////////////////////////////////////////////////////////////////////////
  //// std::cout << "Poissionian generation turned ON." << std::endl << std::endl;
  nAtom = rng.get_poissonian_int(dN);
  for (int n = 0; n < nAtom; n++) 
  {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
      // debug
      // if (abs(newAtom.external.P[2]) <= 0.5*param.tau/param.lambda)
      //   ensemble.atoms.push_back(newAtom);
      // debug
  }
}

void removeAtomsAtWalls(Ensemble& ensemble, const Param& param) 
{
  // generateInitialSpinFinal(ensemble, param);
  std::vector<Atom> newAtoms;

  // int nLeaving = 0;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++) 
  {
    if (a->external.X[1] < param.yWall) 
    {
      newAtoms.push_back(*a);
    } 
    // else 
    // {
    //   nLeaving += 1;
    //   ensemble.spinFinalVector.sxFinalVector += a->internal.sx;
    //   ensemble.spinFinalVector.syFinalVector += a->internal.sy;
    //   ensemble.spinFinalVector.szFinalVector += a->internal.sz;
    // }
  }
  ensemble.atoms = newAtoms;
  // ensemble.spinFinalVector.sxFinalVector /= nLeaving;
  // ensemble.spinFinalVector.syFinalVector /= nLeaving;
  // ensemble.spinFinalVector.szFinalVector /= nLeaving;
}

void advanceExternalStateOneTimeStep(Ensemble& ensemble, const Param& param) 
{
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    a->external.X += param.dt * a->external.P;
}

void getDriftVector(VectorXd& drift, const VectorXd& s_total, const VectorXd& rabiMult, const Param& param) 
{
  //For convenience
  const int nAtom = rabiMult.size();
  double gc, gd, g0;
  getGammaValues(gc, gd, g0, param);

  //Define Jx and Jy
  double Jx = 0, Jy = 0;
  for (int i = 0; i < nAtom; i++)
  {
    Jx += s_total(3 * i) * rabiMult(i);
    Jy += s_total(3 * i + 1) * rabiMult(i);
  }

  //Drift vector terms. 
   for (int i = 0; i < nAtom; i++) 
   {
      double c_i = rabiMult(i);
      double sx_i = s_total(3 * i);
      double sy_i = s_total(3 * i + 1);
      double sz_i = s_total(3 * i + 2);
      drift(3 * i) = gc / 2 * c_i * (Jx * sz_i - c_i * sx_i * (sz_i + 1))
                   - gd / 2 * c_i * (Jy * sz_i - c_i * sy_i * (sz_i + 1));// - invT2*sx(i);
      drift(3 * i + 1) = gc / 2 * c_i * (Jy * sz_i - c_i * sy_i * (sz_i + 1))
                       - gd / 2 * c_i * (Jx * sz_i - c_i * sx_i * (sz_i + 1));// - invT2*sy(i);
      drift(3 * i + 2) = - gc * pow(c_i, 2) * (sz_i + 1)
                         - gc / 2 * c_i * (Jx * sx_i + Jy * sy_i - c_i * (pow(sx_i, 2) + pow(sy_i, 2)))
                         + gd / 2 * c_i * (Jy * sx_i - Jx * sy_i);
   }
}

void getDiffusionMatrix(Matrix<double, Dynamic, 2>& diffusionMatrix, const VectorXd& s_total, const VectorXd& rabiMult, const Param& param)
{
  //Useful parameters
  int nAtom = rabiMult.size();
  double gc, gd, g0;
  getGammaValues(gc, gd, g0, param);

  // diffusion for the atoms
  for (int i = 0; i < nAtom; i++) 
  {
    double c_i = rabiMult[i];
    double sx_i = s_total(3 * i);
    double sy_i = s_total(3 * i + 1);
    double sz_i = s_total(3 * i + 2);
    diffusionMatrix(3 * i, 0) = - c_i / sqrt(g0) * sz_i * gd;
    diffusionMatrix(3 * i, 1) = - c_i / sqrt(g0) * sz_i * gc;
    diffusionMatrix(3 * i + 1, 0) = c_i / sqrt(g0) * sz_i * gc;
    diffusionMatrix(3 * i + 1, 1) = - c_i / sqrt(g0) * sz_i * gd;
    diffusionMatrix(3 * i + 2, 0) = c_i / sqrt(g0) * (sx_i * gd - sy_i * gc);
    diffusionMatrix(3 * i + 2, 1) = c_i / sqrt(g0) * (sx_i * gc + sy_i * gd);
  }
}

void stochasticIntegration_order_1(VectorXd& s_total, const VectorXd& rabiMult, const Param& param)
{
  //Useful parameters
  const int nAtom = rabiMult.size();

  //Drift
  VectorXd drift;
  drift.setZero(3 * nAtom);
  getDriftVector(drift, s_total, rabiMult, param);

  //Diffusion
  //Diffusion matrix
  Matrix<double, Dynamic, 2> diffusionMatrix;
  diffusionMatrix.setZero(3 * nAtom, 2);
  getDiffusionMatrix(diffusionMatrix, s_total, rabiMult, param);
  //Define the 2 wiener processes
  Vector2d dW_wiener;
  for (int i = 0; i < 2; i++) 
  {
    dW_wiener(i) = rng.get_gaussian_rn(sqrt(param.dt));
  }
  //Diffusion vector
  VectorXd dW;
  dW.setZero(3 * nAtom);
  dW = diffusionMatrix * dW_wiener;

  //1st order
  s_total += drift * param.dt + dW;
}

void stochasticIntegration_order_2(VectorXd& s_total, const VectorXd& rabiMult, const Param& param)
{
  //Useful parameters
  const int nAtom = rabiMult.size();
  const double dt = param.dt;

  //===============================================================================================================================================
  //=========================================================stochasticIntegration_order_1=========================================================
  //Drift
  VectorXd drift_1;
  drift_1.setZero(3 * nAtom);
  getDriftVector(drift_1, s_total, rabiMult, param);

  //Diffusion
  //Diffusion matrix
  Matrix<double, Dynamic, 2> diff_1;
  diff_1.setZero(3 * nAtom, 2);
  getDiffusionMatrix(diff_1, s_total, rabiMult, param);
  //Define the 2 wiener processes
  Vector2d dW_wiener;
  for (int i = 0; i < 2; i++) 
  {
    dW_wiener(i) = rng.get_gaussian_rn(sqrt(dt));
  }

  //1st order
  VectorXd s_total_temp;
  s_total_temp = s_total + drift_1 * dt + diff_1 * dW_wiener;
  //=========================================================stochasticIntegration_order_1=========================================================
  //===============================================================================================================================================
  
  //Drift_temp=====================================================================================================================================
  VectorXd drift_temp;
  drift_temp.setZero(3 * nAtom);
  getDriftVector(drift_temp, s_total_temp, rabiMult, param);
  //Define drift term
  VectorXd drift_term;
  drift_term = 0.5 * (drift_temp + drift_1) * dt;

  //Diffusion_temp=================================================================================================================================
  //Set up R and U
  Matrix<double, Dynamic, 2> R_p, R_m, U_p, U_m;
  R_p.setZero(3 * nAtom, 2);
  R_m.setZero(3 * nAtom, 2);
  U_p.setZero(3 * nAtom, 2);
  U_m.setZero(3 * nAtom, 2);
  for (int i = 0; i < 2; i++)
  { 
    //R_p and R_m
    R_p.col(i) = s_total + drift_1 * dt + diff_1.col(i) * sqrt(dt);
    R_m.col(i) = s_total + drift_1 * dt - diff_1.col(i) * sqrt(dt);
    //U_p and U_m
    U_p.col(i) = s_total + diff_1.col(i) * sqrt(dt);
    U_m.col(i) = s_total - diff_1.col(i) * sqrt(dt);
  }

  //Set up V(i, j)
  Matrix<double, 2, 2> V;
  for (int i = 0; i < 2; i++)
  {
    //Diagonal
    V(i, i) = - dt;
    //Off-diagonal
    for (int j = 0; j < i; j++)
    {
      V(i, j) = dt * (double(rng.get_binomial_int(0.5, 1)) * 2 - 1); //50% to get dt or - dt
    }
  }
  //The other half
  for (int i = 0; i < 2; i++)
  {
    for (int j = i + 1; j < 2; j++)
    {
      V(i, j) = - V(j, i);
    }
  }

  //Define diffusion term
  VectorXd diffusion_term;
  diffusion_term.setZero(3 * nAtom);

  //Get diffusion by double loop
  for (int i = 0; i < 2; i++)
  { 
    //diff_R_p, diff_R_m;
    Matrix<double, Dynamic, 2> diff_R_p, diff_R_m;
    diff_R_p.setZero(3 * nAtom, 2);
    diff_R_m.setZero(3 * nAtom, 2);
    getDiffusionMatrix(diff_R_p, R_p.col(i), rabiMult, param);
    getDiffusionMatrix(diff_R_m, R_m.col(i), rabiMult, param);
    //Diffusion term increment
    diffusion_term += 0.25 * (diff_R_p.col(i) + diff_R_m.col(i) + 2 * diff_1.col(i)) * dW_wiener(i)
                    + 0.25 * (diff_R_p.col(i) - diff_R_m.col(i)) * (pow(dW_wiener(i), 2) - dt) / sqrt(dt);
    for(int j = 0; j < 2; j++)
    {
      if (j != i)
      {
        //diff_U_p, diff_U_m;
        Matrix<double, Dynamic, 2> diff_U_p, diff_U_m;
        diff_U_p.setZero(3 * nAtom, 2);
        diff_U_m.setZero(3 * nAtom, 2);
        getDiffusionMatrix(diff_U_p, R_p.col(j), rabiMult, param);
        getDiffusionMatrix(diff_U_m, R_m.col(j), rabiMult, param);
        //Diffusion term increment
        diffusion_term += 0.25 * (diff_U_p.col(i) + diff_U_m.col(i) - 2 * diff_1.col(i)) * dW_wiener(i)
                        + 0.25 * (diff_U_p.col(i) - diff_U_m.col(i)) * (dW_wiener(i) * dW_wiener(j) + V(i, j)) / sqrt(dt);
      }
    }
  }

  //2nd order
  s_total += drift_term + diffusion_term;
}

void stochasticIntegration(VectorXd& s_total, const VectorXd& rabiMult, const Param& param, const int& num_order) 
{
  if (num_order == 1)
  {
    //Using Euler's method to integrate the stochastic differential equations.
    stochasticIntegration_order_1(s_total, rabiMult, param);
  }
  else if (num_order == 2)
  {
    //Using weak RK-2 method to integrate the stochastic differential equations.
    stochasticIntegration_order_2(s_total, rabiMult, param);
  }
  else//Here allows higher orders
  {
    std::cout << "Wrong input for order of stochastic integration." << std::endl;
    exit(-1);
  } 
}

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param)
{
  //Useful parameters
  const int nAtom = ensemble.atoms.size();

  //effective rabi
  VectorXd rabiMult;
  rabiMult.setZero(nAtom);
  getRabimult(rabiMult, ensemble, param);

  //Loop over all trajectories.
  for (int n = 0; n < param.nTrajectory; n++) 
  {
    //s_total
    VectorXd s_total;
    s_total.setZero(3 * nAtom);;//a vector of spins for all the atoms 
    for (int i = 0; i < nAtom; i++) 
    {
      s_total[3 * i] = ensemble.atoms[i].internal.sx(n);
      s_total[3 * i + 1] = ensemble.atoms[i].internal.sy(n);
      s_total[3 * i + 2] = ensemble.atoms[i].internal.sz(n);
    }

    //Stochastic integration
    int num_order = 2;
    stochasticIntegration(s_total, rabiMult, param, num_order);

    //Put back
    for (int i = 0; i < nAtom; i++) 
    {
      ensemble.atoms[i].internal.sx(n) = s_total[3 * i];
      ensemble.atoms[i].internal.sy(n) = s_total[3 * i + 1];
      ensemble.atoms[i].internal.sz(n) = s_total[3 * i + 2];
    }
  }
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param)
{
  //When Doppler effects are considered, should advance internal state first.
  advanceInternalStateOneTimeStep(ensemble, param); //Including both atoms and cavity
  advanceExternalStateOneTimeStep(ensemble, param);
}

void advanceInterval(Ensemble& ensemble, const Param& param)
{
  //Newly added atoms are in the tail of the "atoms" vector, so the first atoms 
  //  in the "atoms" vector will be the first to be removed.
  if (param.pois == 0) {
    advanceAtomsOneTimeStep(ensemble, param);
    addAtomsFromVanillaSource(ensemble, param);
    removeAtomsAtWalls(ensemble, param);
  }
  if (param.pois == 1) {
    advanceAtomsOneTimeStep(ensemble, param);
    addAtomsFromPoissonianSource(ensemble, param);
    removeAtomsAtWalls(ensemble, param);
  }
}

void storeObservables(Observables& observables, const Ensemble& ensemble, const Param& param, int nstore)
{
  //Useful parameters
  const int nTrajectory = param.nTrajectory;
  const int nTimeStep = param.tMax / param.dt;
  const int nBin = param.nBin;
  const int nAtom = ensemble.atoms.size();
  const int mAtom = param.mAtom;
  double gc, gd, g0;
  getGammaValues(gc, gd, g0, param);

  //Useful variable vectors
  VectorXd Jx, Jy, rabiMult;
  rabiMult.setZero(nAtom);
  Jx.setZero(nTrajectory);
  Jy.setZero(nTrajectory);
  getRabimult(rabiMult, ensemble, param);
  for (int i = 0; i < nAtom; i++)
  {
    Jx += ensemble.atoms[i].internal.sx;
    Jy += ensemble.atoms[i].internal.sy;
  }
  //field observables//////////////////////////////////////////////////////////////////////////////////
  
  //intensity
  observables.intensity(nstore) = gc * 0.25 * (Jx.array().square().sum()/nTrajectory + Jy.array().square().sum()/nTrajectory);

  //JxMatrix, JyMatrix, and JzMatrix
  observables.JxMatrix.col(nstore) = Jx;
  observables.JyMatrix.col(nstore) = Jy;
  // observables.JzMatrix.col(s) = ensemble.cavity.Jz.col(nStep);
  
  //Atomic observables//////////////////////////////////////////////////////////////////////////////////

  //nAtom
  observables.nAtom(nstore) = nAtom * mAtom;

  // //inversionAve
  // observables.inversionAve(s) = ensemble.cavity.Jz.col(nStep).sum()/nTrajectory/nAtom/mAtom;

  // //spinFinalMatrices
  // observables.sxFinalMatrix.col(s) = ensemble.spinFinalVector.sxFinalVector;
  // observables.syFinalMatrix.col(s) = ensemble.spinFinalVector.syFinalVector;
  // observables.szFinalMatrix.col(s) = ensemble.spinFinalVector.szFinalVector;

  //For quick runs, COMMENT THESE OUT/////////////////////////////////////////////////////////////////////
  if (param.fast == 0) 
  {
    //internal matrices and bin indices
    MatrixXd SX, SY, SZ; 
    //Note: cannot do Matrix SX, SY, SZ = MatrixXd::Zero(nAtom, nTrajectory); have to be this way
    SX = MatrixXd::Zero(nAtom, nTrajectory);
    SY = MatrixXd::Zero(nAtom, nTrajectory);
    SZ = MatrixXd::Zero(nAtom, nTrajectory);
    VectorXd binIndex = VectorXd::Zero(nBin);
    VectorXd binSum = VectorXd::Zero(nBin);
    double binSize = param.yWall*2/nBin;
    
    //Define rabiMultiplier???
    // double k = 2*M_PI/param.lambda;
    // VectorXd rabiEff = VectorXd::Zero(nAtom);
    
    //assign sx_j*cos(kx_j) to SX, etc.
    for (int i = 0; i < nAtom; i++) 
    {
      // //get rabiEff
      // rabiEff[i] = cos(k*ensemble.atoms[i].external.X[2]);???
      //internal indices
      SX.row(nAtom-i-1) = ensemble.atoms[i].internal.sx;// * rabiEff[i];//???
      //Put atoms in SX in such order that new atoms are in the top rows.
      SY.row(nAtom-i-1) = ensemble.atoms[i].internal.sy;// * rabiEff[i];//???
      SZ.row(nAtom-i-1) = ensemble.atoms[i].internal.sz;
      //bin indices
      int binNumber = (ensemble.atoms[i].external.X[1]+param.yWall)/binSize;
      if (binNumber > nBin-1)
        binNumber = nBin-1;
      binIndex[binNumber] += 1;
    }
    for (int i = 1; i < nBin; i++) 
    {
      binSum[i] = binSum[i-1] + binIndex[i-1];
    }
  
    //sxMatrix, syMatrix, szMatrix
    for (int i = 0; i < nBin; i++) 
    {
      observables.sxMatrix(i, nstore) = SX.middleRows(binSum[i], binIndex[i]).sum()/binIndex[i]/nTrajectory;
      observables.syMatrix(i, nstore) = SY.middleRows(binSum[i], binIndex[i]).sum()/binIndex[i]/nTrajectory;
      observables.szMatrix(i, nstore) = SZ.middleRows(binSum[i], binIndex[i]).sum()/binIndex[i]/nTrajectory;
    }

    // //spinSpinCorAve
    // MatrixXd SX2, SY2, SXSY, SYSX = MatrixXd::Zero(nAtom, nAtom);
    // SX2 = SX*SX.transpose();
    // SY2 = SY*SY.transpose();
    // SXSY = SX*SY.transpose();
    // SYSX = SY*SX.transpose();
    // observables.spinSpinCorAve_re(s) = 
    //   0.25*((SX2.sum()*mAtom-SX2.diagonal().sum())
    //       +(SY2.sum()*mAtom-SY2.diagonal().sum()))/nAtom/(mAtom*nAtom-1)/nTrajectory;
    // observables.spinSpinCorAve_im(s) = 
    //   0.25*((SYSX.sum()*mAtom-SYSX.diagonal().sum())
    //       -(SXSY.sum()*mAtom-SXSY.diagonal().sum()))/nAtom/(mAtom*nAtom-1)/nTrajectory;

    // //spinSpinCor between y = y1 and y = y2
    // for (int i = 0; i < nBin; i++) { //Can be optimized to half diagonal, but testing on symmetry first???
    //   MatrixXd SX_1, SX_2, SY_1, SY_2;
    //   SX_1 = SX.middleRows(binSum[i], binIndex[i]);
    //   SY_1 = SY.middleRows(binSum[i], binIndex[i]);
    //   //diagonal terms
    //   MatrixXd SX_1sq, SY_1sq, SY_1_SX_1, SX_1_SY_1;
    //   SX_1sq = SX_1*SX_1.transpose();
    //   SY_1sq = SY_1*SY_1.transpose();
    //   SY_1_SX_1 = SY_1*SX_1.transpose();
    //   SX_1_SY_1 = SX_1*SY_1.transpose();
    //   observables.spinSpinCor_re(i*nBin+i, s) = 
    //     0.25*((SX_1sq.sum()*mAtom-SX_1sq.diagonal().sum())+(SY_1sq.sum()*mAtom-SY_1sq.diagonal().sum()))
    //       /binIndex[i]/(mAtom*binIndex[i]-1)/nTrajectory;
    //   observables.spinSpinCor_im(i*nBin+i, s) = 
    //     0.25*((SY_1_SX_1.sum()*mAtom-SY_1_SX_1.diagonal().sum())+(SX_1_SY_1.sum()*mAtom-SX_1_SY_1.diagonal().sum()))
    //       /binIndex[i]/(mAtom*binIndex[i]-1)/nTrajectory;
    //   //off-diagonal terms
    //   for (int j = i+1; j < nBin; j++) {
    //     SX_2 = SX.middleRows(binSum[j], binIndex[j]);
    //     SY_2 = SY.middleRows(binSum[j], binIndex[j]);
    //     observables.spinSpinCor_re(i*nBin+j, s) = 
    //       0.25*((SX_1*SX_2.transpose()).sum()+(SY_1*SY_2.transpose()).sum())
    //         /binIndex[i]/binIndex[j]/nTrajectory;
    //     observables.spinSpinCor_im(i*nBin+j, s) = 
    //       0.25*((SY_1*SX_2.transpose()).sum()-(SX_1*SY_2.transpose()).sum())
    //         /binIndex[i]/binIndex[j]/nTrajectory;
    //     //The other half diagonal terms
    //     observables.spinSpinCor_re(j*nBin+i, s) = observables.spinSpinCor_re(i*nBin+j, s);
    //     observables.spinSpinCor_im(j*nBin+i, s) = observables.spinSpinCor_im(i*nBin+j, s);
    //   }
  }
  //For quick runs, COMMENT THESE OUT/////////////////////////////////////////////////////////////////////
}

void evolve(Ensemble& ensemble, const Param& param, Observables& observables)
{
  //evolve
  int nTime = param.tMax / param.dt; 

  //For "nTimeStep" number of data, keep "nStore" of them. 
  for (int ntime = 0, nstore = 0; ntime <= nTime; ntime++) 
  {
    if ((long)(ntime + 1) * param.nStore / (nTime + 1) > nstore) 
    {
      storeObservables(observables, ensemble, param, nstore++);
      std::cout << "Data " << nstore << "/" << param.nStore << " stored." << std::endl << std::endl;
    }
    if (ntime != nTime) 
    {
      advanceInterval(ensemble, param);
    }
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables, const Ensemble& ensemble, const Param& param)
{
  std::cout << "Writing data... (This may take several minutes.)" << std::endl << std::endl;
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  // observableFiles.inversionAve << observables.inversionAve << std::endl;
  observableFiles.JxMatrix << observables.JxMatrix << std::endl;
  observableFiles.JyMatrix << observables.JyMatrix << std::endl;
  //debug
  // observableFiles.JzMatrix << observables.JzMatrix << std::endl;
  // observableFiles.sxFinalMatrix << observables.sxFinalMatrix << std::endl;
  // observableFiles.syFinalMatrix << observables.syFinalMatrix << std::endl;
  // observableFiles.szFinalMatrix << observables.szFinalMatrix << std::endl; 
  //debug
  if (param.fast == 0) {
    // observableFiles.spinSpinCorAve_re << observables.spinSpinCorAve_re << std::endl;
    // observableFiles.spinSpinCorAve_im << observables.spinSpinCorAve_im << std::endl;
    // observableFiles.spinSpinCor_re << observables.spinSpinCor_re << std::endl;
    // observableFiles.spinSpinCor_im << observables.spinSpinCor_im << std::endl;
    observableFiles.sxMatrix << observables.sxMatrix << std::endl;
    observableFiles.syMatrix << observables.syMatrix << std::endl;
    observableFiles.szMatrix << observables.szMatrix << std::endl;   
  } 
}

void mkdir(Param& param) 
{ 
  std::string dirName = "./" + param.controlType + "/";
  std::string mkdir_1 = "mkdir "+ dirName; //make a new catogory to store data
  system(mkdir_1.c_str());
  std::string mkdir_2 = "mkdir "+ dirName + param.name; //make a new directory to store data
  system(mkdir_2.c_str());
  std::string cpInput = "cp input.txt " + dirName + param.name;
  system(cpInput.c_str());  
  std::string moveparam = "mv *.dat " + dirName + param.name;
  system(moveparam.c_str());
}

int main(int argc, char *argv[])
{
  //Count time
  clock_t t1,t2;
  t1=clock();
/////////////////////////////////////////////////////////////////////////////

  //Configuration. Calling functions from "config.hpp".
  CmdLineArgs config;
  getOptions(argc, argv, &config);

  //Set up parameters
  Param param;
  getParam (config.configFile, &param);
	
  //Set up initial conditions
  Ensemble ensemble;
  // generateInitialField(ensemble, param);
  Observables observables(param.nStore, param.nTrajectory, param.nBin);

  //Start simulation
  evolve(ensemble, param, observables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables, ensemble, param);
  
  //Move .dat files into the directory named "name"
  mkdir(param);

/////////////////////////////////////////////////////////////////////////////
  //Count time
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  std::cout << "\nThis program takes " << diff << " seconds." << std::endl << std::endl;

  return 0;
}