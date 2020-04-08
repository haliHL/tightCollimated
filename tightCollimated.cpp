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

  while (!configInput.eof()) 
  {
    configInput >> dummy;
    if (configInput.eof()) break;
    if (!configInput.good()) 
    {
      std::cout << "Bad read in input file" << std::endl;
      exit(-1);
    }
    else if (dummy.compare("tMax") == 0)
      configInput >> param->tMax;
    else if (dummy.compare("nStore") == 0)
      configInput >> param->nStore;
    else if (dummy.compare("nTrajectory") == 0)
      configInput >> param->nTrajectory;
    else if (dummy.compare("nBin") == 0)
      configInput >> param->nBin;
    else if (dummy.compare("density") == 0)
      configInput >> param->density;
    else if (dummy.compare("gc") == 0)
      configInput >> param->gc;
    else if (dummy.compare("controlType") == 0)
      configInput >> param->controlType;    
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else if (dummy.compare("fast") == 0)
      configInput >> param->fast;
    else 
    {
      std::cout << "Error: invalid label " << dummy << " in " << filename << std::endl;
      exit(-1);
    }
  }
}

void generateExternalState(Atom& newAtom, const Param& param)
{
  Vector3d X (0, 0, 0);
  Vector3d P (0, 1, 0); //2*ywall = 1, tau = 1
  //Complete initiation
  External newExternal = {X, P};
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
  
  //Random initialization for sx and sy.
  //+1, -1 approach//////////////////////////////////////////////////
  // for (int j = 0; j < nTrajectory; j++)
  // {
  //   newSx(j) += double(rng.get_binomial_int(0.5, 1)) * 2 - 1;//50percent giving 1 or -1
  //   newSy(j) += double(rng.get_binomial_int(0.5, 1)) * 2 - 1;//50percent giving 1 or -1
  // } 
  //random phase approach//////////////////////////////////////////////////
  for (int j = 0; j < nTrajectory; j++)
  {
    double phi = rng.get_uniform_rn(0, 2 * M_PI);
    newSx(j) = sqrt(2) * cos(phi); 
    newSy(j) = sqrt(2) * sin(phi); 
  }
  // std::cout << "Using random phase initialization." << std::endl;


  //Complete initiation
  Internal newInternal = {newSx, newSy, newSz};
  newAtom.internal = newInternal;
}

void addAtomsFromVanillaSource(Ensemble& ensemble, const Param& param)
{
  Atom newAtom; //Create a new atom
  generateExternalState(newAtom, param);    
  generateInternalState(newAtom, param);           
  ensemble.atoms.push_back(newAtom);
}

void removeAtomsAtWalls(Ensemble& ensemble, const Param& param) 
{
  std::vector<Atom> newAtoms;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++) 
  {
    if (a->external.X(1) < 1) 
    {
      newAtoms.push_back(*a);
    }
      
  } 
  ensemble.atoms = newAtoms;
}

void advanceExternalStateOneTimeStep(Ensemble& ensemble, const Param& param) 
{
  double dt = 1.0 / param.density;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    a->external.X += dt * a->external.P;
}

void getDriftVector(VectorXd& drift, const VectorXd& s_total, const Param& param) 
{
  //For convenience
  const int nAtom = s_total.size() / 3;
  const double gc = param.gc;

  //Define Jx and Jy
  double Jx = 0, Jy = 0;
  for (int i = 0; i < nAtom; i++)
  {
    Jx += s_total(3 * i);
    Jy += s_total(3 * i + 1);
  }

  //Drift vector terms. 
   for (int i = 0; i < nAtom; i++) 
   {
      double sx_i = s_total(3 * i);
      double sy_i = s_total(3 * i + 1);
      double sz_i = s_total(3 * i + 2);
      drift(3 * i) = gc / 2 * (Jx * sz_i - sx_i * (sz_i + 1));      
      drift(3 * i + 1) = gc / 2 * (Jy * sz_i - sy_i * (sz_i + 1));
      drift(3 * i + 2) = - gc * (sz_i + 1) - gc / 2 * (Jx * sx_i + Jy * sy_i - (pow(sx_i, 2) + pow(sy_i, 2)));
   }
}

void getDiffusionVector(VectorXd& dW, const VectorXd& s_total, const Param& param)
{
  //For convenience
  const int nAtom = s_total.size() / 3;
  const double gc = param.gc;
  const double dt = 1.0 / param.density;

  //Define the 2 wiener processes
  double dW_q, dW_p;
  dW_q= rng.get_gaussian_rn(sqrt(dt));
  dW_p= rng.get_gaussian_rn(sqrt(dt));

  // diffusion for the atoms
  for (int i = 0; i < nAtom; i++) 
  {
    double sx_i = s_total(3 * i);
    double sy_i = s_total(3 * i + 1);
    double sz_i = s_total(3 * i + 2);
    dW(3 * i) = - sqrt(gc) * sz_i * dW_p;
    dW(3 * i + 1) = sqrt(gc) * sz_i * dW_q;
    dW(3 * i + 2) = sqrt(gc) * (sx_i * dW_p - sy_i * dW_q);
  }
}

void stochasticIntegration(VectorXd& s_total, const Param& param)
{
  //Useful parameters
  const int nAtom = s_total.size() / 3;
  const double dt = 1.0 / param.density;

  //Drift
  VectorXd drift;
  drift.setZero(3 * nAtom);
  getDriftVector(drift, s_total, param);

  //Diffusion
  VectorXd dW;
  dW.setZero(3 * nAtom);
  getDiffusionVector(dW, s_total, param);

  //1st order
  s_total += drift * dt + dW;
}

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param)
{
  //Useful parameters
  const int nAtom = ensemble.atoms.size();

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
    stochasticIntegration(s_total, param);

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
  advanceAtomsOneTimeStep(ensemble, param);
  addAtomsFromVanillaSource(ensemble, param);
  removeAtomsAtWalls(ensemble, param);
}

void storeObservables(Observables& observables, const Ensemble& ensemble, const Param& param, int nstore)
{
  //Useful parameters
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  const double gc = param.gc;
  const int nBin = param.nBin;

  //Useful variable vectors
  VectorXd Jx, Jy;
  Jx.setZero(nTrajectory);
  Jy.setZero(nTrajectory);
  for (int i = 0; i < nAtom; i++)
  {
    Jx += ensemble.atoms[i].internal.sx;
    Jy += ensemble.atoms[i].internal.sy;
  }
  //field observables//////////////////////////////////////////////////////////////////////////////////
  
  //intensity
  observables.intensity(nstore) = gc * 0.25 * (Jx.array().square().sum()/nTrajectory + Jy.array().square().sum()/nTrajectory);

  //JxMatrix, JyMatrix
  observables.JxMatrix.col(nstore) = Jx;
  observables.JyMatrix.col(nstore) = Jy;
  
  //Atomic observables//////////////////////////////////////////////////////////////////////////////////

  //nAtom
  observables.nAtom(nstore) = nAtom;

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
    double binSize = 1.0 / nBin;
    for (int i = 0; i < nAtom; i++) 
    {
      //internal indices
      SX.row(nAtom-i-1) = ensemble.atoms[i].internal.sx;//
      //Put atoms in SX in such order that new atoms are in the top rows.
      SY.row(nAtom-i-1) = ensemble.atoms[i].internal.sy;//
      SZ.row(nAtom-i-1) = ensemble.atoms[i].internal.sz;
      //bin indices
      int binNumber = (ensemble.atoms[i].external.X(1)) / binSize;
      if (binNumber > nBin - 1)
      {
        binNumber = nBin - 1;
      }
      binIndex(binNumber) += 1;
    }
    for (int i = 1; i < nBin; i++) 
    {
      binSum(i) = binSum(i - 1) + binIndex(i - 1);
    }
  
    //sxMatrix, syMatrix, szMatrix
    for (int i = 0; i < nBin; i++) 
    {
      observables.sxMatrix(i, nstore) = SX.middleRows(binSum(i), binIndex(i)).sum() / binIndex(i) / nTrajectory;
      observables.syMatrix(i, nstore) = SY.middleRows(binSum(i), binIndex(i)).sum() / binIndex(i) / nTrajectory;
      observables.szMatrix(i, nstore) = SZ.middleRows(binSum(i), binIndex(i)).sum() / binIndex(i) / nTrajectory;
      observables.szSqMatrix(i, nstore) = SZ.middleRows(binSum(i), binIndex(i)).squaredNorm() / binIndex(i) / nTrajectory;
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
  int nTime = param.tMax * param.density; //dt = 1 / density

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

void writeObservables(ObservableFiles& observableFiles, Observables& observables, const Ensemble& ensemble, const Param& param)
{
  std::cout << "Writing data... (This may take several minutes.)" << std::endl << std::endl;
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.JxMatrix << observables.JxMatrix << std::endl;
  observableFiles.JyMatrix << observables.JyMatrix << std::endl;
  if (param.fast == 0) {
    // observableFiles.spinSpinCorAve_re << observables.spinSpinCorAve_re << std::endl;
    // observableFiles.spinSpinCorAve_im << observables.spinSpinCorAve_im << std::endl;
    // observableFiles.spinSpinCor_re << observables.spinSpinCor_re << std::endl;
    // observableFiles.spinSpinCor_im << observables.spinSpinCor_im << std::endl;
    observableFiles.sxMatrix << observables.sxMatrix << std::endl;
    observableFiles.syMatrix << observables.syMatrix << std::endl;
    observableFiles.szMatrix << observables.szMatrix << std::endl;
    observableFiles.szSqMatrix << observables.szSqMatrix << std::endl;      
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