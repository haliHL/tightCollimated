#ifndef __TIGHTCOLLIMATED__HPP__
#define __TIGHTCOLLIMATED__HPP__
//This program is used to simulate the tight collimated beam laser
//  --using the cNumber method 
//  --using the Adiabatic Elimination Approximation

//Include Eigen package
//Work in Eigen namespace
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
using namespace Eigen;

//Include standard packages
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <stdlib.h> 
//Define pi
#define _USE_MATH_DEFINES
#include <cmath>

//Include and define RNG
#include "RNG.hpp"
RNG rng(time(NULL));

//Define data structure

//Atom external states
typedef struct 
{
  Vector3d X;     //position
  Vector3d P;     //momentum. We suppose mass is one, so momentum is velocity.
} External;

//Atom internal states
typedef struct 
{
  VectorXd sx;       //sigma_x. Dim: nTrajectory
  VectorXd sy;       //sigma_y. Dim: nTrajectory
  VectorXd sz;       //sigma_z. Dim: nTrajectory
} Internal;

//Atom total states
typedef struct 
{
  External external;  //The position and velocity of an atom.
  Internal internal;  //The internal states of an atom. We keep track of sx, sy, and sz
                        //of a single atom at all times for all trajectories.
} Atom;

// //Cavity field enslaved to collective variables
// typedef struct {
//   MatrixXd Jx;        // \sum \cos{kz_j}s_j^x Dim: nTrajectory*(nTimeStep+1)
//   MatrixXd Jy;        // \sum \cos{kz_j}s_j^y Dim: nTrajectory*(nTimeStep+1)
//   MatrixXd Jz;        // \sum s_j^z Dim: nTrajectory*(nTimeStep+1)
// } Cavity;

// //Final spin matrices
// typedef struct 
// {
//   VectorXd sxFinalVector;        //Dim: nTrajectory
//   VectorXd syFinalVector;        //Dim: nTrajectory
//   VectorXd szFinalVector;        //Dim: nTrajectory
// } SpinFinalVector;

//Ensemble of atoms
typedef struct 
{
  std::vector<Atom> atoms;
  // Cavity cavity;
  // SpinFinalVector spinFinalVector; 
} Ensemble;

//Simulation parameters
typedef struct Param 
{
  //parameter definitions
  double dt; 
  double tMax;
  int nStore; // number of times to store observables
  int nTrajectory; //number of trajectories
  int nBin; //number of bins along the y direction
  //beam parameters
  double yWall; //position of the wall where atoms are destroyed
                //The coordinated are chosen s.t. atoms are created at -yWall.
                //The walls are assumed to be in xz plane.
  double lambda;  //The wavelength of the laser light
  double deltaZ;  //z direction standard deviation. 
                    //y direction is taken care of by the Poisson distribution. 
                    //x direction is irrelevant, but we still keep it here.
  double deltaPz;  //standard deviation of z velocity
  double tau;     //the transit time tau1>0, unit 1/gammaC.
  double density;   //the mean number of atoms per unit time;
  int mAtom; //the number of atoms for each collective variable
  double rabi; //Single-photon rabi
  double kappa;  //Cavity decay rate
  double detuning; //Cavity detuning from the atomic frequency, defined by omegac-omegaa, i.e., cavity frequency - atomic frequency
  double T2; //T2 dephasing.

  //other parameters
  std::string controlType; //name of the parent directory for a certain controlle variable
  std::string name; //name of the directory to store results
  int pois; //If the Poissonian noise is included; 0 = no, 1 = yes
  int fast; //If we do not want any bin-related calculations in order to save time; 0 = not fast, 1 = fast

  //Constructor; initial values
  Param() : dt(0.01), tMax(10), 
            nStore(10), nTrajectory(1), nBin(10), 
            yWall(10.0), lambda(1.0),
            deltaZ(0.0), deltaPz(0.0), 
            tau(1.0), density(1.0), mAtom(1),
            rabi(10), kappa(1000), detuning(0), T2(0.0),
            controlType("test"), name("aProgramHasNoName"), pois(0), fast(0)
  {}

} Param;

///Observables; n is the nTimeStep
typedef struct Observables 
{
  VectorXi nAtom; 
  VectorXd intensity;
  // VectorXd inversionAve;
  MatrixXd JxMatrix;
  MatrixXd JyMatrix;
  // MatrixXd JzMatrix;
  // MatrixXd sxFinalMatrix;
  // MatrixXd syFinalMatrix;
  // MatrixXd szFinalMatrix;
  MatrixXd sxMatrix;
  MatrixXd syMatrix;
  MatrixXd szMatrix;
  // VectorXd spinSpinCorAve_re;
  // VectorXd spinSpinCorAve_im;
  // MatrixXd spinSpinCor_re;
  // MatrixXd spinSpinCor_im;

  Observables(const int nStore, const int nTrajectory, const int nBin)
  {
    nAtom = VectorXi(nStore); 
    intensity = VectorXd(nStore);
    // inversionAve = VectorXd(nStore);
    JxMatrix = MatrixXd(nTrajectory, nStore);
    JyMatrix = MatrixXd(nTrajectory, nStore);
    // JzMatrix = MatrixXd(nTrajectory, nStore);
    // sxFinalMatrix = MatrixXd(nTrajectory, nStore);
    // syFinalMatrix = MatrixXd(nTrajectory, nStore);
    // szFinalMatrix = MatrixXd(nTrajectory, nStore);
    sxMatrix = MatrixXd(nBin, nStore);
    syMatrix = MatrixXd(nBin, nStore);
    szMatrix = MatrixXd(nBin, nStore);
    // spinSpinCorAve_re = VectorXd(nStore);
    // spinSpinCorAve_im = VectorXd(nStore);
    // int nBinSquare = nBin * nBin;
    // spinSpinCor_re = MatrixXd(nBinSquare, nStore);
    // spinSpinCor_im = MatrixXd(nBinSquare, nStore);
  } 

} Observables;

typedef struct ObservableFiles 
{
  //Definition
  std::ofstream nAtom, 
                intensity, 
                // inversionAve, 
                JxMatrix, 
                JyMatrix,
                // JzMatrix,
                // sxFinalMatrix,
                // syFinalMatrix,
                // szFinalMatrix,
                // spinSpinCorAve_re, 
                // spinSpinCorAve_im,
                // spinSpinCor_re, 
                // spinSpinCor_im,
                sxMatrix,
                syMatrix,
                szMatrix;

  //Constructor              
  ObservableFiles() : nAtom("nAtom.dat"), 
                      intensity("intensity.dat"), 
                      // inversionAve("inversionAve.dat"),
                      JxMatrix("JxMatrix.dat"),
                      JyMatrix("JyMatrix.dat"),
                      // JzMatrix("JzMatrix.dat"),                      
                      // sxFinalMatrix("sxFinalMatrix.dat"),
                      // syFinalMatrix("syFinalMatrix.dat"),
                      // szFinalMatrix("szFinalMatrix.dat"),
                      // spinSpinCorAve_re("spinSpinCorAve_re.dat"),
                      // spinSpinCorAve_im("spinSpinCorAve_im.dat"),
                      // spinSpinCor_re("spinSpinCor_re.dat"),
                      // spinSpinCor_im("spinSpinCor_im.dat"),
                      sxMatrix("sxMatrix.dat"),
                      syMatrix("syMatrix.dat"),
                      szMatrix("szMatrix.dat")
  {}
  
  //Deconstructor
  ~ObservableFiles() 
  {
    nAtom.close();
    intensity.close();
    // inversionAve.close();
    JxMatrix.close();
    JyMatrix.close();
    // JzMatrix.close();
    // sxFinalMatrix.close();
    // syFinalMatrix.close();
    // szFinalMatrix.close();
    // spinSpinCorAve_re.close();
    // spinSpinCorAve_im.close();
    // spinSpinCor_re.close();
    // spinSpinCor_im.close();
    sxMatrix.close();
    syMatrix.close();
    szMatrix.close();
  }
  
} ObservableFiles;

#endif