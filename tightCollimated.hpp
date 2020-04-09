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
  Vector3d X;     //position, (x is irrelevant, y is longgitudinal, z is transverse)
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

//Ensemble of atoms
typedef struct 
{
  std::vector<Atom> atoms;
} Ensemble;

//Simulation parameters
typedef struct Param 
{
  //parameter definitions
  double tMax;
  int nStore; // number of times to store observables
  int nTrajectory; //number of trajectories
  int nBin; //number of bins along the y direction
  //beam parameters
  double density;   //the mean number of atoms per unit time;
  double gc; //gammac

  //other parameters
  std::string controlType; //name of the parent directory for a certain controlle variable
  std::string name; //name of the directory to store results
  int fast; //If we do not want any bin-related calculations in order to save time; 0 = not fast, 1 = fast

  //Constructor; initial values
  Param() : tMax(10), nStore(10), nTrajectory(1), nBin(10), density(1.0), 
            gc(0.1), controlType("test"), name("aProgramHasNoName"), fast(0)
  {}
} Param;

///Observables; n is the nTimeStep
typedef struct Observables 
{
  VectorXi nAtom; 
  VectorXd intensity;
  MatrixXd JxMatrix;
  MatrixXd JyMatrix;
  MatrixXd sxMatrix;
  MatrixXd syMatrix;
  MatrixXd szMatrix;
  // MatrixXd szSqMatrix;
  // VectorXd spinSpinCorAve_re;
  // VectorXd spinSpinCorAve_im;
  // MatrixXd spinSpinCor_re;
  // MatrixXd spinSpinCor_im;

  Observables(const int nStore, const int nTrajectory, const int nBin)//, int density)
  {
    nAtom = VectorXi(nStore); 
    intensity = VectorXd(nStore);
    JxMatrix = MatrixXd(nTrajectory, nStore);
    JyMatrix = MatrixXd(nTrajectory, nStore);
    sxMatrix = MatrixXd(nBin, nStore);
    syMatrix = MatrixXd(nBin, nStore);
    szMatrix = MatrixXd(nBin, nStore);
    // szSqMatrix = MatrixXd(density, nStore);
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
                JxMatrix, 
                JyMatrix,
                // spinSpinCorAve_re, 
                // spinSpinCorAve_im,
                // spinSpinCor_re, 
                // spinSpinCor_im,
                sxMatrix,
                syMatrix,
                szMatrix;
                // szSqMatrix;

  //Constructor              
  ObservableFiles() : nAtom("nAtom.dat"), 
                      intensity("intensity.dat"), 
                      JxMatrix("JxMatrix.dat"),
                      JyMatrix("JyMatrix.dat"),
                      // spinSpinCorAve_re("spinSpinCorAve_re.dat"),
                      // spinSpinCorAve_im("spinSpinCorAve_im.dat"),
                      // spinSpinCor_re("spinSpinCor_re.dat"),
                      // spinSpinCor_im("spinSpinCor_im.dat"),
                      sxMatrix("sxMatrix.dat"),
                      syMatrix("syMatrix.dat"),
                      szMatrix("szMatrix.dat")
                      // szSqMatrix("szSqMatrix.dat")
  {}
  
  //Deconstructor
  ~ObservableFiles() 
  {
    nAtom.close();
    intensity.close();
    JxMatrix.close();
    JyMatrix.close();
    // spinSpinCorAve_re.close();
    // spinSpinCorAve_im.close();
    // spinSpinCor_re.close();
    // spinSpinCor_im.close();
    sxMatrix.close();
    syMatrix.close();
    szMatrix.close();
    // szSqMatrix.close();
  }
  
} ObservableFiles;

#endif