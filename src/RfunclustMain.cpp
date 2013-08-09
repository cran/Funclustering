#include <Rcpp.h>
#include <map>
#include "Eigen/Eigen"
#include "Funclustering/src/IModel/IModel.h"
#include "Funclustering/src/IModel/Model.h"
#include "Funclustering/src/IAlgo/IAlgo.h"
#include "Funclustering/src/IAlgo/EMAlgo.h"
#include "Funclustering/src/IStrategy/IStrategy.h"
#include "Funclustering/src/IStrategy/SimpleStrategy.h"
#include "conversion.h"

using namespace Eigen;
using namespace std;
using namespace Rcpp;


RcppExport SEXP RfunclustMain (SEXP inpobj, SEXP outpobj)
{
  BEGIN_RCPP
  //Initialize Rcpp objects
  Rcpp::S4 inputObj(inpobj);
  Rcpp::S4 outputObj(outpobj);

  //convert parameters

  SEXP coefsInit=inputObj.slot("coefs");
  NumericMatrix coefR(coefsInit);
  MatrixXd coefs=convertMatrix<Eigen::MatrixXd,Rcpp::NumericMatrix>(coefR);

  SEXP basis=inputObj.slot("basisProd");
  NumericMatrix basisR(basis);
  MatrixXd basisProd=convertMatrix<Eigen::MatrixXd,Rcpp::NumericMatrix>(basisR);

  int K=as<int>(inputObj.slot("K"));
  double thd=as<double>(inputObj.slot("thd"));
  double epsilon=as<double>(inputObj.slot("epsilon"));
  int nbInit=as<int>(inputObj.slot("nbInit"));
  int nbIterInit=as<int>(inputObj.slot("nbIterInit"));
  int nbIteration=as<int>(inputObj.slot("nbIteration"));
  bool increaseDimension=as<bool>(inputObj.slot("increaseDimension"));
  bool hard=as<bool>(inputObj.slot("hard"));

  SEXP fixedDimension=inputObj.slot("fixedDimension");
  NumericVector fixedDimensionR=(fixedDimension);
  VectorXi fixedDim=convertvector<Eigen::VectorXi,Rcpp::NumericVector>(fixedDimensionR);

  // initialize the C++ classes objects
  IModel* model= new Model(coefs, basisProd, K, thd);
  IAlgo* algo=new EMAlgo(model);
  IModel* best=new Model();
  IStrategy* strategy=new SimpleStrategy(algo,best,epsilon,nbInit);

  //Initialize the model
  strategy->initializeModel(increaseDimension,nbIterInit,hard);

  // run the appropriate algorithm
  	  // if the user gives the dimensions, then we run the version with fixed dimensions
  if (fixedDim.size() > 0){
	  strategy -> run(fixedDim,nbIteration);
  }
  // else we run the version with dimensions varying according to the scree test
  	  // then the parameter increaseDimension specify if or not we constraint the dim to increase
  else
	  strategy -> run(increaseDimension,nbIteration);

  //convert output
  outputObj.slot("loglikelihood")=Rcpp::wrap(model->getLoglik());
  outputObj.slot("loglikTotal")=convertvector<Rcpp::NumericVector,Eigen::VectorXd>(strategy->getLoglikTotal());
  outputObj.slot("aic")=Rcpp::wrap(model->getAic());
  outputObj.slot("bic")=Rcpp::wrap(model->getBic());
  outputObj.slot("icl")=Rcpp::wrap(model->getIcl());
  outputObj.slot("cls")=convertvector<Rcpp::NumericVector,Eigen::VectorXi>(model->getClusters());
  outputObj.slot("proportions")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getParamobj().prop);
  outputObj.slot("dimensions")=convertvector<Rcpp::NumericVector,Eigen::VectorXi>(model->getParamobj().r);
  outputObj.slot("dimTotal")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXi>(strategy->getRtotal());
  outputObj.slot("tik")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getWeights());
  outputObj.slot("V")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getParamobj().V);
  outputObj.slot("empty")=Rcpp::wrap(model->getEmpty());

  // free memory
	delete strategy;
	strategy=0;
	delete best;
	best=0;
	delete algo;
	algo=0;
	delete model;
	model=0;


  return outputObj;
  END_RCPP
}


