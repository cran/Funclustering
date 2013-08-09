/*
 * mfpca.cpp
 *
 *  Created on: 29 juil. 2013
 *      Author: soueidat
 */

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


RcppExport SEXP mfpcaCpp (SEXP inpobj, SEXP nharmR, SEXP tikR/*, SEXP mpcaobj*/)
{
  BEGIN_RCPP
  //Initialize Rcpp objects
  Rcpp::S4 inputObj(inpobj);
  //Rcpp::S4 pcaobj(mpcaobj);


  //convert parameters (input)
  int nharm=as<int>(nharmR);

  NumericVector tikVect(tikR);
  VectorXd tik=convertvector<Eigen::VectorXd,Rcpp::NumericVector>(tikVect);

  SEXP coefsInit=inputObj.slot("coefs");
  NumericMatrix coefR(coefsInit);
  MatrixXd coefs=convertMatrix<Eigen::MatrixXd,Rcpp::NumericMatrix>(coefR);

  SEXP basis=inputObj.slot("basisProd");
  NumericMatrix basisR(basis);
  MatrixXd basisProd=convertMatrix<Eigen::MatrixXd,Rcpp::NumericMatrix>(basisR);

  int K=as<int>(inputObj.slot("K"));
  double thd=as<double>(inputObj.slot("thd"));


  // initialize the C++ classes objects
  IModel* model= new Model(coefs, basisProd, K, thd);


  // running functional pca
  model->mfpca(nharm,tik);

  // converting output
 // pcaobj.slot("values")= convertvector<Rcpp::NumericVector,Eigen::VectorXd>(model->getPca().values);
  //pcaobj.slot("varprop")= convertvector<Rcpp::NumericVector,Eigen::VectorXd>(model->getPca().varprop);
  //pcaobj.slot("scores")= convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getPca().scores);

  Rcpp::List pca=List::create(
		  Named("harmonics")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getPca().harmonics),
		  Named("values")=convertvector<Rcpp::NumericVector,Eigen::VectorXd>(model->getPca().values),
		  Named("varprop")=convertvector<Rcpp::NumericVector,Eigen::VectorXd>(model->getPca().varprop),
		  Named("scores")=convertMatrix<Rcpp::NumericMatrix,Eigen::MatrixXd>(model->getPca().scores)
		  	  	  	  	  	  	  );

  delete model;
  model=0;

  //return pcaobj;
  return pca;
  END_RCPP
}



