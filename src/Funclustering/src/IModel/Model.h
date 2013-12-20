/*
 * Model.h
 *
 *  Created on: 21 juin 2013
 *      Author: soueidat
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "IModel.h"

/**
 * @file Model.h
 * @class Model
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Thursday, July 11, 2013
 * @brief This class inherit from the IModel class and implements its virtual methods.
 */
class Model: public IModel {

public:
	/// Default constructor which calls the default constructor of the IModel class.
	Model();

	/**
	 *The initialization constructor which calls the initialization constructor of the IModel class.
	 * @param coefs : The coefficients matrix.
	 * @param basisProd : The matrix of the inner product between the basis functions.
	 * @param nbClust : The number of clusters.
	 * @param thd : The threshold  to consider in the scree test.
	 */
	Model(MatrixXd coefs,MatrixXd basisProd, int nbClust, double thd);
	/**
	 * The copy constructor which calls the copy constrctor of the IModel class.
	 * @param modelToCopy the model to copy in *this.
	 */
	Model(Model const& modelToCopy);
	virtual ~Model();
  virtual void functionalPca(int nHarm, VectorXd weights);
  virtual void mfpca(int nHarm, VectorXd weights);
	virtual void screeTest(int k);
	virtual void screeTestDimIncrease(int k);
	virtual void mStep(bool dimIncrease);
	virtual void mStep(VectorXi dimFixed);
	virtual void eStep();
	virtual VectorXd normalDistribution(VectorXd X, double mean, double std) const;
	virtual void logliklihood();
	virtual void aic();
	virtual void bic();
	virtual void icl();
	virtual void updateClusters();
	virtual void generateHardWeights();
	virtual void generateSoftWeights();
	virtual void setEmpty(bool empty);


};

#endif /* MODEL_H_ */
