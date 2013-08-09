
#ifndef EMALGO_H_
#define EMALGO_H_

#include "IAlgo.h"
#include "../IModel/IModel.h"
#include "../IModel/Model.h"

/**
 * @file EMAlgo.h
 * @class EMAlgo
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Friday, June 21, 2013
 * @brief This inherit from the IAlgo class and made to carry out an EM algorithm.
 */
class EMAlgo: public IAlgo {
public:
	/// Default constructor
	EMAlgo();

	/**
	 * We initialize the EMAlgo object with the model on which it must run.
	 * @param model : a pointer to an IModel object.
	 */
	EMAlgo(IModel* model);
	virtual ~EMAlgo();

	/**
	 * This method is to run the EM algorithm. It's a simple method and calls successively on the member m_model
	 * its methods mStep(dimIncrease) and eStep(). It also calls the method logliklihood() to update the likelihood.
	 * @param dimIncrease a boolean paramater to specify if the dimension must constrain to increase.
	 */
	virtual void run(bool dimIncrease);

	virtual void run(VectorXi dimFixed);

};

#endif /* EMALGO_H_ */
