#ifndef SIMPLESTRATEGY_H_
#define SIMPLESTRATEGY_H_

#include "IStrategy.h"

/**
 * Our strategy is to run the EM algorithm after an initialization step in which we take a good starting point
 * in order to increase the chance to converge to the maximum of the logliklihood.
 * In the initialization step we run several small (with a small number of iteration) EM and store the best model
 * (with the grater logliklihood) among these small EM. This best model will be the starting point of the long EM.
 * @file SimpleStrategy.h
 * @class SimpleStrategy
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Thursday, July 11, 2013
 * @brief This inherit from the IStrategy class and made to carry out a our strategy to run a EM algorithm.
 */

class SimpleStrategy: public IStrategy {
public:
	/// Default constuctor
	SimpleStrategy();

	/**
	 * This constructor calls the initialization constructor of the IStrategy to set up the parameters algo, epsilon
	 * and nbInit. bestModel is the specific member of the class SimpleStrategy.
	 * @param algo : IAlgo pointer
	 * @param bestModel : IModel pointer and stores the best model among the iterations of the algorithm.
	 * @param epsilon : convergence criterion.
	 * @param nbInit : number of the small EM in the initialization step.
	 */
	SimpleStrategy(IAlgo* algo, IModel* bestModel, double epsilon, int nbInit);
	virtual ~SimpleStrategy();
	virtual void run(bool dimIncrease, int nbIter);
	virtual void run(VectorXi dimFixed, int nbIter);
	virtual void initializeModel(bool increaseDimInit, int nbIterInit, bool hard);



protected:
	/**
	 * Since the logliklihood is not necessary increasing among the iterations, we store in this member
	 * the state of the best model (with the grater logliklihood).
	 */
	IModel* m_bestModel;

};

#endif /* SIMPLESTRATEGY_H_ */
