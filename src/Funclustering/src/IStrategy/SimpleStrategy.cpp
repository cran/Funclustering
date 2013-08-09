/*
 * SimpleStrategy.cpp
 *
 *  Created on: 21 juin 2013
 *      Author: soueidat
 */

#include "SimpleStrategy.h"

SimpleStrategy::SimpleStrategy() :IStrategy(){
	m_bestModel=0;
	// TODO Auto-generated constructor stub

}

/*
 * This constructor calls the initialization constructor of the IStrategy to set up the parameters algo, epsilon
 * and nbInit. bestModel is the specific member of the class SimpleStrategy.
 * @param algo : IAlgo pointer
 * @param bestModel : IModel pointer and stores the best model among the iterations of the algorithm.
 * @param epsilon : convergence criterion.
 * @param nbInit : number of the small EM in the initialization step.
 */
SimpleStrategy::SimpleStrategy(IAlgo* algo,IModel* bestModel,double epsilon,int nbInit)
:IStrategy(algo,epsilon,nbInit){
	m_bestModel=bestModel;
}

SimpleStrategy::~SimpleStrategy() {
	// TODO Auto-generated destructor stub
}



/*
 * This is a first version of the run method.
 * Here we run the algorithm with a variable dimensions according to the screeTest selection.
 * In addition in this method we can impose to the dimensions to increase after each iterations.
 * @param dimIncrease a boolean, if true the dimensions must increase after each iteration.
 * @param nbIter the maximum number of iteration.
 */
void SimpleStrategy::run(bool dimIncrease, int nbIter){
	// iter will count the number of iterations of the algorithm.
	int iter=1;

	//as its name indicates, logliklihood_prev will store the previous liklihood; ie logliklihood befor runnig algo
	double logliklihood_prev=m_algo->getModel()->getLoglik();
	//logliklihood_max stores the maximum liklihood among the iterations of the algorithm.
	double logliklihood_max=m_algo->getModel()->getLoglik();

	// m_logliktotal and m_Rtotal stores respectively loglik and dimensions among the iteration of the algorithm
	m_loglikTotal=VectorXd(nbIter+1);
	m_loglikTotal(0)=m_algo->getModel()->getLoglik();
	m_Rtotal=MatrixXi(m_algo->getModel()->getNbClust(),nbIter+1);
	m_Rtotal.col(0)=m_algo->getModel()->getParamobj().r;


	//criterion stop: small increase in consecutive loglikelihood or max iterations
	do {
		// storing the loglik befor running algo.
		logliklihood_prev=m_algo->getModel()->getLoglik();

		/*
		 Here is the first version of runnig run algo.
		 In fact we use the version of running algo that permit to specify if the dimensions must increase
		 or let them change with the scree test
		 */
		m_algo -> run(dimIncrease);
/*
		//update clusters and chck if there is an empty class. here is the only place that the clusters membership is calculate
		m_algo->getModel()->updateClusters();
*/
		//stores the new loglikelihood and the news dimensions
		m_loglikTotal(iter)=m_algo->getModel()->getLoglik();
		m_Rtotal.col(iter)=m_algo->getModel()->getParamobj().r;

		//if we have increased the loglikelihood, this models is copied in bestModels
		if (m_algo->getModel()->getLoglik() > logliklihood_max){
			//update the logliklihood_max.
			logliklihood_max=m_algo->getModel()->getLoglik();

			// we stores the current model which make increasing the logliklihood in the bestModel.
			*m_bestModel=*(m_algo->getModel());
		}
/*
		// this to fix the empty class error
		if (m_algo->getModel()->getLoglik() != m_algo->getModel()->getLoglik()) {
			break;
		}
*/
		iter+=1;
	} while (fabs(m_algo->getModel()->getLoglik()-logliklihood_prev) > m_epsilon &&  iter <= nbIter
			&& m_algo->getModel()->getEmpty()==false);

	// we stores only the liklihood and dimension which have been calculated.
	m_loglikTotal.conservativeResize(iter);
	m_Rtotal.conservativeResize(m_algo->getModel()->getNbClust(),iter);

	//we update the actual model with the best model
	*(m_algo->getModel())=*m_bestModel;

	//compute the new criterion
	m_algo->getModel()->aic();
	m_algo->getModel()->bic();
	m_algo->getModel()->icl();


}



/*
 * Here we have another version of the method run.
 * In contrast of the above one, here we run the algorithm with a fixed dimensions among all iterations.
 * @param dimFixed vector of integer with length the number of cluster.
 * 		  Each component contains the number of dimensions of the corresponding cluster.
 * @param nbIter The maximum number of iterations.
 */
void SimpleStrategy::run(VectorXi dimFixed, int nbIter){
	// iter will count the number of iterations of the algorithm.
	int iter=1;

	//as its name indicates, logliklihood_prev will store the previous liklihood; ie logliklihood befor runnig alg
	double logliklihood_prev=m_algo->getModel()->getLoglik();
	//logliklihood_max stores the maximum liklihood among the iterations of the algorithm.
	double logliklihood_max=m_algo->getModel()->getLoglik();
	// m_logliktotal and m_Rtotal stores respectively loglik and dimensions among the iteration of the algorithm
	m_loglikTotal=VectorXd(nbIter+1);
	m_loglikTotal(0)=m_algo->getModel()->getLoglik();
	m_Rtotal=MatrixXi(m_algo->getModel()->getNbClust(),nbIter+1);
	m_Rtotal.col(0)=m_algo->getModel()->getParamobj().r;

	//criterion stop: small increase in consecutive loglikelihood or max iterations
	do {
		// storing the loglik befor running algo.
		logliklihood_prev=m_algo->getModel()->getLoglik();

		/*
		 run algo. here we use the version of running algo with a fixed dimensions among all iterations.
		 */
		m_algo -> run(dimFixed);

		//update clusters and chck if there is an empty class. here is the only place that the clusters membership is calculate
		m_algo->getModel()->updateClusters();

		//stores the new loglikelihood and the news dimensions
		m_loglikTotal(iter)=m_algo->getModel()->getLoglik();
		m_Rtotal.col(iter)=m_algo->getModel()->getParamobj().r;

		//if we have increased the loglikelihood, this models is copied in bestModels
		if (m_algo->getModel()->getLoglik() > logliklihood_max){
			//update the logliklihood_max.
			logliklihood_max=m_algo->getModel()->getLoglik();

			//update clusters. here is the only place that the clusters membership is calculate
			m_algo->getModel()->updateClusters();

			// we stores the current model which make increasing the logliklihood in the bestModel.
			*m_bestModel=*(m_algo->getModel());
		}

		iter+=1;
	} while (fabs(m_algo->getModel()->getLoglik()-logliklihood_prev) > m_epsilon &&  iter <= nbIter
			&& m_algo->getModel()->getEmpty()==false);

	// we stores only the liklihood and dimension which have been calculated.
	m_loglikTotal.conservativeResize(iter);
	m_Rtotal.conservativeResize(m_algo->getModel()->getNbClust(),iter);

	//we update the actual model with the best model
	*(m_algo->getModel())=*m_bestModel;

	//compute the new criterion
	m_algo->getModel()->aic();
	m_algo->getModel()->bic();
	m_algo->getModel()->icl();
}

/*
 * Since the EM algorithm made the logliklihood converge to a local maximum, we make an initialization
 * step in which we run m_nbInit small EM, ie with small number of iterations.
 * From each EM we take a bestModel (with a greater logliklihood), and finally we take the best Model
 * around the all algorithms that have been ran.
 * At the end of the initialization step we assign to the members of the bestModel, their default values
 * except the weights matrix which will be the starting point for the long algorithm that will be run
 * after initialization.
 * The logliklihood will set to -Infinity, the dimensions to zeros. And the other are reinitialize in the mStep.
 * @param increaseDimInit
 * @param nbIterInit the number of iteration in the initialization part of the algorithm,
 * 			this number must be less than the number of iteration of the long EM.
 * @param hard to specify which kind of initialization we use.
 * 		  if true we will generate a weights matrix according to generateHardWeights,
 * 	      if false we use generateSoftWeights.
 */
void SimpleStrategy::initializeModel(bool increaseDimInit, int nbIterInit, bool hard){

	/* weightsInit will stores the weights matrix of the model with the maximum liklihood
	among the m_nbInit initialization.
	*/
	MatrixXd weightsInit;
	for (int i = 0; i < m_nbInit; ++i) {
		// l_max is suppose to store the maximum liklihood among the m_nbInit initialization.
		double l_max=m_algo->getModel()->getLoglik();
		//chose one of two possibilities: hardweights or softweigths (cf. the documentation of these two methods)
		if (hard==true){
			m_algo->getModel()->generateHardWeights();
		}
		else {m_algo->getModel()->generateSoftWeights();}

		//running strategy with a small number of iteration.
		run(increaseDimInit,nbIterInit);

		/*
		 Since the logliklihood is not an increasing function, we have to store the weights matrix from
		 the bestModel (with the grater logliklihood).
		 */
		if (m_algo->getModel()->getLoglik() > l_max) {
			//update l_max.
			l_max=m_algo->getModel()->getLoglik();

			//update weightsInit
			weightsInit=m_algo->getModel()->getWeights();
		}
	}
	// set weights with the weights of the best model among the all previous initializations.
	m_algo->getModel()->setWeights(weightsInit);

	// set logliklihood to it's default value (-Infinity) in order to stars the big algorithm with a standard model.
	m_algo->getModel()->setLoglikToInfinity();

	// set dimensions to it's default value (zero) in order to stars the big algorithm with a standard model.
	m_algo->getModel()->setDimensionsToZero();
}
