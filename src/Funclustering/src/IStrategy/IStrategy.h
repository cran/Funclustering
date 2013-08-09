#ifndef ISTRATEGY_H_
#define ISTRATEGY_H_

#include "../IModel/IModel.h"
#include "../IModel/Model.h"
#include "../IAlgo/IAlgo.h"
#include "../IAlgo/EMAlgo.h"

/**
 * A strategy is the manner with which we carry out our classification and running our algorithm.
 * For example in a strategy class the members can be the number of iteration of the algorithm or the parameters
 * that we suppose to stores among the iterations.
 * @file IStrategy.h
 * @class IStrategy
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Friday, June 21, 2013
 * @brief This is a basic class for the strategy.
 */
class IStrategy {
public:
	/// Default constructor
	IStrategy();

	/**
	 * This a initialization constructor.
	 * @param algo : a pointer to an IAlgo object specifying which algorithm must be use.
	 * @param epsilon : this is a convergence criterion.
	 * @param nbInit : This is the number of small EM in the initialization.
	 */
	IStrategy(IAlgo* algo, double epsilon, int nbInit);

	/// Destructor
	virtual ~IStrategy();

	/**
	 * Selector to get the member m_loglikTotal.
	 * @return  The logliklihood of the model among all iterations of the algorithm.
	 */
	inline VectorXd getLoglikTotal() const {return m_loglikTotal;};

	/**
	 * Selector to get the member m_Rtotal.
	 * @return The dimensions of the curves retained among all iterations of the algorithm.
	 */
	inline MatrixXi getRtotal() const {return m_Rtotal;};

	/**
	 * This is a first version of the run method.
	 * Here we run the algorithm with a variable dimensions according to the screeTest selection.
	 * In addition in this method we can impose to the dimensions to increase after each iterations.
	 * @param dimIncrease a boolean, if true the dimensions must increase after each iteration.
	 * @param nbIter the maximum number of iteration.
	 */
	virtual void run(bool dimIncrease, int nbIter)=0;

	/**
	 * Here we have another version of the method run.
	 * In contrast of the above one, here we run the algorithm with a fixed dimensions among all iterations.
	 * @param dimFixed vector of integer with length the number of cluster.
	 * 		  Each component contains the number of dimensions of the corresponding cluster.
	 * @param nbIter The maximum number of iterations.
	 */
	virtual void run(VectorXi dimFixed, int nbIter)=0;

	/**
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
	virtual void initializeModel(bool increaseDimInit, int nbIterInit, bool hard)=0;


protected:
	/**
	 * A pointer to an IAlgo object, storing the algorithm to be use.
	 */
	IAlgo* m_algo;

	/**
	 * This member specify the number of small EM in the initialization step.
	 */
	int m_nbInit;

	/**
	 * Convergence criterion, we stop the algorithm if the difference between two consecutive logliklihood
	 * is less than epsilon.
	 */
	double m_epsilon;

	/**
	 * Here we store the logliklihood values among the iterations of the algorithm.
	 */
	VectorXd m_loglikTotal;

	/**
	 * In this member we store the dimensions of the curves among the iterations.
	 * Each column of the matrix m_Rtotal corresponds to the dimensions at a given iterations.
	 */
	MatrixXi m_Rtotal;
};

#endif /* ISTRATEGY_H_ */
