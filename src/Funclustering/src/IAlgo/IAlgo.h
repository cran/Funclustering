#ifndef IALGO_H_
#define IALGO_H_

#include "../IModel/IModel.h"
#include "../IModel/Model.h"

/**
 * @file IAlgo.h
 * @class IAlgo
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Friday, June 21, 2013
 * @brief This is a basic class for the algorithms.
 */
class IAlgo {
public:

	///Default constructor
	IAlgo();

	/**
	 * This constructor initialize an algorithm with a model. In fact the algorithm will run with the parameters
	 * of the model.
	 * @param model a pointer to an object of type IModel. *model is the model on which the algorithm must run.
	 */
	IAlgo(IModel* model);

	/**
	 * Selector to get the model.
	 * @return a pointer to the model.
	 */
	inline IModel* getModel() const {return m_model;};

	/// Destructor
	virtual ~IAlgo();

	/**
	 * This is the main function of this class. Here we run the algorithm.
	 * @param dimIncrease a parameter to specify if the dimensions must be constrain to increase.
	 */
	virtual void run(bool dimIncrease)=0;

	/**
	 * The only difference between this version of running algorithm and the previous one, is that in this version
	 * we run the mStep with fixed dimensions.
	 * @param dimFixed
	 */
	virtual void run(VectorXi dimFixed)=0;

protected:
	/**
	 * A pointer to the model on which we have to running the algorithm.
	 */
	IModel* m_model;

};

#endif /* IALGO_H_ */
