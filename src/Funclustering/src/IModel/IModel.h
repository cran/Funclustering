#ifndef IMODEL_H_
#define IMODEL_H_

#include "../../../Eigen/Eigen"
#include <math.h>
#include <limits>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
using namespace Eigen;
using namespace std;

/**
 * This structure named param is created to store the parameters of the model updated in the mStep().
 */
struct param{
	/**
	 * prop correspond to the proportion of each cluster in the classification.
	 */
	MatrixXd prop;

	/**
	 * m is the mean matrix of the curves in each cluster.
	 */
	MatrixXd m;

	/**
	 * V is the covariance operator matrix of the curves in each cluster.
	 */
	MatrixXd V;

	/**
	 * scores is an c++ array that contains the matrixs of scores in each cluster.
	 */
	MatrixXd* scores;

	/**
	 * r is the vector of the number of the principal component of the pca inside each cluster.
	 */
	VectorXi r;
};

/**
 * This structure named funPca stores the results of the functional pca.
 */
struct funPca{
	/**
	 * scores is the scores matrix.
	 */
	MatrixXd scores;

	/**
	 * values is the eigen values.
	 */
	VectorXd values;

	/**
	 * varprop corresponds to the eigenvalues divided by their sum.
	 */
	VectorXd varprop;

	/**
	 * The harmonics are the coefficients of the eigen function in the basis expansion
	 */
	MatrixXd harmonics;
};

/**
 * This is a basic class for models. In this class we construct the model through its principals attributes. \n
 * These attributes are the members of the class. By the way, a model is the functional data object with all
 * methods and parameters permitting it's clustering by the methods made by J.Jacques & all using the EM algorithm.
 * @file IModel.h
 * @class IModel
 * @author Mohamed Soueidatt
 * @version 1.0
 * @date Monday, May 13, 2013
 * @brief This class is created to manage the parameters of the model for clustering functional data.
 */
class IModel {

protected:
	/**
	 * In the univariate case this is the transposed matrix of the coefficients matrix \n
	 * in the basis expansion of the functional data. \n
	 * In the multivariate case this member will be the concatenation of the transposed matrix coefficients
	 * of each dimension of the data.
	 */
	MatrixXd m_coefs;

	/**
	 * In the univariate case this is the matrix of the inner product between the basis functions \n
	 * In the multivariate case the m_basisProd member will be a block diagonal matrix and each block \n
	 * will be the matrix of the inner product between the basis functions of each dimension of the data.
	 */
	MatrixXd m_basisProd;

	/**
	 * The number of cluster to be achieve
	 */
	int m_nbClust;

	/**
	 * This member stores the parameters updated in the mStep() and is described above in the description \n
	 * of the structure param.
	 */
	param m_paramobj;

	/**
	 * In this member we stores the results of the functional pca and is described above in the description \n
	 * of the structure funPca.
	 */
	funPca m_pca;

	/**
	 * Matrix of curves weights such as in each column we have the weights of the curves \n
	 * in the corresponding cluster
	 */
	MatrixXd m_weights;

	/**
	 * The threshold in the scree test to determine the number principal component \n
	 * to be retain in the functional pca
	 */
	double m_thd;

	/**
	 * The logliklihood of the model.
	 */
	double m_loglik;

	/**
	 * The aic criterion of the model.
	 */
	double m_aic;

	/**
	 * The bic criterion of the model.
	 */
	double m_bic;

	/**
	 * The icl criterion of the model.
	 */
	double m_icl;
	/**
	 * This vector assigns to each curve the index of the corresponding cluster.
	 */
	VectorXi m_clusters;

	/**
	 * The coefficient Z(i,g) of the matrix m_Z is equal to 1 if the i'th curve belong \n
	 * in the cluster g and it's equal to zero otherwise.
	 */
	MatrixXi m_Z;

	/**
	 * check empty class
	 */
	bool m_empty;


public:
	/// Default constructor, here we construct a default model.
	IModel();

	/**
	 * This constructor construct a model by initializing it's members by the parameters given in input.
	 * In fact this two first input parameters will be provided by R using the fda package.
	 * @param coefs : the matirx which sets the m_coefs member.
	 * @param basisProd : sets the matrix of the basis function inner product.
	 * @param nbClust : the number of clusters.
	 * @param thd : the threshold for the scree test.
	 */
	IModel(MatrixXd coefs, MatrixXd basisProd, int nbClust, double thd);

	/**
	 * Selector to set a weights matrix to the matrix given in the input.
	 * @param weights : a matrix which will set m_weights.
	 */
	void setWeights(MatrixXd const& weights);

	/**
	 * This method set the liklihood to -inf. we call it after the initialization to restart the EM algorithm
	 * with the smallest liklihood possible and thus we make sure that the liklihood will increase.
	 */
	void setLoglikToInfinity();

	/**
	 * This method set the dimensions to zero in each cluster. we call it after the initialization
	 * to restart the EM algorithm with the basic parameters.
	 */
	void setDimensionsToZero();

	/**
	 * This is copy constructor and it lead in particular to make a deep copy of the object
	 * when it's needed. we assign to each member of *this it's value in modelToCopy.
	 * @param modelToCopy an model to be copied in *this.
	 */
	IModel(IModel const& modelToCopy);

	/**
	 * This is the overload of the assignment operator and it help to make a deep copy when it's needed.
	 * @param modelToCopy : model to be copied in *this.
	 * @return *this.
	 */
	IModel& operator= (IModel const& modelToCopy);

	/**
	 *  Here we destruct the IModel object by making a delete statement on the only pointer in members
	 *  of the IModel object, ie: m_paramobj.scores. The other members are eigen object and the eigen destructor
	 *  will be called automatically.
	 */

	virtual ~IModel();

	/**
	 * Selector to get the coefficient matrix.
	 * @return m_coefs.
	 */
	inline MatrixXd getCoefs() const {return m_coefs;};

	/**
	 * Selector to get the inner product between the basis elements.
	 * @return m_basisProd.
	 */
	inline MatrixXd getBasisProd() const {return m_basisProd;};

	/**
	 * Selector to get the number of cluster to be achieve.
	 * @return int m_nbClust.
	 */
	inline int getNbClust() const {return m_nbClust;};

	/**
	 *  Slector to get the threeshod for the scree test
	 * @return m_thd
	 */
	inline double getThd() const {return m_thd;};

	/**
	 * Selector to get the weights matrix.
	 * @return m_weights.
	 */
	inline MatrixXd getWeights() const {return m_weights;};

	/**
	 * Selector to get the pca results.
	 * @return m_pca.
	 */
	inline funPca getPca() const {return m_pca;};

	/**
	 * Selector to get the parameters set in the mStep().
	 * @return m_paramobj.
	 */
	inline param getParamobj() const {return m_paramobj;};

	/**
	 * Selector to get the logliklihood value of the model.
	 * @return m_loglik.
	 */
	inline double getLoglik() const {return m_loglik;};

	/**
	 * Selector to get the aic criterion of the model.
	 * @return m_aic.
	 */
	inline double getAic() const {return m_aic;};


	/**
	 * Selector to get the bic criterion of the model.
	 * @return m_bic.
	 */
	inline double getBic() const {return m_bic;};


	/**
	 * Selector to get the icl of the model.
	 * @return m_icl.
	 */
	inline double getIcl() const {return m_icl;};

	/**
	 * Selector to get the clustering of the curves.
	 * @return m_clusters.
	 */
	inline VectorXi getClusters() const {return m_clusters;};

	/**
	 * selctor to get the Z matrix.
	 * @return m_Z.
	 */
	inline MatrixXi getZ() const {return m_Z;};

	/**
	 * set the member m_empty with empty in input
	 * @param empty a boolean which indicates if we have at least an empty class
	 */
	virtual void setEmpty(bool empty)=0;

	/**
	 * get the member m_empty
	 * @return m_empty which indicates if we have at least an empty class
	 */
	inline bool getEmpty() const {return m_empty;};


	/**
	 * In this function we made an weighted functional PCA.
	 * In addition, this PCA will exclude any curve whose weight is less than a threshold, which is set at 10-5.
	 * The functional PCA is one of the main methods in the IModel class.
	 * Here, we update the parameters scores and egeinvalue that permits to update in the mStep most of the members
	 * of the model.
	 * @param nHarm : number of harmonics or principal component to be retain
	 * @param weights : the weights of the curves, an example of this is in the mstep() we set the
	 * vector weights to one of the columns of m_weights (the weights matrix of the curves).
	 */
	virtual void functionalPca(int nHarm, VectorXd weights)=0;

	/**
	 * This function run the method functionalPca and in addition update the parameter m_pca.harmonics
	 * which is the coefficients of the eigen functions in the basis expansion.
	 * @param nHarm : number of harmonics or principal component to be retain
	 * @param weights : the weights of the curves, an example of this is in the mstep() we set the
	 * vector weights to one of the columns of m_weights (the weights matrix of the curves).
	 */
	virtual void mfpca(int nHarm, VectorXd weights)=0;

	/**
	 * Function providing the scree test of Cattel to determine the number of principal component \n
	 * to be retain in a functional PCA. This function update m_paramobj.r \n
	 * @param k The index of cluster
	 */
	virtual void screeTest(int k)=0;

	/**
	 * This method in addition to calling the screeTest method in order to compute the dimensions of the curves,
	 * guarantee (when it's needed) that the dimensions of the curves increases between two
	 * consecutive iterations of the algorithm.
	 * @param k the index of the cluster.
	 */
	virtual void screeTestDimIncrease(int k)=0;

	/**
	 * This function is to carry out the mStep in order to getting a parameters estimation.
	 * In the mStep() we update the member m_paramobj.
	 * This is a first version of the mStep, in this version the dimensions may vary by screeTest between two
	 * successive iteration (dimIncrease=false). We can also force them to only increase (dimIncrease=true)
	 * @param dimIncrease this is a boolean parameter to specify if we constrain the dimensions to increase or not.
	 */
	virtual void mStep(bool dimIncrease)=0;

	/**
	 * This function compute the mStpe in the case where the dimensions are fixed. So the only difference
	 * with the other version of the mStep is that here we don't make a screeTest for the dimensions we
	 * take the value inside dimFixed.
	 * @param dimFixed A vector of integer with size equal to the number of clusters.
	 * 		  The kth component of dimFixed	corresponds to the dimensions of curves in the kth cluster.
	 */
	virtual void mStep(VectorXi dimFixed)=0;

	/**
	 * This method is to calculate a normal distribution within the parameters given in input, by computing
	 * the formula of the normal distribution on a vector X with mean m and standard déviation std.
	 * @param X : the data matrix.
	 * @param mean : the mean of the normal distribution.
	 * @param std : the standard deviation of the normal distribution. Std must be positive.
	 * @return a matrix which is normally distributed with mean and standard deviation described above.
	 */
	virtual VectorXd normalDistribution(VectorXd X, double mean, double std) const=0;

	/**
	 * This function compute the eStep of the model.
	 * In the eStep() we update the member m_weights.
	 */
	virtual void eStep()=0;

	/**
	 * Here we compute the logliklihood function of the model.
	 */
	virtual void logliklihood()=0;

	/**
	 * This function compute the aic criterion.
	 */
	virtual void aic()=0;

	/**
	 * Computation of the bic criterion.
	 */
	virtual void bic()=0;

	/**
	 * The icl computation.
	 */
	virtual void icl()=0;

	/**
	 * This method updates the member "m_clusters" which represents \n
	 * the classification of curves in the differents clusters.
	 * Here we update also the member m_Z.
	 */
	virtual void updateClusters()=0;

	/**
	 * This is a first way to generate a random matrix for the weigths of the curves in order to initialize
	 * the algorithm. We generate here a matrix from the multinomial distribution with equal probabilities.
	 * Each row is built in two steps.
	 * If we denote by K the number of clusters, then the first step is to cut out the interval [0,1[
	 * into K equidistant intervals.
	 *@example if K = 3, the three intervals will be [0,1/3[, [1/3,2/3[ and [2/3,1[.
	 * The second step is to generate a number in [0,1[ from the uniform law.
	 * Thus the current row of the matrix will take one at the column with index which corresponds
	 * to the index of the interval in which the generated number belongs.
	 * @example if K=3 and the uniform number is equal to 0.4, then 1 will be in the second column because 0.4
	 * is in the second interval [1/3,2/3[.
	 */
	virtual void generateHardWeights()=0;

	/**
	 * This is another version to generate a random weight matrix.
	 * Here each row of this weight matrix will contain a probabilities, and we built each row in two steps.
	 * First we generate, in [0,1[, K-1 uniform number and we make K intervals with this numbers.
	 * The second step is to assign to each column the longer of the corresponding interval
	 * (ie column 1=longer of interval 1)
	 * @example if K=3, and for row index i, we generate 0.1 and 0.6. So we make the 3 intervals
	 * [0,0.1[, [0.1,0.6[ and [0.6,1[, finally if we denote by T the weight matrix,
	 *  then T(i,1)=0.1, T(i,2)=0.5, T(i,3)=0.4.
	 */
	virtual void generateSoftWeights()=0;



};

#endif /* IMODEL_H_ */
