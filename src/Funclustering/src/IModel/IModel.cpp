#include "IModel.h"


/// Default constructor
IModel::IModel() {
	// TODO Auto-generated constructor stub
	m_nbClust=0;
	m_thd=0;

	// this MatrixXd pointer is initialized to null pointer.
	m_paramobj.scores=0;

	//m_loglik initialized into -Inf.
	m_loglik=- numeric_limits<double>::infinity();
	// The dimensions are initialized to zero

	m_paramobj.r=VectorXi::Constant(m_nbClust,0);

	/*
	 the three below criterion are also initialized to -Inf beacause all of them depends on loglik and loglik
	 is initialized to -Inf.
	 */
	m_aic=- numeric_limits<double>::infinity();
	m_bic=- numeric_limits<double>::infinity();
	m_icl=- numeric_limits<double>::infinity();


	m_empty=false;
}

/*
 * This constructor construct a model by initializing it's members by the parameters given in input.
 * In fact this two first input parameters will be provided by R using the fda package.
 * @param coefs : the matirx which sets the m_coefs member.
 * @param basisProd : sets the matrix of the basis function inner product.
 * @param nbClust : the number of clusters.
 * @param thd : the threshold for the scree test.
 */
IModel::IModel(MatrixXd coefs,MatrixXd basisProd, int nbClust, double thd) {
	m_coefs=coefs;
	m_basisProd=basisProd;
	m_nbClust=nbClust;
	m_paramobj.scores = new MatrixXd[nbClust];
	m_thd=thd;
	m_loglik = - numeric_limits<double>::infinity();/* m_loglik initialized into -Inf */
	m_paramobj.r=VectorXi::Constant(m_nbClust,0);
	m_aic=- numeric_limits<double>::infinity();
	m_bic=- numeric_limits<double>::infinity();
	m_icl=- numeric_limits<double>::infinity();

	m_empty=false;
}

/*
 * This is copy constructor and it lead in particular to make a deep copy of the object
 * when it's needed. we assign to each member of *this it's value in modelToCopy.
 * @param modelToCopy an model to be copied in *this.
 */
IModel::IModel(IModel const& modelToCopy){
	// assign all not pointer members
	m_coefs=modelToCopy.m_coefs;
	m_basisProd=modelToCopy.m_basisProd;
	m_nbClust=modelToCopy.m_nbClust;
	m_paramobj.prop=modelToCopy.m_paramobj.prop;
	m_paramobj.m=modelToCopy.m_paramobj.m;
	m_paramobj.V=modelToCopy.m_paramobj.V;
	m_paramobj.r=modelToCopy.m_paramobj.r;
	m_paramobj.scores=0;
	m_pca=modelToCopy.m_pca;
	m_weights=modelToCopy.m_weights;
	m_loglik=modelToCopy.m_loglik;
	m_aic=modelToCopy.m_aic;
	m_bic=modelToCopy.m_bic;
	m_icl=modelToCopy.m_icl;
	m_thd=modelToCopy.m_thd;
	m_clusters=modelToCopy.m_clusters;
	m_empty=modelToCopy.m_empty;

	/*
	  allocate memory for the pointer "m_paramobj.scores" to make a deep copy of the member
	  modelToCopy.m_paramobj.scores if it is a not null pointer.
	 */
	if (modelToCopy.m_paramobj.scores){
		// allocate memory for m_paramobj.scores
		m_paramobj.scores=new MatrixXd [modelToCopy.m_nbClust];
		//deep copy
		for (int i=0;i<modelToCopy.m_nbClust;i++){
			m_paramobj.scores[i]=modelToCopy.m_paramobj.scores[i];
		}
	}
}

/*
 * This is the overload of the assignment operator and it help to make a deep copy when it's needed.
 * @param modelToCopy : model to be copied in *this.
 * @return *this.
 */
IModel& IModel::operator =(IModel const& modelToCopy){
	//testing if the to models are different, to avoid make computation when we make statement like model=model.
	if (this != &modelToCopy){
		// assign all not pointer members
		m_coefs=modelToCopy.m_coefs;
		m_basisProd=modelToCopy.m_basisProd;
		m_nbClust=modelToCopy.m_nbClust;
		m_paramobj.prop=modelToCopy.m_paramobj.prop;
		m_paramobj.m=modelToCopy.m_paramobj.m;
		m_paramobj.V=modelToCopy.m_paramobj.V;
		m_paramobj.r=modelToCopy.m_paramobj.r;
		m_pca=modelToCopy.m_pca;
		m_weights=modelToCopy.m_weights;
		m_loglik=modelToCopy.m_loglik;
		m_aic=modelToCopy.m_aic;
		m_bic=modelToCopy.m_bic;
		m_icl=modelToCopy.m_icl;
		m_thd=modelToCopy.m_thd;
		m_clusters=modelToCopy.m_clusters;
		m_empty=modelToCopy.m_empty;

		//delete the pointer m_paramobj.scores befor making copy of modelToCopy.m_paramobj.scores.
		if (m_paramobj.scores){
			delete [] m_paramobj.scores;
			m_paramobj.scores=0;
		}
		// allocate new memory and make a deep copy of the pointer modelToCopy.m_paramobj.scores.
		m_paramobj.scores=new MatrixXd [modelToCopy.m_nbClust];
		for (int i=0;i<modelToCopy.m_nbClust;i++){
			m_paramobj.scores[i]=modelToCopy.m_paramobj.scores[i];
		}
	}
	return *this;
}

/*
 *  Here we destruct the IModel object by making a delete statement on the only pointer in members
 *  of the IModel object, ie: m_paramobj.scores. The other members are eigen object and the eigen destructor
 *  will be called automatically.
 */
IModel::~IModel() {
	// TODO Auto-generated destructor stub
	if (m_paramobj.scores){
		delete [] m_paramobj.scores;
		m_paramobj.scores=0;
	}
}


// Set weights matrix to the matrix given in input.
void IModel::setWeights(MatrixXd const& weights){
	m_weights=weights;
}

/*
 * This method set the liklihood to -inf. we call it after the initialization to restart the EM algorithm
 * with the smallest liklihood possible and thus we make sure that the liklihood will increase.
 */
void IModel::setLoglikToInfinity(){
	m_loglik=-numeric_limits<double>::infinity();
}

/*
 * This method set the dimensions to zero in each cluster. we call it after the initialization
 * to restart the EM algorithm with the basic parameters.
 */
void IModel::setDimensionsToZero(){
	m_paramobj.r=VectorXi::Zero(m_paramobj.r.size());
}
