/*
 * Model.cpp
 *
 *  Created on: 21 juin 2013
 *      Author: soueidat
 */

#include "Model.h"


Model::Model():IModel() {
	// TODO Auto-generated constructor stub


}

Model::~Model() {
	// TODO Auto-generated destructor stub
}

/*
 *The initialization constructor which calls the initialization constructor of the IModel class.
 * @param coefs : The coefficients matrix.
 * @param basisProd : The matrix of the inner product between the basis functions.
 * @param nbClust : The number of clusters.
 * @param thd : The threshold  to consider in the scree test.
 */
Model::Model(MatrixXd coefs, MatrixXd basisProd, int nbClust, double thd)
:IModel(coefs,basisProd,nbClust,thd){}

/*
 * The copy constructor which calls the copy constrctor of the IModel class.
 * @param modelToCopy the model to copy in *this.
 */
Model::Model(Model const& modelToCopy):IModel(modelToCopy){


}

typedef std::pair<double, int> mypair;
bool comparator (const mypair& l, const mypair& r)
   {return l.first > r.first;}

/*
 * In this function we made an weighted functional PCA.
 * In addition, this PCA will exclude any curve whose weight is less than a threshold, which is set at 10-5.
 * The functional PCA is one of the main methods in the IModel class.
 * Here, we update the parameters scores and egeinvalue that permits to update in the mStep most of the members
 * of the model.
 * @param nHarm : number of harmonics or principal component to be retain
 * @param weights : the weights of the curves, an example of this is in the mstep() we set the
 * vector weights to one of the columns of m_weights (the weights matrix of the curves).
 */
void Model::functionalPca(int nHarm, VectorXd weights) {
/*
  // fix empty class erreur if the weigths is null for one class,
	//then we make an equal weights for the curves in this class
	if((weights.array() < 1e-5).all()){
		int n=weights.rows();
		double tk=double(1)/double(n);
		weights=VectorXd::Constant(n,tk);
		//cout<<"weights const"<<endl<<weights.transpose()<<endl;
	}
*/
	//normalization
	if (weights.sum() != double(1)) {
		weights=(weights)/(weights.sum());
	}
	// first we keep only the curves with acceptable weights, ie grater than 10e-5.

	//coefsAccept and weightsAccept stores the coefs and weights of curves with acceptable weights.
	MatrixXd coefsAccept(m_coefs.rows(),m_coefs.cols());
	VectorXd weightsAccept(weights.rows());

	//AcceptIndex will take the last index of the curves with acceptable weights
	int AcceptIndex=0;
	for (int i=0; i<weights.rows(); i++){
		// test for each curves if we have an acceptable weights
		if (weights(i) > 10e-5){
			//stores the weights and coefs of the curves with acceptable weights
			coefsAccept.row(AcceptIndex)=m_coefs.row(i);
			weightsAccept(AcceptIndex)=weights(i);
			//update AcceptIndex
			AcceptIndex +=1;
		}
	}
	// Keep only curves with acceptable weights
	if (AcceptIndex==weights.rows()) {
		weightsAccept=weights;
		coefsAccept=m_coefs;
	}
	else {
		coefsAccept.conservativeResize(AcceptIndex,m_coefs.cols());
		weightsAccept.conservativeResize(AcceptIndex);
	}

	int n=coefsAccept.rows();
	MatrixXd matcent;
	matcent=MatrixXd::Identity(n,n)-(MatrixXd::Constant(n,1,1)*weightsAccept.transpose());

	MatrixXd A=matcent*coefsAccept;

	MatrixXd V;
	V=A*m_basisProd*A.transpose()*weightsAccept.asDiagonal();

	// eigen decomposition
	EigenSolver<MatrixXd> eV(V);
	MatrixXd c=eV.eigenvectors().real();//eigen vectors
  VectorXd val=eV.eigenvalues().real();//eigenvalue
  VectorXd vp=val/val.sum();//varprop
  
  std::vector<mypair> indexVar(val.size());
  // sorting eigenvalues
  for (int i = 0; i < vp.size(); ++i)
  {
    indexVar[i] = mypair(vp(i), i);
  }
  std::sort(indexVar.begin(), indexVar.end(), comparator);
  // computing reordered eigenvalues and eigenvectors
  VectorXd val_ord(val.size());
  VectorXd vp_ord(vp.size());
  MatrixXd c_ord(c.rows(), c.cols());
  for (int i = 0; i < indexVar.size(); ++i)
  {
    val_ord  [i] = val  [indexVar[i].second];
    vp_ord   [i] = vp   [indexVar[i].second];
    c_ord.col(i) = c.col(indexVar[i].second);
  }

  // statement val.head(nHarm) requires that val.size() > nHarm
  if(val.size()>nHarm){
  	//we keep the number eigenvalues requested
  	m_pca.values =val_ord.head(nHarm);//values
  	m_pca.varprop=vp_ord .head(nHarm);
  }
  //else: we take all eigenvalue (so in this case we will have less than the eigenvalues requested)
  else{
  	m_pca.values =val_ord;
  	m_pca.varprop=vp_ord ;
  }

  // if the number of columns of c1 is greater than nHarm we keep the required nHarm scores
   if(c_ord.cols() > nHarm){
 		m_pca.scores = c_ord.block(0, 0, c_ord.rows(), nHarm);
 	}
 	else{ // we take all c1 and in this case we can have a number of harmonics less than the requested
 		m_pca.scores = c_ord;
 	}

/* Since we keep in the eigen decomposition only curves with acceptable weights,
   the dimension of the score matrix may be incoherent with the subsequent computation,
   especially the number of rows. In fact in the mStep we work with the columns of the matrix scores
   to make some product with another vector or columns with size equal to m_coefs.rows().
   So we have to resize the number of rows of scores matrix by adding the required rows=zeros.
  */
  if (c_ord.rows() != m_coefs.rows()){
  	m_pca.scores.conservativeResize(m_coefs.rows(),c.cols());
  	for (int i=c_ord.rows()-1; i<m_coefs.rows(); i++){
  		m_pca.scores.row(i).setZero();
  	}
  }
}

/*
 * This function run the method functionalPca and in addition update the parameter m_pca.harmonics
 * which is the coefficients of the eigen functions in the basis expansion.
 * @param nHarm : number of harmonics or principal component to be retain
 * @param weights : the weights of the curves, an example of this is in the mstep() we set the
 * vector weights to one of the columns of m_weights (the weights matrix of the curves).
 */
void Model::mfpca(int nHarm, VectorXd weights){
//  from here to the harmonincs calculation; this is the same code of functionalPca  //

  //normalization
	if (weights.sum() != double(1)) {
		weights=(weights)/(weights.sum());
	}
	// first we keep only the curves with acceptable weights, ie grater than 10e-5.

	//coefsAccept and weightsAccept stores the coefs and weights of curves with acceptable weights.
	MatrixXd coefsAccept(m_coefs.rows(),m_coefs.cols());
	VectorXd weightsAccept(weights.rows());

	//AcceptIndex will take the last index of the curves with acceptable weights
	int AcceptIndex=0;
	for (int i=0; i<weights.rows(); i++){
		// test for each curves if we have an acceptable weights
		if (weights(i) > 10e-5){
			//stores the weights and coefs of the curves with acceptable weights
			coefsAccept.row(AcceptIndex)=m_coefs.row(i);
			weightsAccept(AcceptIndex)=weights(i);
			//update AcceptIndex
			AcceptIndex +=1;
		}
	}
	// Keep only curves with acceptable weights
	if (AcceptIndex==weights.rows()) {
		weightsAccept=weights;
		coefsAccept=m_coefs;
	}
	else {
		coefsAccept.conservativeResize(AcceptIndex,m_coefs.cols());
		weightsAccept.conservativeResize(AcceptIndex);
	}

	int n=coefsAccept.rows();
	MatrixXd matcent;
	matcent=MatrixXd::Identity(n,n)-(MatrixXd::Constant(n,1,1)*weightsAccept.transpose());

	MatrixXd A=matcent*coefsAccept;

	MatrixXd V;
	V=A*m_basisProd*A.transpose()*weightsAccept.asDiagonal();

	// eigen decomposition
	EigenSolver<MatrixXd> eV(V);
	MatrixXd c=eV.eigenvectors().real();//eigen vectors
  VectorXd val=eV.eigenvalues().real();//eigenvalue
  VectorXd vp=val/val.sum();//varprop
  
  std::vector<mypair> indexVar(val.size());
  // sorting eigenvalues
  for (int i = 0; i < vp.size(); ++i)
  {
    indexVar[i] = mypair(vp(i), i);
  }
  std::sort(indexVar.begin(), indexVar.end(), comparator);
  // computing reordered eigenvalues and eigenvectors
  VectorXd val_ord(val.size());
  VectorXd vp_ord(vp.size());
  MatrixXd c_ord(c.rows(), c.cols());
  for (int i = 0; i < indexVar.size(); ++i)
  {
    val_ord  [i] = val  [indexVar[i].second];
    vp_ord   [i] = vp   [indexVar[i].second];
    c_ord.col(i) = c.col(indexVar[i].second);
  }

  // statement val.head(nHarm) requires that val.size() > nHarm
  if(val.size()>nHarm){
  	//we keep the number eigenvalues requested
  	m_pca.values =val_ord.head(nHarm);//values
  	m_pca.varprop=vp_ord .head(nHarm);
  }
  //else: we take all eigenvalue (so in this case we will have less than the eigenvalues requested)
  else{
  	m_pca.values =val_ord;
  	m_pca.varprop=vp_ord ;
  }

  // if the number of columns of c1 is greater than nHarm we keep the required nHarm scores
 	if(c_ord.cols() > nHarm){
 		m_pca.scores = c_ord.block(0, 0, c_ord.rows(), nHarm);
 	}
 	else{ // we take all c1 and in this case we can have a number of harmonics less than the requested
 		m_pca.scores = c_ord;
 	}

/* Since we keep in the eigen decomposition only curves with acceptable weights,
   the dimension of the score matrix may be incoherent with the subsequent computation,
   especially the number of rows. In fact in the mStep we work with the columns of the matrix scores
   to make some product with another vector or columns with size equal to m_coefs.rows().
   So we have to resize the number of rows of scores matrix by adding the required rows=zeros.
  */
  if (c_ord.rows() != m_coefs.rows()){
  	m_pca.scores.conservativeResize(m_coefs.rows(),c.cols());
  	for (int i=c_ord.rows()-1; i<m_coefs.rows(); i++){
  		m_pca.scores.row(i).setZero();
  	}
  }

//      End of the code of functional pca           //
  // harmonics caluculation
    MatrixXd V1;
    V1=A.transpose()*weightsAccept.asDiagonal()*A*m_basisProd;

    // eigen decomposition
    EigenSolver<MatrixXd> eV1(V1);
   	MatrixXd c1=eV1.eigenvectors().real();//eigen vectors
   	VectorXd vp1=eV1.eigenvalues().real();

    std::vector<mypair> indexVar1(vp1.size());
    // sorting eigenvalues
    for (int i = 0; i < vp1.size(); ++i)
    {
      indexVar1[i] = mypair(vp1(i), i);
    }
    std::sort(indexVar1.begin(), indexVar1.end(), comparator);
    // computing reordered eigenvalues and eigenvectors
    VectorXd vp1_ord(vp1.size());
    MatrixXd c1_ord(c1.rows(), c1.cols());
    for (int i = 0; i < indexVar1.size(); ++i)
    {
      vp1_ord   [i] = vp1   [indexVar1[i].second];
      c1_ord.col(i) = c1.col(indexVar1[i].second);
    }    

   	// if the number of columns of c1 is greater than nHarm we keep the required nHarm harmonics
   	if(c1_ord.cols() > nHarm){
   		m_pca.harmonics = c1_ord.block(0, 0, c1_ord.rows(), nHarm);
   	}
   	else{ // we take all c1 and in this case we can have a number of harmonics less than the requested
   		m_pca.harmonics = c1_ord;
   	}
}

/*
 * Function providing the scree test of Cattel to determine the number of principal component \n
 * to be retain in a functional PCA. This function update m_paramobj.r \n
 * @param k The index of cluster
 */
void Model::screeTest(int k){

	VectorXd val(m_pca.values.size());
	VectorXd sc(m_pca.values.size()-1);

	// the eigenvalues
	val=m_pca.values;

	// replace the negatives values by 0.
	for (int i=0;i<val.size();i++){
		if (val(i)<0) {val(i)=0;}
	}

	// sc is the vector of the differences between two consecutive eigenvalues
	for (int i=0;i<sc.size();i++){
		sc(i)=val(i+1)-val(i);
	}

	// The bellow statement replace each value of sc with its absolute value
	sc=sc.cwiseAbs();

	// r is the result of the screeTest, and by the way sc.size() is the maximum number of principal component
		//(or dimensions of the curves) to be selected.
	int r=sc.size();
	VectorXd tmp;
	for (int j = 1; j <= sc.size()-1; ++j) {
		//the last (sc.size()-j) values
		tmp=sc.tail(sc.size()-j);

		// if all coefficients of tmp, are less than the threshold the r will take the current value of j
		if ((tmp.maxCoeff() < m_thd*sc.maxCoeff())){
			r=j;
			break;
		}
	}
	// the dimensions of the curves in the kth cluster, m_paramobj.r(k) equal to the result of the scrreTest
	m_paramobj.r(k)=r;
}

/*
 * This method in addition to calling the screeTest method in order to compute the dimensions of the curves,
 * guarantee (when it's needed) that the dimensions of the curves increases between two
 * consecutive iterations of the algorithm.
 * @param k the index of the cluster.
 */
void Model::screeTestDimIncrease(int k){
	// the previous values of the dimension
	// fix dim, the dimensions can never be grater than m_pca.scores.cols().
	int rmax=m_pca.scores.cols();
	int rk_prev=min(m_paramobj.r(k),rmax);

	// update dimensions
	screeTest(k);

	// Constraint the dimensions to increase, if the previous dimension is greater than the current one,
		// then keep the previous dimensions
	if (rk_prev > m_paramobj.r(k)){m_paramobj.r(k)=rk_prev;}

}


/*
 * This function is to carry out the mStep in order to getting a parameters estimation.
 * In the mStep() we update the member m_paramobj.
 * This is a first version of the mStep, in this version the dimensions may vary by screeTest between two
 * successive iteration (dimIncrease=false). We can also force them to only increase (dimIncrease=true)
 * @param dimIncrease this is a boolean parameter to specify if we constrain the dimensions to increase or not.
 */
void Model::mStep(bool dimIncrease) {
	// an initial number of principal components, and will take the minimum between
		// the number of basis functions (m_coefs.cols(), and the number of curves (m_coefs.rows()).
	int r_init=min(m_coefs.cols(),m_coefs.rows());

	// the number of cluster
	int K=m_nbClust;

	// Proportions of each cluster in classification (prop) is the mean of the corresponding column in the
		// weights matrix
	m_paramobj.prop=m_weights.colwise().mean();

	//Initialisation of m and V
	m_paramobj.V.setZero(K,r_init);
	m_paramobj.m.setZero(K,r_init);
	for (int k=0;k<K;k++){
		// Functional PCA
		functionalPca(r_init, m_weights.col(k));

		//Dimension estimation, if dimIncrease=True call the version of screeTest with allow to increase the
			//dimension, else call the standard version of screeTest
		if (dimIncrease){
			screeTestDimIncrease(k);
		}
		else
			screeTest(k);

		// Scores Update
		m_paramobj.scores[k]=m_pca.scores.block(0,0,m_coefs.rows(),m_paramobj.r(k));

		//Weighted mean and weighted covariance operator calculation
		VectorXd tmp;
		for (int j=0;j<m_paramobj.r(k);j++){
			/* m(k,j) is the mean of the jth column of the kth score matrix weighted
			   by the kth column of the weight matrix
			 */

			// here we will weight the jth column of the score[k] by the kth column of the weights matrix
			tmp=m_paramobj.scores[k].col(j).cwiseProduct(m_weights.col(k));

			// m(k,j) is the sum of the previous product divided by the sum of the kth column of the weights
			m_paramobj.m(k,j)=(tmp.sum())/(m_weights.col(k).sum());


			/*
			 V(k,j) is the covariance of the jth column of the scores[k] weighted by the kth column of the weight
			 matrix, m.cwiseAbs2()= square value of all coefficients, the array() function is called
			 to make possible operation like (matrix - double)
			 */

			// we take square of the jth column of the score matrix minus it's mean.
			tmp=((m_paramobj.scores[k].col(j).array())-m_paramobj.m(k,j)).cwiseAbs2();

			// here we will weight the previous tmp by the kth column of the weights matrix
			tmp=tmp.cwiseProduct(m_weights.col(k));

			// weighted covariance is the sum of the weighted vector divided by the sum of the weights
			m_paramobj.V(k,j)=(tmp.sum())/(m_weights.col(k).sum());
		}
	}
}


/*
 * This function compute the mStpe in the case where the dimensions are fixed. So the only difference
 * with the other version of the mStep is that here we don't make a screeTest for the dimensions we
 * take the value inside dimFixed.
 * @param dimFixed A vector of integer with size equal to the number of clusters.
 * 		  The kth component of dimFixed	corresponds to the dimensions of curves in the kth cluster.
 */
void Model::mStep(VectorXi dimFixed){
	//Dimension are fixed with values in dimFixed
	m_paramobj.r=dimFixed;

	// an initial number of principal components, and will take the minimum between
		// the number of basis functions (m_coefs.rows(), and the number of curves (m_coefs.cols()).
	int r_init=min(m_coefs.cols(),m_coefs.rows());

	// the number of cluster
	int K=m_nbClust;

	// Proportions of each cluster in classification (prop) is the mean of the corresponding column in the
		// weights matrix
	m_paramobj.prop=m_weights.colwise().mean();

	//Initialisation of m and V
	m_paramobj.V.setZero(K,r_init);
	m_paramobj.m.setZero(K,r_init);
	for (int k=0;k<K;k++){
		// Functional PCA
		functionalPca(r_init, m_weights.col(k));

		// Scores Update
		m_paramobj.scores[k]=m_pca.scores.block(0,0,m_coefs.rows(),m_paramobj.r(k));

		//Weighted mean and weighted covariance operator calculation
		VectorXd tmp;
		for (int j=0;j<m_paramobj.r(k);j++){
			/* m(k,j) is the mean of the jth column of the kth score matrix weighted
			   by the kth column of the weight matrix
			 */
			// here we will weight the jth column of the score[k] by the kth column of the weights matrix
			tmp=m_paramobj.scores[k].col(j).cwiseProduct(m_weights.col(k));
			// m(k,j) is the sum of the previous product divided by the sum of the kth column of the weights
			m_paramobj.m(k,j)=(tmp.sum())/(m_weights.col(k).sum());


			/*
			 V(k,j) is the covariance of the jth column of the scores[k] weighted by the kth column of the weight
			 matrix, m.cwiseAbs2()= square value of all coefficients, the array() function is called
			 to make possible operation like (matrix - double)
			 */
			// we take square of the jth column of the score matrix minus it's mean.
			tmp=((m_paramobj.scores[k].col(j).array())-m_paramobj.m(k,j)).cwiseAbs2();
			// here we will weight the previous tmp by the kth column of the weights matrix
			tmp=tmp.cwiseProduct(m_weights.col(k));
			// weighted covariance is the sum of the weighted vector divided by the sum of the weights
			m_paramobj.V(k,j)=(tmp.sum())/(m_weights.col(k).sum());
		}
	}
}


/*
 * This method is to calculate a normal distribution within the parameters given in input, by computing
 * the formula of the normal distribution on a vector X with mean m and standard déviation std.
 * @param X : the data matrix.
 * @param mean : the mean of the normal distribution.
 * @param std : the standard deviation of the normal distribution. Std must be positive.
 * @return a matrix which is normally distributed with mean and standard deviation described above.
 */
VectorXd Model::normalDistribution(VectorXd X, double mean, double std) const{
	VectorXd normalData(X.rows());
	for (int i=0;i<X.rows();i++){
		double x=X(i);
		// comuting the formula of normal distribution
		normalData(i)=(1/(std*sqrt(2*M_PI)))*exp(-pow(x-mean,2)/(2*pow(std,2)));
	}
	return normalData;
}


/*
 * This function compute the eStep of the model.
 * In the eStep() we update the member m_weights.
 */
void Model::eStep(){
	m_weights=m_paramobj.prop.colwise().replicate(m_weights.rows());//supprime
	VectorXd tmp;//(m_paramobj.scores[k].col(j).rows());//allouer une 1ere à l'extérieur
	for (int k=0;k<m_nbClust;k++){

		for (int j=0;j<m_paramobj.r(k);j++){//j=1
			double v=m_paramobj.V(k,j);

			// fix the probleme of empty class : we have to compute the formulas in the following "else" but
			//if v is nul then the std in input of normalDistribution will be =0, hence we will have a nan,
			// so we make tmp=mean (variance=0 => v.a constante=mean), or mean=0 => tmp=0 => m_weights.col(k)=0
			if (v==0.0){
				m_weights.col(k)=VectorXd::Constant(m_coefs.rows(),0.0);
				break;
			}
			else{
				tmp=normalDistribution(m_paramobj.scores[k].col(j),m_paramobj.m(k,j),sqrt(v));
				m_weights.col(k)=m_weights.col(k).cwiseProduct(tmp);
			}
		}
		// if tik is equal to zero for column(k), the the k-th class is empty
		if((m_weights.col(k).array()==0.0).all()){m_empty=true;}
	}

	// check for infinty tik : if tik=Inf then we replace by tik=1 and tij=0 for j != k
	for(int i=0;i<m_weights.rows();i++){
		for(int k=0;k<m_nbClust;k++){
			if(fabs(m_weights(i,k))==numeric_limits<double>::infinity()){
				for(int j=0;j<m_nbClust;j++){m_weights(i,j)=0.0;}
				m_weights(i,k)=1.0;
				break;
			}
		}
	}

	//normalization
	for (int i=0;i<m_weights.rows();i++){
		m_weights.row(i)=m_weights.row(i).array()/m_weights.rowwise().sum()(i);
	}
	// update proportions
	m_paramobj.prop=m_weights.colwise().mean();
}


/// logliklihood

void Model::logliklihood(){
	double l_prev=m_loglik;
	MatrixXd L(m_coefs.rows(),m_nbClust);
	VectorXd tmp;
	L=m_paramobj.prop.colwise().replicate(L.rows());
	for (int k=0;k<m_nbClust;k++){
		for (int j=0;j<m_paramobj.r(k);j++){
			double v=m_paramobj.V(k,j);

			tmp=normalDistribution(m_paramobj.scores[k].col(j),m_paramobj.m(k,j),sqrt(v));
			L.col(k)=L.col(k).cwiseProduct(tmp);
	/*
			// fix the probleme of empty class : we have to compute the formulas in the following else but
			//if v is nul then the std in input of normalDistribution will be =0, hence we will have a nan,
			// so we make tmp=mean (variance=0 => v.a constante=mean), or mean=0 => tmp=0
			if (v==0){
				L.col(k)=VectorXd::Constant(L.rows(),0.0);
				break;
			}
			else{
				tmp=normalDistribution(m_paramobj.scores[k].col(j),m_paramobj.m(k,j),sqrt(v));
				L.col(k)=L.col(k).cwiseProduct(tmp);
			}
	*/
		}
	}
	ArrayXd l;
	l=(L.rowwise().sum()).array().log();
	m_loglik=l.sum();


	// the bellow test will return true in the special case of m_loglik is not a number
		// in this case we reinitialize m_loglik with -infinity
	if (m_loglik !=m_loglik) {
		m_loglik=- numeric_limits<double>::infinity() ;
	}
	// if m_loglik take = +Inf then we replace m_loglik by the last m_loglik calculated
	if (m_loglik==numeric_limits<double>::infinity()){m_loglik=l_prev;}

}


void Model::aic(){
	int nbParams=m_nbClust -1 + m_paramobj.r.sum();
	m_aic=2*(m_loglik - nbParams);
}

void Model::bic(){
	int nbParams=m_nbClust -1 + m_paramobj.r.sum();
	m_bic=2*m_loglik - (nbParams*log(m_coefs.rows()));
}

void Model::icl(){
	bic();
	MatrixXd Z=m_Z.cast<double>();
	Z=m_weights.cwiseProduct(Z) - Z;
	Z=Z.array()+ 1;
	Z=Z.array().log();
	MatrixXd Z1=m_Z.cast<double>();
	Z=Z1.cwiseProduct(Z);
	m_icl = m_bic + Z.sum();
}

/// Clustering
void Model::updateClusters(){
	//Update m_cls
	MatrixXf::Index maxIndex;
	m_clusters=VectorXi(m_coefs.rows());
	for (int i=0;i<m_weights.rows();i++){
		m_weights.row(i).colwise().sum().maxCoeff(&maxIndex);
		m_clusters(i)=maxIndex+1;
	}

	// check for empty class : a k-th cluster is empty when it's index k is an absent value in m_clusters
	for (int k = 1; k <= m_nbClust; ++k) {
		if ((m_clusters.array() != k).all()) {
			m_empty=true;
			break;
		}
	}

	//Update m_Z
	m_Z=MatrixXi::Zero(m_coefs.rows(),m_nbClust);
	for (int i = 0; i < m_weights.rows(); i++) {
		for (int k = 0; k < m_nbClust; k++) {
			if (m_clusters(i)==k+1) {
				m_Z(i,k)=1;
			}
		}
	}
	// Update proportions according to the last weights matrix calculated
	m_paramobj.prop=m_weights.colwise().mean();


}

void Model::setEmpty(bool empty){
	m_empty=empty;
}

/*
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
void Model::generateHardWeights(){
	srand(time(NULL));
	int K=m_nbClust;
	int n=m_coefs.rows();
	m_weights=MatrixXd::Zero(n,K);
	int rand1;
	for (int i=0; i< n; i++){
		rand1=rand();
		double uniformNumber=double(rand1)/double(RAND_MAX);
	    for (int j = 0; j < K; ++j) {
		 	double prob=double(j+1)/double(K);
		    if ( uniformNumber < prob) {
		    	m_weights(i,j)=1.0;
		  		break;
			}
		}
	}
}


/*
 * This is another version to generate a random weight matrix.
 * Here each row of this weight matrix will contain a probabilities, and we built each row in two steps.
 * First we generate, in [0,1[, K-1 uniform number and we make K intervals with this numbers.
 * The second step is to assign to each column the longer of the corresponding interval
 * (ie column 1=longer of interval 1)
 * @example if K=3, and for row index i, we generate 0.1 and 0.6. So we make the 3 intervals
 * [0,0.1[, [0.1,0.6[ and [0.6,1[, finally if we denote by T the weight matrix,
 *  then T(i,1)=0.1, T(i,2)=0.5, T(i,3)=0.4.
 */
void Model::generateSoftWeights(){
	srand(time(NULL));
	int K=m_nbClust;
	int n=m_coefs.rows();
	int rand1;
	m_weights=MatrixXd(n,K);
	vector<double> unifNumbers(K+1);
	unifNumbers[0]=0.0;
	unifNumbers[K]=1.0;
	for (int i = 0; i < n; ++i) {
		for (int j = 1; j < K; ++j) {
			rand1=rand();
			unifNumbers[j]=((double) rand1 / (RAND_MAX));
		}
		stable_sort(unifNumbers.begin(),unifNumbers.end());
		for (int j = 0; j < K; ++j) {
			m_weights(i,j)=unifNumbers[j+1]-unifNumbers[j];
		}
	}
}



