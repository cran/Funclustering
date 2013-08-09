#include "EMAlgo.h"

EMAlgo::EMAlgo() {
	// TODO Auto-generated constructor stub
}

EMAlgo::EMAlgo(IModel* model) :IAlgo(model) {}


EMAlgo::~EMAlgo() {
	// TODO Auto-generated destructor stub
}

void EMAlgo::run(bool dimIncrease){
	m_model -> mStep(dimIncrease);
	m_model -> eStep();
	m_model -> updateClusters();
	m_model -> logliklihood();
}

void EMAlgo::run(VectorXi dimFixed){
	m_model -> mStep(dimFixed);
	m_model -> eStep();
	m_model -> updateClusters();
	m_model -> logliklihood();
}
