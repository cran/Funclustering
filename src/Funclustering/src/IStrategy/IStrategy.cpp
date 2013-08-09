#include "IStrategy.h"

IStrategy::IStrategy() {
	// TODO Auto-generated constructor stub
	m_algo=0;
	m_epsilon=0.0;
	m_nbInit=0;
}

IStrategy::IStrategy(IAlgo* algo, double epsilon, int nbInit){
	m_algo=algo;
	m_epsilon=epsilon;
	m_nbInit=nbInit;
}
IStrategy::~IStrategy() {
	// TODO Auto-generated destructor stub
}
