#include "IAlgo.h"

IAlgo::IAlgo() {
	// TODO Auto-generated constructor stub
	m_model=0;

}

IAlgo::IAlgo(IModel* model){
	m_model=model;
}

IAlgo::~IAlgo() {
	// TODO Auto-generated destructor stub
}

