/*
 * ACParticle.cpp
 *
 *  Created on: 17/mar/2013
 *      Author: jonah
 */

#include "ACParticle.h"

template<typename number>
ACParticle<number>::ACParticle() : BaseParticle<number>()  {

}

template<typename number>
ACParticle<number>::~ACParticle() {

}

template<typename number>
void ACParticle<number>::add_bonded_neighbor(ACParticle<number> *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		ParticlePair<number> new_pair(this, nn);
		this->affected.push_back(new_pair);
		nn->affected.push_back(new_pair);
	}
}

template<typename number>
bool ACParticle<number>::is_bonded(BaseParticle<number> *q) {
	ACParticle<number> *Cq = static_cast<ACParticle<number> *>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

template class ACParticle<double>;
template class ACParticle<float>;
