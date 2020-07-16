/*
 * ACTParticle.cpp
 *
 *  Created on: 17/mar/2013
 *      Author: jonah
 */

#include "ACTParticle.h"

template <typename number> LR_vector<number> const ACTParticle<number>::principal_axis(1, 0, 0);
template <typename number> LR_vector<number> const ACTParticle<number>::second_axis(0, 0, 1);
template <typename number> LR_vector<number> const ACTParticle<number>::third_axis(0, 1, 0);

template<typename number>
ACTParticle<number>::ACTParticle() : BaseParticle<number>()  {

}

template<typename number>
ACTParticle<number>::~ACTParticle() {

}

template<typename number>
void ACTParticle<number>::add_bonded_neighbor(ACTParticle<number> *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		ParticlePair<number> new_pair(this, nn);
		this->affected.push_back(new_pair);
		nn->affected.push_back(new_pair);
	}
}


template<typename number>
bool ACTParticle<number>::is_bonded(BaseParticle<number> *q) {
	ACTParticle<number> *Cq = static_cast<ACTParticle<number> *>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

template class ACTParticle<double>;
template class ACTParticle<float>;
