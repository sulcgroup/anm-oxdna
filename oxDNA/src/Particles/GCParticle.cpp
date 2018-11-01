/*
 * GCParticle.cpp
 *
 *  Created on: 17/mag/2013
 *      Author: jonah
 */

#include "GCParticle.h"

template<typename number>
GCParticle<number>::GCParticle() : BaseParticle<number>()  {

}

template<typename number>
GCParticle<number>::~GCParticle() {

}

template<typename number>
void GCParticle<number>::add_bonded_neighbor(GCParticle<number> *nn) {
	if(!is_bonded(nn)) {
		bonded_neighs.insert(nn);
		nn->bonded_neighs.insert(this);

		//ParticlePair<number> new_pair(this, nn);
		//this->affected.push_back(new_pair);
		//nn->affected.push_back(new_pair);
	}
}

template<typename number>
bool GCParticle<number>::is_bonded(BaseParticle<number> *q) {
	GCParticle<number> *Cq = static_cast<GCParticle<number> *>(q);
	return !(bonded_neighs.find(Cq) == bonded_neighs.end());
}

template class GCParticle<double>;
template class GCParticle<float>;
