/*
 * GCParticle.h
 *
 *  Created on: 17/mag/2013
 *      Author: lorenzo
 */

#ifndef GCPARTICLE_H_
#define GCPARTICLE_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A customisable particle. Used by GCInteraction.
 */
template<typename number>
class GCParticle: public BaseParticle<number> {
protected:

public:
	GCParticle();
	virtual ~GCParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void add_bonded_neigh(GCParticle<number> *nn);

	std::set<GCParticle<number> *> bonded_neighs;
};

#endif /* CUSTOMPARTICLE_H_ */
