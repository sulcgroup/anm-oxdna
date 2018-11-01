/*
 * GCParticle.h
 *
 *  Created on: 29/oct/2018
 *      Author: jonah
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
	virtual void add_bonded_neighbor(GCParticle<number> *nn);

	std::set<GCParticle<number> *> bonded_neighs;
};

// TODO: Can I delete this??
#endif /* CUSTOMPARTICLE_H_ */
