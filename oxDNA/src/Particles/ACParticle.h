/*
 * ACParticle.h
 *
 *  Created on: 29/oct/2018
 *      Author: jonah
 */

#ifndef ACParticle_H_
#define ACParticle_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
template<typename number>
class ACParticle: public BaseParticle<number> {
protected:

public:
	ACParticle();
	virtual ~ACParticle();

	virtual bool is_rigid_body() { return false; }

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void add_bonded_neighbor(ACParticle<number> *nn);

	std::set<ACParticle<number> *> bonded_neighs;
};

#endif /* ACPARTICLE_H_ */
