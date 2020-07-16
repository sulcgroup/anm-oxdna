/*
 * ACTParticle.h
 *
 *  Created on: 29/oct/2018
 *      Author: jonah
 */

#ifndef ACTParticle_H_
#define ACTParticle_H_

#include "BaseParticle.h"
#include <set>

/**
 * @brief A customisable particle. Used by ACInteraction.
 */
template<typename number>
class ACTParticle: public BaseParticle<number> {
protected:

public:
	ACTParticle();
	virtual ~ACTParticle();

    const static LR_vector<number> principal_axis;
    const static LR_vector<number> second_axis;
    const static LR_vector<number> third_axis;

	virtual bool is_rigid_body() { return True; }

	virtual bool is_bonded(BaseParticle<number> *q);
	virtual void add_bonded_neighbor(ACTParticle<number> *nn);
	std::set<ACTParticle<number> *> bonded_neighs;
};

#endif /* ACTPARTICLE_H_ */
