/*
 * GCInteraction.h
 *
 *  Created on: 13/Nov/18
 *      Author: lorenzo
 */

#ifndef GC_INTERACTION_H_
#define GC_INTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Handles (generalised) Gaussian-Chain interactions between spheres of size .1573 simulation units
 *
 *
 * This interaction is selected with
 * interaction_type = GC
 *
 * Input options:
 *
 * @verbatim
 * none
 */
template <typename number>
class GCInteraction: public BaseInteraction<number, GCInteraction<number> > {
protected:

	number _k; //stiffness of the spring
	number _r; //radius of alpha carbon of amino acid

	number _sigma, _rstar, _b, _rcut, STRENGTH;

	//inline number _exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	inline number _spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces );
	//inline number _repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);

public:
	enum {
		SPRING_POTENTIAL = 0,
		EXC_VOL = 1

	};

	GCInteraction();
	virtual ~GCInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
};

/*
template<typename number>
number GCInteraction<number>::_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(_rcut)) {
		if(rnorm > SQR(_rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - _rcut;
			energy = EXCL_EPS * _b * SQR(rrc);
			if(update_forces) force = -r * (2 * EXCL_EPS * _b * rrc / rmod);
		}
		else {
			number tmp = SQR(_sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}
*/
/*
template<typename number>
number GCInteraction<number>::_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	LR_vector<number> force(0,0,0);

	number energy =  GCInteraction<number>::_repulsive_lj(*r, force, update_forces);

	if(update_forces)
	{
		p->force -= force;
		q->force += force;
	}

	return energy;
}
*/

template<typename number>
number GCInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {



	number rnorm = SQR(r->x) + SQR(r->y) + SQR(r->z);
	number energy = 0.5 * _k * rnorm;

	if (update_forces) {
		LR_vector<number> force(*r) ;
		force *= (-1.0f * _k);

		p->force -= force;
		q->force += force;
	}

	return energy;
}



// TODO: Delete This Bottom Line? It's commented out anyway??

#endif /* LJINTERACTION_H_ */
