/*
 * GCInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef GC_INTERACTION_H_
#define GC_INTERACTION_H_

#include "BaseInteraction.h"
#include "../Particles/GCParticle.h"


/**
 * @brief Handles (generalised) Lennard-Jones interactions between spheres of size 1 or a Kob-Andersen interaction.
 *
 * TODO: the Kob-Andersen mixture should be better implemented.
 *
 * This interaction is selected with
 * interaction_type = LJ
 *
 * Input options:
 *
 * @verbatim
LJ_rcut = <float> (interaction cutoff)
[LJ_kob_andersen = <bool> (Simulate a Kob-Andersen mixture. Defaults to false.)]
[LJ_n = <int> (Generalised LJ exponent. Defaults to 6, which is the classic LJ value.)]
@endverbatim
 */
template <typename number>
class GCInteraction: public BaseInteraction<number, GCInteraction<number> > {
protected:

	/*
	number _E_cut[3];
	bool _is_ka_mixture;
	int _N_A, _N_B;
	int _n[3];
	number _sigma[3];
	number _sqr_sigma[3];
	number _epsilon[3];
	number _sqr_LJ_rcut[3];
    */

	number _k; //stiffness of the spring
	number _r; //radius of the amino acid

	number _sigma, _rstar, _b, _rc, STRENGTH; //TODO initialize in init() !!!!

	inline number _exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	inline number _spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces );
	inline number _repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces);

public:
	enum {
		SPRING_POTENTIAL = 0,
		EXC_VOL = 1

	};

	GCInteraction();
	virtual ~GCInteraction();

	virtual void get_settings(input_file &inp);
	/* TODO: Figure out these values */
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


template<typename number>
number GCInteraction<number>::_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
			if(update_forces) force = -r * (2 * EXCL_EPS * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}


template<typename number>
number GCInteraction<number>::_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	LR_vector<number> force(0,0,0);

	number energy =  GCInteraction<number>::_repulsive_lj(*r, force, this->_sigma, this->_rstar, this->_b, this->_rc,update_forces);

	if(update_forces)
	{
		p->force -= force;
		q->force += force;
	}

	return energy;
}


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
