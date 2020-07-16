/*
 * DNANMInteraction_relax.cpp
 *
 *  Created on: June 4, 2020
 *      Author: Jonah
 */

#include <fstream>

#include "DNANMInteraction_relax.h"
#include "../Particles/DNANucleotide.h"
#include "../Particles/ACParticle.h"

template<typename number>
DNANMInteraction_relax<number>::DNANMInteraction_relax() : DNANMInteraction<number>() {
	OX_LOG(Logger::LOG_INFO, "Using non-diverging backbone potential (DNA_relax2 interaction)");
	_fmax = 10.0f;
}

template<typename number>
DNANMInteraction_relax<number>::~DNANMInteraction_relax() {

}

template<typename number>
number DNANMInteraction_relax<number>::_dna_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!this->_check_bonded_neighbour(&p, &q, r)) {
	    return (number) 0.f;
	}

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - FENE_R0_OXDNA;
	number energy = 0.;

	// compute the x for which the potential becomes a straight line
	number xmax = (-FENE_EPS + sqrt(FENE_EPS * FENE_EPS + 4.f * _fmax * _fmax * FENE_DELTA2)) / (2.f * _fmax);
	number fenemax = - (FENE_EPS / 2.f) * log (1.f - SQR(xmax) / FENE_DELTA2);
	//OX_LOG(Logger::LOG_INFO, "relax: xmax = %g ", xmax);

	LR_vector<number> force;
	
	if (fabs(rbackr0) < xmax) {
		// we use the standard FENE
		energy = - (FENE_EPS / 2.f) * log(1.f - SQR(rbackr0) / FENE_DELTA2);

		if (update_forces) force = rback * (-(FENE_EPS * rbackr0  / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
	}
	else {
		// we use the straight potential
		energy = fenemax + _fmax * (fabs(rbackr0) - xmax);  

		if (update_forces) force = rback * (- _fmax * copysign(1., rbackr0) / rbackmod);
	}

	if (update_forces) {
		p->force -= force;
		q->force += force;

		p->torque -= p->orientationT * p->int_centers[DNANucleotide<number>::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[DNANucleotide<number>::BACK].cross(force);
	}

	return energy;
}

template<typename number>
void DNANMInteraction_relax<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
	}
}

template<typename number>
void DNANMInteraction_relax<number>::get_settings(input_file &inp) {
	DNANMInteraction<number>::get_settings(inp);

	getInputNumber (&inp, "relax_fmax", &_fmax, 0);

	OX_LOG(Logger::LOG_INFO, "relax_fmax = %g ", _fmax);
}

template<typename number>
number DNANMInteraction_relax<number>::_protein_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
    if(p->btype >= 0 || q->btype >= 0)    //this function is only for proteins
    {
        return 0.f;
    }
    number eqdist;
    pair <int, int> keys (std::min(p->index, q->index), std::max(p->index, q->index));
    eqdist = this->_rknot[keys];
    if (eqdist != 0.0) { //only returns number if eqdist is in .par file
            //Harmonic Spring Potential
            if ((eqdist < 0.0) || (eqdist > 3.0))  //ensures r0 is less than 7 Angstrom cutoff and nonnegative
            {
                if (keys.first + 1 != keys.second) {
                    throw oxDNAException("No rknot or invalid rknot value for particle %d and %d rknot was %f",
                                         q->index, p->index, eqdist);
                }
            }
            number _k = 20.f; //stiffness of the spring

            number rnorm = r->norm();
            number rinsta = sqrt(rnorm);
            number energy = 0.5 * _k * SQR((rinsta - eqdist));

            if (update_forces) {
                LR_vector<number> force(*r);
                force *= (-1.0f * _k) * (rinsta - eqdist) / rinsta;
                p->force -= force;
                q->force += force;
            }
            return energy;

    } else {
        return (number) 0.f; //returns 0 if no rknot value in parameter value aka they aren't bonded
    }
}





template class DNANMInteraction_relax<float>;
template class DNANMInteraction_relax<double>;
