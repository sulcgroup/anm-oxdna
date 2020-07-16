/*
 * ACInteraction.h
 *
 *  Created on: 13/Nov/18
 *      Author: Jonah
 */

#ifndef ACTINTERACTION_H_
#define ACTINTERACTION_H_

#include "BaseInteraction.h"
#include "../Particles/ACTParticle.h"
/**
 * @brief Handles (generalised) Anisotropic-Chain interactions between spheres of size .1573 simulation units
 *
 *
 * This interaction is selected with
 * interaction_type = AC
 *
 * Input options:
 *
 * @verbatim
 * none
 */
template <typename number>
class ACTInteraction: public ACTInteraction<number> {
protected:

//	map<pair<int, int>, double> _rknot; //eqdist of each bond of psuedobonds
    map<int, tuple <double, double, double> > _ang_vals; //eqdist of each bond of psuedobonds
//	map<pair<int, int>, pair<char, double> > _potential; //switch to tell lj, FENE or spring as well as strength for each pair of particles
//
	number _sigma_bond_sqr, _sigma_tor_sqr;
//
//	inline number _exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
//	inline number _spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces );
//	inline number _repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);
//    inline number _repulsive_lj_quart(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);

public:
	enum {
//		SPRING_POTENTIAL = 0,
//		EXC_VOL = 1
        ANGULAR_POTENTIAL = 2

	};

	ACTInteraction();
	virtual ~ACTInteraction();

	virtual void get_settings(input_file &inp); // Need to add parameters for angular potential
	virtual void init();  //will add sigma 1 and 2 to this

	virtual void allocate_particles(BaseParticle<number> **particles, int N); //change to ACT particles
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles); //change to ACT particles

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	//above 3 will need the most work

	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);  //Still don't think this does much
};
//
////The Below function is currently unused
//template<typename number>
//number ACInteraction<number>::_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces) {
//
//	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
//	number energy = (number) 0;
//	if(rnorm < SQR(_rc)) {
//		if(rnorm > SQR(_rstar)) {
//			number rmod = sqrt(rnorm);
//			number rrc = rmod - _rc;
//			energy = EXCL_EPS * _b * SQR(rrc);
//			if(update_forces) force = -r * (2 * EXCL_EPS * _b * rrc/ rmod);
//		}
//		else {
//			number tmp = SQR(_sigma) / rnorm;
//			number lj_part = tmp * tmp * tmp;
//			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
//			if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
//		}
//	}
//
//	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;
//
//	return energy;
//}


//template<typename number>
//number ACInteraction<number>::_repulsive_lj_quart(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces) {
//
//    number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
//    number energy = (number) 0;
//    if(rnorm < SQR(_rc)) {
//        if(rnorm > SQR(_rstar)) {
//            number rmod = sqrt(rnorm);
//            number rrc = rmod - _rc;
//            energy = EXCL_EPS * _b * SQR(SQR(rrc));
//            if(update_forces) force = -r * (EXCL_EPS * 4.f * _b * CUB(rrc)/ rmod);
//        }
//        else {
//            number tmp = SQR(_sigma) / rnorm;
//            number lj_part = tmp * tmp * tmp;
//            energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
//            if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
//        }
//    }
//
//    if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;
//
//    return energy;
//}


//template<typename number>
//number ACInteraction<number>::_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//	if (p->index != q->index && abs(q->index - p->index) != 1){
//		LR_vector<number> force(0,0,0);
//
//		number energy =  this->_repulsive_lj_quart(*r, force, update_forces);
//
//		if(update_forces)
//		{
//			p->force -= force;
//			q->force += force;
//		}
//
//		return energy;
//	} else {
//		return (number) 0.f;
//	}
//}



//template<typename number>
//number ACInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//	number eqdist;
//	char interactiontype;
//    pair <int,int> keys (std::min(p->index, q->index), std::max(p->index, q->index));
//	if (p->index != q->index)
//	{
//		eqdist = _rknot[keys];
//		interactiontype = _potential[keys].first;
//		if (eqdist != 0.0){ //only returns number if eqdist is in .par file
//			switch (interactiontype){
//				case 's':
//					{
//						//Harmonic Spring Potential
//						if ((eqdist < 0.0) || (eqdist > 3.0))  //ensures r0 is less than 7 Angstrom cutoff and nonnegative
//						{
//							if (keys.first+1 != keys.second){
//								throw oxDNAException("No rknot or invalid rknot value for particle %d and %d rknot was %f", q->index, p->index, eqdist);
//							}
//						}
//						number _k = _potential[keys].second; //stiffness of the spring
//						if ((_k == 0) || (_k < 0)){
//							throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
//						}
//						number rnorm = r->norm();
//						number rinsta = sqrt(rnorm);
//						number energy = 0.5 * _k * SQR(rinsta-eqdist);
//
//						if (update_forces)
//						{
//							LR_vector<number> force(*r );
//							force *= (-1.0f * _k ) * (rinsta-eqdist)/rinsta;
////							printf("p %d, q%d | f.x %.5f f.y %.5f f.z %.5f\n", p->index, q->index, force.x, force.y, force.z);
//							p->force -= force;
//							q->force += force;
//						//printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
//						//printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
//						}
//						return energy;
//					} break;
//                case 'i':
//                    {
//                        //Every Possible Pair of Particles Needs to be Calculated
//                        number _k = _potential[keys].second; //stiffness of the spring
//                        if ((_k == 0) || (_k < 0)){
//                            throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
//                        }
//                        number rnorm = r->norm();
//                        number rinsta = sqrt(rnorm);
//                        number energy = 0.5 * _k * SQR(rinsta-eqdist)*(1/eqdist);
//
//                        if (update_forces)
//                        {
//                            LR_vector<number> force(*r );
//                            force *= (-1.0f * _k ) * ((rinsta-eqdist)/rinsta) * (1/eqdist);
//                            p->force -= force;
//                            q->force += force;
//                        //printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
//                        //printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
//                        }
//                        return energy;
//                    } break;
//                case 'e':
//                    {
//                        //Every Possible Pair of Particles Needs to be Calculated
//                        number _k = _potential[keys].second; //stiffness of the spring
//                        if ((_k == 0) || (_k < 0)){
//                            throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
//                        }
//                        number rnorm = r->norm();
//                        number rinsta = sqrt(rnorm);
//                        number energy = 0.5 * _k * SQR(rinsta-eqdist)*(1/(eqdist*eqdist));
//
//                        if (update_forces)
//                        {
//                            LR_vector<number> force(*r );
//                            force *= (-1.0f * _k ) * ((rinsta-eqdist)/rinsta) * (1/(eqdist*eqdist));
//                            p->force -= force;
//                            q->force += force;
//                        //printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
//                        //printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
//                        }
//                        return energy;
//                    } break;
//                default:
//                    {
//                        throw oxDNAException("Interaction type specified in .par file is Invalid, particles %d and %d, switch %c",p->index,q->index,interactiontype);
//                    }
//				}
//		} else {
//			return (number) 0.f; //returns 0 if no rknot value in parameter value aka they aren't bonded
//		}
//	} else {
//		return (number) 0.f; //returns 0 if particle pair consists of particle and itself
//	}
//
//}


#endif /* ACINTERACTION_H_ */
