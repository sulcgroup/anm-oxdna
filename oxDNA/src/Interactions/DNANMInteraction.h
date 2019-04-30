/*
 * DNANMInteraction.h
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 */

#ifndef DNANM_INTERACTION_H
#define DNANM_INTERACTION_H

#include "DNANMInteraction.h"
#include "DNA2Interaction.h"

template<typename number>
class DNANMInteraction: public DNA2Interaction<number> {

protected:
	int ndna; //How many particles of DNA type: Used in allocate_particles
	int npro; //How many particles of AC type: Used in allocate_particles


public:
	DNANMInteraction();
	virtual ~DNANMInteraction();
	virtual void get_settings(input_file &inp);
	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	number _exlcudedvolume(BaseParticle<number> *p,BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual void DNANMInteraction<number>::init();
};


template class DNANMInteraction<float>;
template class DNANMInteraction<double>;
#endif /* SRC_INTERACTIONS_DNANMINTERACTION_H_ */
