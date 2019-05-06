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
#include "ACInteraction.h"

template<typename number>
class DNANMInteraction: public DNA2Interaction<number>, public ACInteraction<number> {

protected:
	int ndna; //How many particles of DNA type: Used in allocate_particles
	int npro; //How many particles of AC type: Used in allocate_particles
	number _pro_dna_sigma, _pro_dna_rstar, _pro_dna_b, _pro_dna_rcut, _pro_dna_stiffness;
	map<pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
	map<pair<int, int>, pair<char, double> > _potential;


public:
	DNANMInteraction();
	virtual ~DNANMInteraction();
	virtual void get_settings(input_file &inp);
	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	number _protein_dna_exc_volume(BaseParticle<number> *p,BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _protein_dna_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
	virtual void init();
};


template class DNANMInteraction<float>;
template class DNANMInteraction<double>;
#endif /* SRC_INTERACTIONS_DNANMINTERACTION_H_ */
