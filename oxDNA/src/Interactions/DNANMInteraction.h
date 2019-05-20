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
#include "BaseInteraction.h"

template<typename number>
class DNANMInteraction: public DNA2Interaction<number> {

protected:
	int ndna;//How many particles of DNA type: Used in allocate_particles
	int npairs;//How many bonds b/t different particle types
	int ndnas;//Number of Strands that are dna
	number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
	map<pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
	map<pair<int, int>, pair<char, double> > _potential;
	number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut;
	std::vector< std::pair<int, int> > _ippairs; //interparticle pairs



public:
	enum {
			SPRING = 8,
			EXC_VOL = 9
			//Assigned 8 and 9 so it won't overwrite the already existing DNA function pointers in the _int_map
	};

	DNANMInteraction();
	virtual ~DNANMInteraction();
	virtual void get_settings(input_file &inp);
	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	number _protein_dna_exc_volume(BaseParticle<number> *p,BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _protein_dna_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
	virtual void init();
	number _protein_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);
	number _protein_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _protein_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
};


template class DNANMInteraction<float>;
template class DNANMInteraction<double>;
#endif /* SRC_INTERACTIONS_DNANMINTERACTION_H_ */
