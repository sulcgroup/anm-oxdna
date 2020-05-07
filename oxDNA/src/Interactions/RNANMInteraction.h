/**
 * @brief Handles interactions between RNA nucleotides and Proteins.
 * Implements RNA2 model and ANM Protein Model
 *
 * Jonah (May 2020)
 *
 */

#ifndef RNANM_INTERACTION_H
#define RNANM_INTERACTION_H

#include "RNANMInteraction.h"
#include "BaseInteraction.h"
#include "RNAInteraction2.h"

template<typename number>
class RNANMInteraction: public RNA2Interaction<number> {

protected:
	//parameters of the interaction

	int nrna, npro, nrnas, _firststrand;

    map<pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
    map<pair<int, int>, pair<char, double> > _potential;
    number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut;
    number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;

    virtual void allocate_particles(BaseParticle<number> **particles, int N);
    virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

public:
	enum {
            SPRING = 8,
            PRO_EXC_VOL = 9,
            PRO_RNA_EXC_VOL = 10
	};

    char _parameterfile[500];

	RNANMInteraction();    // Constructor
	virtual ~RNANMInteraction() {} // Destructor

	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
    virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual void get_settings(input_file &inp); //get settings from input file
	virtual void init(); // initialisation


    number _protein_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces);
    number _protein_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
    number _protein_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

    number _protein_rna_exc_volume(BaseParticle<number> *p,BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
    number _protein_rna_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
};

#endif /* RNANM_INTERACTION_H_ */
