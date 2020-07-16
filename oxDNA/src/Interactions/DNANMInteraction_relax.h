/*
 * DNANMInteraction_relax.h
 *
 *  Created on: June 4th, 2020
 *      Author: Jonah
 */

#ifndef DNANM_INTERACTION_RELAX_H
#define DNANM_INTERACTION_RELAX_H

#include "BaseInteraction.h"
#include "DNANMInteraction.h"

// HARMONIC BONDED BACKBONE-BACKBONE
#define HARMONIC_R0 FENE_R0_OXDNA

/**
 * @brief
 * All forces in DNANM kept same with the exception of the dna-backbone
 *
 * Modified version of DNAInteraction which modifies the bonded backbone-backbone potential so it does not diverge
 *
 * Replaces the bonded backbone-backbone FENE potential with a harmonic potential. This is to allow very stressed initial
 * configurations, which might otherwise cause the simulation to fail, to relax to a sensible structure
 *
 * This interaction takes 3 compulsory arguments:
 *
 * This interaction is selected with
 * interaction_type = DNANM_relax
 *
@verbatim
relax_type = <string> (Possible values: constant_force, harmonic_force; Relaxation algorithm used)
relax_strength = <float> (Force constant for the replacement of the FENE potential)
@endverbatim
 */
template <typename number>
class DNANMInteraction_relax : public DNANMInteraction<number> {
protected:

	virtual number _dna_backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _protein_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _fmax;

public:
	DNANMInteraction_relax();
	virtual ~DNANMInteraction_relax();

	void check_input_sanity(BaseParticle<number> **particles, int N);
	void get_settings(input_file &inp);
};

#endif /*DNA_INTERACTION_RELAX2_H*/
