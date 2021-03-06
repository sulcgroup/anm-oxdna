/*
 * CUDADNANMInteraction.h
 *
 *  Created on: 6/11/19
 *      Author: jonah
 */

#ifndef CUDADNANMINTERACTION_H_
#define CUDADNANMINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNANMInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model with ANM protein model as implemented in DNANMInteraction.
 */
template<typename number, typename number4>
class CUDADNANMInteraction: public CUDABaseInteraction<number, number4>, public DNANMInteraction<number> {
public:

    bool _use_debye_huckel;
    bool _use_oxDNA2_coaxial_stacking;
    bool _use_oxDNA2_FENE;
    // copied from DNA2Interaction.h (CPU) (and change number -> float), the least bad way of doing things
    float _salt_concentration;
    bool _debye_huckel_half_charged_ends;
    float _debye_huckel_prefactor;
    float _debye_huckel_lambdafactor;

    float _debye_huckel_RC; // this is the maximum interaction distance between backbones to interact with DH
    float _debye_huckel_RHIGH; // distance after which the potential is replaced by a quadratic cut-off
    float _debye_huckel_B; // prefactor of the quadratic cut-off
    float _minus_kappa;

    float _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut;
    float _pro_base_sigma, _pro_base_rstar, _pro_base_b, _pro_base_rcut;
    float _pro_sigma, _pro_rstar, _pro_b, _pro_rcut;

    number *_spring_eqdist, *_spring_potential; //Temp arrays for parameter storage

    //compressed parameter arrays
    number *_h_aff_gamma, *_d_aff_gamma;
    number *_h_aff_eqdist, *_d_aff_eqdist;
    int *_h_affected, *_d_affected;
    int *_affected_len, *_h_affected_indx, *_d_affected_indx;
    size_t _spring_param_size_number;

    bool _read_par;

    int offset; //Will only come into play if proteins are after dna in topology file (particle id wise). Adjusts proteins index for the spring parameter arrays
	CUDADNANMInteraction();
	virtual ~CUDADNANMInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDADNANMINTERACTION_H_ */
