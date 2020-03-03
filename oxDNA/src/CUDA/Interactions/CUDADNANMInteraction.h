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
    //COPIED from CUDADNAInteraction
    enum {
        DEBYE_HUCKEL = 7,
        SPRING = 8,
        PRO_EXC_VOL = 9,
        PRO_DNA_EXC_VOL = 10
    };
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

    float _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    float _pro_base_sigma, _pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
    float _pro_sigma, _pro_rstar, _pro_b, _pro_rcut;

    char *_d_spring_pottype, *_h_spring_pottype;
    number *_d_spring_potential, *_h_spring_potential;
    number *_d_spring_eqdist, *_h_spring_eqdist;
    size_t _spring_param_size_double, _spring_param_size_char;
    int offset; //Will only come into play if proteins are after dna in topology file (particle id wise). Adjusts proteins index for the spring parameter arrays
	CUDADNANMInteraction();
	virtual ~CUDADNANMInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDADNANMINTERACTION_H_ */
