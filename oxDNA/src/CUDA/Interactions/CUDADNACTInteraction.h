/*
 * CUDADNACTInteraction.h
 *
 *  Created on: 6/11/19
 *      Author: jonah
 */

#ifndef CUDADNACTINTERACTION_H_
#define CUDADNACTINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNACTInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model with ANM protein model as implemented in DNANMInteraction.
 */
template<typename number, typename number4>
class CUDADNACTInteraction: public CUDABaseInteraction<number, number4>, public DNACTInteraction<number> {
public:
    //COPIED from CUDADNAInteraction
//    enum {
//        DEBYE_HUCKEL = 7,
//        SPRING = 8,
//        PRO_EXC_VOL = 9,
//        PRO_DNA_EXC_VOL = 10
//    };
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
    float _ktor, _kbend;

    char *_d_spring_pottype, *_h_spring_pottype;
    number *_d_spring_potential, *_h_spring_potential;
    number *_d_spring_eqdist, *_h_spring_eqdist;
    number *_d_ang_params, *_h_ang_params;
    size_t _spring_param_size_number, _spring_param_size_char, _ang_param_size;
    int offset; //Will only come into play if proteins are after dna in topology file (particle id wise). Adjusts proteins index for the spring parameter arrays
	CUDADNACTInteraction();
	virtual ~CUDADNACTInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDADNACTINTERACTION_H_ */
