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
class CUDADNANMINTERACTION: public CUDABaseInteraction<number, number4>, public DNANMInteraction<number> {
public:
    char* _d_spring_pottype, _h_spring_pottype;
    number* _d_spring_potential, _h_spring_potential;
    number* _d_spring_eqdist, _h_spring_eqdist;
    size_t _spring_param_size;
    int offset; //Basically will only come into play if proteins are after dna in topology file. Just adjust proteins index for the spring parameter arrays
	CUDADNANMINTERACTION();
	virtual ~CUDADNANMINTERACTION();

	virtual void get_settings(input_file &inp);
	virtual void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDADNANMINTERACTION_H_ */
