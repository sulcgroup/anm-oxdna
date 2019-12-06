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
    char* _spring_pottype;
    double* _spring_potential;
    double* _spring_eqdist;
	CUDADNANMINTERACTION();
	virtual ~CUDADNANMINTERACTION();

	virtual void get_settings(input_file &inp);
	virtual void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDADNANMINTERACTION_H_ */
