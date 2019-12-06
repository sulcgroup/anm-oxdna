/*
 * CUDADNAInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDADNANMInteraction.h"
#include "CUDADNAInteraction.h"

#include "CUDA_DNA.cuh"
#include "CUDA_ANM.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/DNA2Interaction.h"
#include "../../Interactions/DNANMInteraction.h"

template<typename number, typename number4>
CUDADNANMInteraction<number, number4>::CUDADNANMInteraction() {

}

template<typename number, typename number4>
CUDADNANMInteraction<number, number4>::~CUDADNANMInteraction() {

}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::get_settings(input_file &inp) {
	CUDADNAInteraction<number, number4>::get_settings(inp);
    // get_settings for DNANM
    char parameterfile[500];
    getInputString(&inp, "PARFILE", parameterfile, 0);

    //Addition of Reading Parameter File
    int key1, key2;
    char potswitch;
    double potential, dist;
    string carbons;
    fstream parameters;
    parameters.open(parameterfile, ios::in);
    getline (parameters,carbons);
    int cc = stoi(carbons);
    pair<char, number>* _spring_params = NULL;
    _spring_params = new pair<char, number>[cc][cc];
    pair<char, number> d_fault ('x', 0.f);
    for(int i = 0; i<cc; i++) {
        for (int j = 0; j < cc; j++) {
            if(i > j) break;
            *_spring_params[i][j] = d_fault;
        }
    }
    if (parameters.is_open())
    {
        while (parameters.good())
        {
            parameters >> key1 >> key2 >> dist >> potswitch >> potential;
            pair <char, double> pot (potswitch, potential);
            _spring_params[key1][key2] = pot;
        }
        parameters.close();
    }
    else
    {
        throw oxDNAException("ParameterFile Could Not Be Opened");
    }
}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::cuda_init(number box_side, int N) {
    CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
    DNANMInteraction<number>::init();
    char* _spring_pottype = new char[this->npro*cc];
    double* _spring_potential = new double[cc*cc];
    double* _spring_eqdist = new double[cc*cc];
    for(std::map<pair<int, int>, double>::iterator it=DNANMInteraction<number>::_rknot.begin(); it!=DNANMInteraction<number>::_rknot.end(); ++it) {
        map<pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
        map<pair<int, int>, pair<char, double> > _potential;
    }


    float f_copy = this->_hb_multiplier;
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

    number tmp[50];
    for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_EPS[i][j][k];

    COPY_ARRAY_TO_CONSTANT(MD_F1_EPS, tmp, 50);

    for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_SHIFT[i][j][k];

    COPY_ARRAY_TO_CONSTANT(MD_F1_SHIFT, tmp, 50);

    COPY_ARRAY_TO_CONSTANT(MD_F1_A, this->F1_A, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RC, this->F1_RC, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_R0, this->F1_R0, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_BLOW, this->F1_BLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_BHIGH, this->F1_BHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RLOW, this->F1_RLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RHIGH, this->F1_RHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RCLOW, this->F1_RCLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RCHIGH, this->F1_RCHIGH, 2);

    COPY_ARRAY_TO_CONSTANT(MD_F2_K, this->F2_K, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RC, this->F2_RC, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_R0, this->F2_R0, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_BLOW, this->F2_BLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_BHIGH, this->F2_BHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RLOW, this->F2_RLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RHIGH, this->F2_RHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RCLOW, this->F2_RCLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RCHIGH, this->F2_RCHIGH, 2);

    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_A, this->F5_PHI_A, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_B, this->F5_PHI_B, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XC, this->F5_PHI_XC, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XS, this->F5_PHI_XS, 4);

    if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
    if (_use_debye_huckel){
        // copied from DNA2Interaction::init() (CPU), the least bad way of doing things
        // We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
        number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
        // RHIGH gives the distance at which the smoothing begins
        _debye_huckel_RHIGH = 3.0 * lambda;
        _minus_kappa = -1.0/lambda;

        // these are just for convenience for the smoothing parameter computation
        number x = _debye_huckel_RHIGH;
        number q = _debye_huckel_prefactor;
        number l = lambda;

        // compute the some smoothing parameters
        _debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(-4.*x*x*x * l * l * q );
        _debye_huckel_RC = x*(q*x + 3. * q* l )/(q * (x+l));

        number debyecut;
        if (this->_grooving){
            debyecut = 2.0f * sqrt((POS_MM_BACK1)*(POS_MM_BACK1) + (POS_MM_BACK2)*(POS_MM_BACK2)) + _debye_huckel_RC;
        }
        else{
            debyecut =  2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
        }
        // the cutoff radius for the potential should be the larger of rcut and debyecut
        if (debyecut > this->_rcut){
            this->_rcut = debyecut;
            this->_sqr_rcut = debyecut*debyecut;
        }
        // End copy from DNA2Interaction

        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );

        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_sigma, this->_pro_sigma, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rstar, this->_pro_rstar, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rc, this->_pro_rcut, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_b, this->_pro_b, sizeof(float)) );

        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_sigma, this->_pro_backbone_sigma, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rstar, this->_pro_backbone_rstar, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rc, this->_pro_backbone_rcut, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_b, this->_pro_backbone_b, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_stiffness, this->_pro_backbone_stiffness, sizeof(float)) );

        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_sigma, this->_pro_base_sigma, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rstar, this->_pro_base_rstar, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rc, this->_pro_base_rcut, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_b, this->_pro_base_b, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_stiffness, this->_pro_base_stiffness, sizeof(float)) );
    }
}
}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				dna_forces_edge_nonbonded<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, d_box);

				this->_sum_edge_forces_torques(d_forces, d_torques);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

				dna_forces_edge_bonded<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf);
			}
			else {
				dna_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box);
				CUT_CHECK_ERROR("forces_second_step simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		dna_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations,  d_forces, d_torques, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box);
		CUT_CHECK_ERROR("forces_second_step no_lists error");
	}
}

template class CUDADNAInteraction<float, float4>;
template class CUDADNAInteraction<double, LR_double4>;
