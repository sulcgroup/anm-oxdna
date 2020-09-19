/*
 * CUDADNAInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDADNANMInteraction.h"
#include "CUDA_DNA.cuh"
#include "CUDA_ANM.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/DNA2Interaction.h"
#include "../../Interactions/DNANMInteraction.h"

template<typename number, typename number4>
CUDADNANMInteraction<number, number4>::CUDADNANMInteraction() {

    _read_par = false;
    //Not copied over to device memory
    _spring_potential = NULL;
    _spring_eqdist = NULL;
    _affected_len = NULL;

    //Copied over to device memory
    _h_affected_indx = NULL;
    _d_affected_indx = NULL;

    _h_affected = NULL;
    _d_affected = NULL;

    _h_aff_eqdist = NULL;
    _d_aff_eqdist = NULL;

    _h_aff_gamma = NULL;
    _d_aff_gamma = NULL;

    _spring_param_size_number = 0;
}

template<typename number, typename number4>
CUDADNANMInteraction<number, number4>::~CUDADNANMInteraction() {
    //Delete All pointers required for spring potential parameters
    //Delete All pointers required for spring potential parameters
    if(_spring_potential != NULL) delete[] _spring_potential;
    if(_spring_eqdist != NULL) delete[] _spring_eqdist;

    if(_affected_len != NULL) delete[] _affected_len;

    if(_d_affected != NULL) CUDA_SAFE_CALL(cudaFree(_d_affected));
    if(_d_aff_gamma != NULL) CUDA_SAFE_CALL(cudaFree(_d_aff_gamma));
    if(_d_aff_eqdist != NULL) CUDA_SAFE_CALL(cudaFree(_d_aff_eqdist));
    if(_d_affected_indx != NULL) CUDA_SAFE_CALL( cudaFree(_d_affected_indx) );

    if(_h_affected != NULL) delete[] _h_affected;
    if(_h_aff_gamma != NULL) delete[] _h_aff_gamma;
    if(_h_aff_eqdist != NULL) delete[] _h_aff_eqdist;
    if(_h_affected_indx != NULL) delete[] _h_affected_indx;
}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::get_settings(input_file &inp) {
    _use_debye_huckel = false;
    _use_oxDNA2_coaxial_stacking = false;
    _use_oxDNA2_FENE = false;
    std::string inter_type;
    if (!getInputString(&inp, "parfile", this->_parameterfile, 0) == KEY_FOUND){
        throw oxDNAException("Key 'PARFILE' not found. Necessary for Protein sims.");
    }

    char s[5] = "none";
    if(strcmp(this->_parameterfile, s) != 0) _read_par = true;

    if (!getInputString(&inp, "topology", this->_topology_filename, 0) == KEY_FOUND){
        throw oxDNAException("Key 'topology_file' not found.");
    }
    if (getInputString(&inp, "interaction_type", inter_type, 0) == KEY_FOUND){
        if (inter_type.compare("DNANM") == 0) {
            _use_debye_huckel = true;
            _use_oxDNA2_coaxial_stacking = true;
            _use_oxDNA2_FENE = true;

            // we don't need the F4_... terms as the macros are used in the CUDA_DNA.cuh file; this doesn't apply for the F2_K term
            this->F2_K[1] = CXST_K_OXDNA2;
            _debye_huckel_half_charged_ends = true;
            this->_grooving = true;
            // end copy from DNA2Interaction

            // copied from DNA2Interaction::get_settings() (CPU), the least bad way of doing things
            getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
            getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

            // lambda-factor (the dh length at T = 300K, I = 1.0)
            _debye_huckel_lambdafactor = 0.3616455f;
            getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

            // the prefactor to the Debye-Huckel term
            _debye_huckel_prefactor = 0.0543f;
            getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
            // End copy from DNA2Interaction

        }
    }
    this->DNANMInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::cuda_init(number box_side, int N) {
    CUDABaseInteraction<number, number4>::cuda_init(box_side, N);

//    Addition of Reading Parameter File -> Moved from get_settings due to needing to fill variables that are filled in the CPU version of DNANMInteraction::read_topology
    fstream top;
    int tmp1, tmp2;
    top.open(this->_topology_filename, ios::in);
    if (top.is_open()){
        top >> tmp1 >> tmp2 >> this->ndna >> this->npro >> this->ndnas;
        top >> this->_firststrand;
        top.close();
    } else {
        throw oxDNAException("Could not open Topology File");
    }


    if(this->_firststrand < 0) offset = 0;
    else if(this->_firststrand > 0) offset = this->ndna;
    else throw oxDNAException("No Strand should have an ID of 0");

    if(_read_par){
        //Initalizing Some Host and Device Arrays for Spring Parameters
        _spring_param_size_number = sizeof(number) * (this->npro*this->npro);
        _spring_potential = new number[this->npro*this->npro]();
        _spring_eqdist = new number[this->npro*this->npro]();

        char potswitch = 'x';
        number potential = 0.f, dist = 0.f;
        for(int i = 0; i< (this->npro*this->npro); i++){
            _spring_eqdist[i] = dist;
            _spring_potential[i] = potential;
        }

        auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
            if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
            if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
            if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
            if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
            if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
        };

        int key1, key2 = 0;
        string carbons;
        fstream parameters;
        parameters.open(this->_parameterfile, ios::in);
        getline (parameters, carbons);

        //total connections
        int spring_connection_num = 0;

        //allocate and declare affected_len vector
        _affected_len = new int[this->npro]();
        for(int i = 0; i < this->npro; i++) _affected_len[i] = 0;

        if (parameters.is_open())
        {
            while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
            {
                valid_spring_params(N, key1, key2, dist, potswitch, potential);
                spring_connection_num += 1;

                if(offset != 0) {
                    key1 -= offset;
                    key2 -= offset;
                }

                _affected_len[key1] += 1;
                _affected_len[key2] += 1;

                _spring_potential[key1*this->npro + key2] = potential;
                _spring_eqdist[key1*this->npro + key2] = dist;

                _spring_potential[key2*this->npro + key1] = potential;
                _spring_eqdist[key2*this->npro + key1] = dist;
            }
            parameters.close();
        }
        else
        {
            throw oxDNAException("ParameterFile Could Not Be Opened");
        }

        //Compressed Parameter Initialization
        _h_affected_indx = new int[this->npro + 1]();
        _h_affected = new int[spring_connection_num*2]();
        _h_aff_gamma = new number[spring_connection_num*2]();
        _h_aff_eqdist = new number[spring_connection_num*2]();
        number zero = (number) 0.f;
        for(int i = 0; i < this->npro+1; i++) _h_affected_indx[i] = 0;
        for(int i = 0; i < spring_connection_num*2; i++){
            _h_affected[i] = 0;
            _h_aff_gamma[i] = zero;
            _h_aff_eqdist[i] = zero;
        }

        //Compressed Index
        int param_indx = 0;
        //For each residue
        for(int i = 0; i < this->npro; i++){
            //Fill _h_affected filtering through larger arrays filled in parameter file reading
            for(int j = i*this->npro; j < i*this->npro+this->npro; j++){
                if(_spring_eqdist[j] != 0.f){
                    //Affected List, Access is controlled with indices in _h_affected_indx
                    _h_affected[param_indx] = j % this->npro;
                    //Stored in same way for easy access, spring constants
                    _h_aff_gamma[param_indx] = _spring_potential[j];
                    //eq_distance
                    _h_aff_eqdist[param_indx] = _spring_eqdist[j];
                    param_indx += 1;
                }
            }
        }

        //Don't need Larger arrays anymore, safe to delete
        if(_spring_eqdist != NULL) delete[] _spring_eqdist;
        _spring_eqdist = NULL; //Otherwise dangling Pointer
        if(_spring_potential != NULL) delete[] _spring_potential;
        _spring_potential = NULL;

        //Allocation and Copying of Compressed Parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_affected, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_affected, _h_affected, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_gamma, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_gamma, _h_aff_gamma, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_eqdist, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_eqdist, _h_aff_eqdist, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        int ind = 0;
        _h_affected_indx[0] = 0;
        //make indx access list where: _h_affected_indx[i] lower bound of i's parameters, _h_affected_indx[i+1] upper bound of i's parameters
        for(int i = 0; i < this->npro; i++){
            ind += _affected_len[i];
            _h_affected_indx[i+1] += ind;
        }

        //Don't need this anymore
        if(_affected_len != NULL) delete[] _affected_len;
        _affected_len = NULL;

        //Allocation and copying of Indice List for accessing compressed parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_affected_indx, (this->npro+1)*sizeof(int)));
        CUDA_SAFE_CALL( cudaMemcpy(_d_affected_indx, _h_affected_indx, (this->npro+1)*sizeof(int), cudaMemcpyHostToDevice));

        //Memory Used by Parameters
        float param_memory_mb = (spring_connection_num * 2 * sizeof(int) + 2 * spring_connection_num * 2 * sizeof(number)
                                 + (this->npro + 1) * sizeof(int) + 4 * this->npro * sizeof(number))/SQR(1024);
        OX_LOG(Logger::LOG_INFO, "Spring Parameters Size: %.2f MB", param_memory_mb);

    } else OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");

    // Copied from CUDADNAINTERACTION
    DNAInteraction<number>::init();

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
    if (_use_debye_huckel) {
        // copied from DNA2Interaction::init() (CPU), the least bad way of doing things
        // We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
        number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
        // RHIGH gives the distance at which the smoothing begins
        _debye_huckel_RHIGH = 3.0 * lambda;
        _minus_kappa = -1.0 / lambda;

        // these are just for convenience for the smoothing parameter computation
        number x = _debye_huckel_RHIGH;
        number q = _debye_huckel_prefactor;
        number l = lambda;

        // compute the some smoothing parameters
        _debye_huckel_B = -(exp(-x / l) * q * q * (x + l) * (x + l)) / (-4. * x * x * x * l * l * q);
        _debye_huckel_RC = x * (q * x + 3. * q * l) / (q * (x + l));

        number debyecut;
        if (this->_grooving) {
            debyecut = 2.0f * sqrt((POS_MM_BACK1) * (POS_MM_BACK1) + (POS_MM_BACK2) * (POS_MM_BACK2)) + _debye_huckel_RC;
        } else {
            debyecut = 2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
        }

        // the cutoff radius for the potential should be the larger of rcut and debyecut
        if (debyecut > this->_rcut){
            this->_rcut = debyecut;
            this->_sqr_rcut = debyecut*debyecut;
        }
        //End copy from DNA2Interaction

        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );
    }
    //Constants for DNA/Protein Interactions
    //NEW VERSION #QuarticExcludedVolume
    //Backbone-Protein Excluded Volume Parameters
    _pro_backbone_sigma = 0.57f;
    _pro_backbone_rstar= 0.569f;
    _pro_backbone_b = 178699253.5f;
    _pro_backbone_rcut = 0.572934f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.36f;
    _pro_base_rstar= 0.359f;
    _pro_base_b = 296866090.0f;
    _pro_base_rcut = 0.362897f;
    //Protein-Protein Excluded Volume Parameters
    _pro_sigma = 0.35f;
    _pro_rstar = 0.349f;
    _pro_b = 306484596.0f;
    _pro_rcut = 0.352894;

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_sigma, &_pro_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rstar, &_pro_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rc, &_pro_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_b, &_pro_b, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_sigma, &_pro_backbone_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rstar, &_pro_backbone_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rc, &_pro_backbone_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_b, &_pro_backbone_b, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_sigma, &_pro_base_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rstar, &_pro_base_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rc, &_pro_base_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_b, &_pro_base_b, sizeof(float)) );

    //Parameters for DNANM book keeping
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_ndna, &this->ndna, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_npro, &this->npro, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_offset, &this->offset, sizeof(int)) );
}

template<typename number, typename number4>
void CUDADNANMInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
        if (_v_lists->use_edge()) {
            dnanm_forces_edge_nonbonded<number, number4>
                    << < (_v_lists->_N_edges - 1) / (this->_launch_cfg.threads_per_block) + 1,
                    this->_launch_cfg.threads_per_block >> >
                    (d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, d_box);

            this->_sum_edge_forces_torques(d_forces, d_torques);

            // potential for removal here
            cudaThreadSynchronize();
            CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

            dnanm_forces_edge_bonded<number, number4>
                    << < this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block >> >
                   (d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box, _d_aff_eqdist, _d_aff_gamma, _d_affected_indx, _d_affected);

        } else throw oxDNAException("Edge Approach is only implemented for DNANM Interaction using CUDA approach. Please add use_edge = 1 to your input file.");

	} else throw oxDNAException("Must Use with Lists to run simulation");
}

template class CUDADNANMInteraction<float, float4>;
template class CUDADNANMInteraction<double, LR_double4>;
