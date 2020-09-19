/*
 * CUDARNACTInteraction.cu
 *
 *  Created on: 25 Aug 2020
 *      Author: jonah
 */

#include "CUDARNACTInteraction.h"
#include "CUDARNAInteraction.h" //questionable
#include "CUDA_RNA.cuh"
#include "CUDA_ACT_R.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/RNAInteraction2.h"
#include "../../Interactions/RNACTInteraction.h"

void copy_Model_to_CUDAModel_RNACT(Model& model_from, CUDAModel& model_to)
{
    model_to.RNA_POS_BACK =  model_from.RNA_POS_BACK ;
    model_to.RNA_POS_STACK =  model_from.RNA_POS_STACK ;
    model_to.RNA_POS_BASE =  model_from.RNA_POS_BASE ;
    model_to.RNA_GAMMA =  model_from.RNA_GAMMA ;
    model_to.RNA_POS_STACK_3_a1 =  model_from.RNA_POS_STACK_3_a1 ;
    model_to.RNA_POS_STACK_3_a2 =  model_from.RNA_POS_STACK_3_a2 ;
    model_to.RNA_POS_STACK_5_a1 =  model_from.RNA_POS_STACK_5_a1 ;
    model_to.RNA_POS_STACK_5_a2 =  model_from.RNA_POS_STACK_5_a2 ;
    model_to.RNA_FENE_EPS =  model_from.RNA_FENE_EPS ;
    model_to.RNA_FENE_R0 =  model_from.RNA_FENE_R0 ;
    model_to.RNA_FENE_DELTA =  model_from.RNA_FENE_DELTA ;
    model_to.RNA_FENE_DELTA2 =  model_from.RNA_FENE_DELTA2 ;
    model_to.RNA_EXCL_EPS =  model_from.RNA_EXCL_EPS ;
    model_to.RNA_EXCL_S1 =  model_from.RNA_EXCL_S1 ;
    model_to.RNA_EXCL_S2 =  model_from.RNA_EXCL_S2 ;
    model_to.RNA_EXCL_S3 =  model_from.RNA_EXCL_S3 ;
    model_to.RNA_EXCL_S4 =  model_from.RNA_EXCL_S4 ;
    model_to.RNA_EXCL_R1 =  model_from.RNA_EXCL_R1 ;
    model_to.RNA_EXCL_R2 =  model_from.RNA_EXCL_R2 ;
    model_to.RNA_EXCL_R3 =  model_from.RNA_EXCL_R3 ;
    model_to.RNA_EXCL_R4 =  model_from.RNA_EXCL_R4 ;
    model_to.RNA_EXCL_B1 =  model_from.RNA_EXCL_B1 ;
    model_to.RNA_EXCL_B2 =  model_from.RNA_EXCL_B2 ;
    model_to.RNA_EXCL_B3 =  model_from.RNA_EXCL_B3 ;
    model_to.RNA_EXCL_B4 =  model_from.RNA_EXCL_B4 ;
    model_to.RNA_EXCL_RC1 =  model_from.RNA_EXCL_RC1 ;
    model_to.RNA_EXCL_RC2 =  model_from.RNA_EXCL_RC2 ;
    model_to.RNA_EXCL_RC3 =  model_from.RNA_EXCL_RC3 ;
    model_to.RNA_EXCL_RC4 =  model_from.RNA_EXCL_RC4 ;
    model_to.RNA_HYDR_EPS =  model_from.RNA_HYDR_EPS ;
    model_to.RNA_HYDR_A =  model_from.RNA_HYDR_A ;
    model_to.RNA_HYDR_RC =  model_from.RNA_HYDR_RC ;
    model_to.RNA_HYDR_R0 =  model_from.RNA_HYDR_R0 ;
    model_to.RNA_HYDR_BLOW =  model_from.RNA_HYDR_BLOW ;
    model_to.RNA_HYDR_BHIGH =  model_from.RNA_HYDR_BHIGH ;
    model_to.RNA_HYDR_RLOW =  model_from.RNA_HYDR_RLOW ;
    model_to.RNA_HYDR_RHIGH =  model_from.RNA_HYDR_RHIGH ;
    model_to.RNA_HYDR_RCLOW =  model_from.RNA_HYDR_RCLOW ;
    model_to.RNA_HYDR_RCHIGH =  model_from.RNA_HYDR_RCHIGH ;
    model_to.RNA_HYDR_THETA1_A =  model_from.RNA_HYDR_THETA1_A ;
    model_to.RNA_HYDR_THETA1_B =  model_from.RNA_HYDR_THETA1_B ;
    model_to.RNA_HYDR_THETA1_T0 =  model_from.RNA_HYDR_THETA1_T0 ;
    model_to.RNA_HYDR_THETA1_TS =  model_from.RNA_HYDR_THETA1_TS ;
    model_to.RNA_HYDR_THETA1_TC =  model_from.RNA_HYDR_THETA1_TC ;
    model_to.RNA_HYDR_THETA2_A =  model_from.RNA_HYDR_THETA2_A ;
    model_to.RNA_HYDR_THETA2_B =  model_from.RNA_HYDR_THETA2_B ;
    model_to.RNA_HYDR_THETA2_T0 =  model_from.RNA_HYDR_THETA2_T0 ;
    model_to.RNA_HYDR_THETA2_TS =  model_from.RNA_HYDR_THETA2_TS ;
    model_to.RNA_HYDR_THETA2_TC =  model_from.RNA_HYDR_THETA2_TC ;
    model_to.RNA_HYDR_THETA3_A =  model_from.RNA_HYDR_THETA3_A ;
    model_to.RNA_HYDR_THETA3_B =  model_from.RNA_HYDR_THETA3_B ;
    model_to.RNA_HYDR_THETA3_T0 =  model_from.RNA_HYDR_THETA3_T0 ;
    model_to.RNA_HYDR_THETA3_TS =  model_from.RNA_HYDR_THETA3_TS ;
    model_to.RNA_HYDR_THETA3_TC =  model_from.RNA_HYDR_THETA3_TC ;
    model_to.RNA_HYDR_THETA4_A =  model_from.RNA_HYDR_THETA4_A ;
    model_to.RNA_HYDR_THETA4_B =  model_from.RNA_HYDR_THETA4_B ;
    model_to.RNA_HYDR_THETA4_T0 =  model_from.RNA_HYDR_THETA4_T0 ;
    model_to.RNA_HYDR_THETA4_TS =  model_from.RNA_HYDR_THETA4_TS ;
    model_to.RNA_HYDR_THETA4_TC =  model_from.RNA_HYDR_THETA4_TC ;
    model_to.RNA_HYDR_THETA7_A =  model_from.RNA_HYDR_THETA7_A ;
    model_to.RNA_HYDR_THETA7_B =  model_from.RNA_HYDR_THETA7_B ;
    model_to.RNA_HYDR_THETA7_T0 =  model_from.RNA_HYDR_THETA7_T0 ;
    model_to.RNA_HYDR_THETA7_TS =  model_from.RNA_HYDR_THETA7_TS ;
    model_to.RNA_HYDR_THETA7_TC =  model_from.RNA_HYDR_THETA7_TC ;
    model_to.RNA_HYDR_THETA8_A =  model_from.RNA_HYDR_THETA8_A ;
    model_to.RNA_HYDR_THETA8_B =  model_from.RNA_HYDR_THETA8_B ;
    model_to.RNA_HYDR_THETA8_T0 =  model_from.RNA_HYDR_THETA8_T0 ;
    model_to.RNA_HYDR_THETA8_TS =  model_from.RNA_HYDR_THETA8_TS ;
    model_to.RNA_HYDR_THETA8_TC =  model_from.RNA_HYDR_THETA8_TC ;
    model_to.RNA_STCK_BASE_EPS =  model_from.RNA_STCK_BASE_EPS ;
    model_to.RNA_STCK_FACT_EPS =  model_from.RNA_STCK_FACT_EPS ;
    model_to.RNA_STCK_A =  model_from.RNA_STCK_A ;
    model_to.RNA_STCK_RC =  model_from.RNA_STCK_RC ;
    model_to.RNA_STCK_R0 =  model_from.RNA_STCK_R0 ;
    model_to.RNA_STCK_BLOW =  model_from.RNA_STCK_BLOW ;
    model_to.RNA_STCK_BHIGH =  model_from.RNA_STCK_BHIGH ;
    model_to.RNA_STCK_RLOW =  model_from.RNA_STCK_RLOW ;
    model_to.RNA_STCK_RHIGH =  model_from.RNA_STCK_RHIGH ;
    model_to.RNA_STCK_RCLOW =  model_from.RNA_STCK_RCLOW ;
    model_to.RNA_STCK_RCHIGH =  model_from.RNA_STCK_RCHIGH ;
    model_to.RNA_STCK_THETA4_A =  model_from.RNA_STCK_THETA4_A ;
    model_to.RNA_STCK_THETA4_B =  model_from.RNA_STCK_THETA4_B ;
    model_to.RNA_STCK_THETA4_T0 =  model_from.RNA_STCK_THETA4_T0 ;
    model_to.RNA_STCK_THETA4_TS =  model_from.RNA_STCK_THETA4_TS ;
    model_to.RNA_STCK_THETA4_TC =  model_from.RNA_STCK_THETA4_TC ;
    model_to.RNA_STCK_THETA5_A =  model_from.RNA_STCK_THETA5_A ;
    model_to.RNA_STCK_THETA5_B =  model_from.RNA_STCK_THETA5_B ;
    model_to.RNA_STCK_THETA5_T0 =  model_from.RNA_STCK_THETA5_T0 ;
    model_to.RNA_STCK_THETA5_TS =  model_from.RNA_STCK_THETA5_TS ;
    model_to.RNA_STCK_THETA5_TC =  model_from.RNA_STCK_THETA5_TC ;
    model_to.RNA_STCK_THETA6_A =  model_from.RNA_STCK_THETA6_A ;
    model_to.RNA_STCK_THETA6_B =  model_from.RNA_STCK_THETA6_B ;
    model_to.RNA_STCK_THETA6_T0 =  model_from.RNA_STCK_THETA6_T0 ;
    model_to.RNA_STCK_THETA6_TS =  model_from.RNA_STCK_THETA6_TS ;
    model_to.RNA_STCK_THETA6_TC =  model_from.RNA_STCK_THETA6_TC ;
    model_to.STCK_THETAB1_A =  model_from.STCK_THETAB1_A ;
    model_to.STCK_THETAB1_B =  model_from.STCK_THETAB1_B ;
    model_to.STCK_THETAB1_T0 =  model_from.STCK_THETAB1_T0 ;
    model_to.STCK_THETAB1_TS =  model_from.STCK_THETAB1_TS ;
    model_to.STCK_THETAB1_TC =  model_from.STCK_THETAB1_TC ;
    model_to.STCK_THETAB2_A =  model_from.STCK_THETAB2_A ;
    model_to.STCK_THETAB2_B =  model_from.STCK_THETAB2_B ;
    model_to.STCK_THETAB2_T0 =  model_from.STCK_THETAB2_T0 ;
    model_to.STCK_THETAB2_TS =  model_from.STCK_THETAB2_TS ;
    model_to.STCK_THETAB2_TC =  model_from.STCK_THETAB2_TC ;
    model_to.RNA_STCK_PHI1_A =  model_from.RNA_STCK_PHI1_A ;
    model_to.RNA_STCK_PHI1_B =  model_from.RNA_STCK_PHI1_B ;
    model_to.RNA_STCK_PHI1_XC =  model_from.RNA_STCK_PHI1_XC ;
    model_to.RNA_STCK_PHI1_XS =  model_from.RNA_STCK_PHI1_XS ;
    model_to.RNA_STCK_PHI2_A =  model_from.RNA_STCK_PHI2_A ;
    model_to.RNA_STCK_PHI2_B =  model_from.RNA_STCK_PHI2_B ;
    model_to.RNA_STCK_PHI2_XC =  model_from.RNA_STCK_PHI2_XC ;
    model_to.RNA_STCK_PHI2_XS =  model_from.RNA_STCK_PHI2_XS ;
    model_to.RNA_CRST_R0 =  model_from.RNA_CRST_R0 ;
    model_to.RNA_CRST_RC =  model_from.RNA_CRST_RC ;
    model_to.RNA_CRST_K =  model_from.RNA_CRST_K ;
    model_to.RNA_CRST_BLOW =  model_from.RNA_CRST_BLOW ;
    model_to.RNA_CRST_RLOW =  model_from.RNA_CRST_RLOW ;
    model_to.RNA_CRST_RCLOW =  model_from.RNA_CRST_RCLOW ;
    model_to.RNA_CRST_BHIGH =  model_from.RNA_CRST_BHIGH ;
    model_to.RNA_CRST_RHIGH =  model_from.RNA_CRST_RHIGH ;
    model_to.RNA_CRST_RCHIGH =  model_from.RNA_CRST_RCHIGH ;
    model_to.RNA_CRST_THETA1_A =  model_from.RNA_CRST_THETA1_A ;
    model_to.RNA_CRST_THETA1_B =  model_from.RNA_CRST_THETA1_B ;
    model_to.RNA_CRST_THETA1_T0 =  model_from.RNA_CRST_THETA1_T0 ;
    model_to.RNA_CRST_THETA1_TS =  model_from.RNA_CRST_THETA1_TS ;
    model_to.RNA_CRST_THETA1_TC =  model_from.RNA_CRST_THETA1_TC ;
    model_to.RNA_CRST_THETA2_A =  model_from.RNA_CRST_THETA2_A ;
    model_to.RNA_CRST_THETA2_B =  model_from.RNA_CRST_THETA2_B ;
    model_to.RNA_CRST_THETA2_T0 =  model_from.RNA_CRST_THETA2_T0 ;
    model_to.RNA_CRST_THETA2_TS =  model_from.RNA_CRST_THETA2_TS ;
    model_to.RNA_CRST_THETA2_TC =  model_from.RNA_CRST_THETA2_TC ;
    model_to.RNA_CRST_THETA3_A =  model_from.RNA_CRST_THETA3_A ;
    model_to.RNA_CRST_THETA3_B =  model_from.RNA_CRST_THETA3_B ;
    model_to.RNA_CRST_THETA3_T0 =  model_from.RNA_CRST_THETA3_T0 ;
    model_to.RNA_CRST_THETA3_TS =  model_from.RNA_CRST_THETA3_TS ;
    model_to.RNA_CRST_THETA3_TC =  model_from.RNA_CRST_THETA3_TC ;
    model_to.RNA_CRST_THETA4_A =  model_from.RNA_CRST_THETA4_A ;
    model_to.RNA_CRST_THETA4_B =  model_from.RNA_CRST_THETA4_B ;
    model_to.RNA_CRST_THETA4_T0 =  model_from.RNA_CRST_THETA4_T0 ;
    model_to.RNA_CRST_THETA4_TS =  model_from.RNA_CRST_THETA4_TS ;
    model_to.RNA_CRST_THETA4_TC =  model_from.RNA_CRST_THETA4_TC ;
    model_to.RNA_CRST_THETA7_A =  model_from.RNA_CRST_THETA7_A ;
    model_to.RNA_CRST_THETA7_B =  model_from.RNA_CRST_THETA7_B ;
    model_to.RNA_CRST_THETA7_T0 =  model_from.RNA_CRST_THETA7_T0 ;
    model_to.RNA_CRST_THETA7_TS =  model_from.RNA_CRST_THETA7_TS ;
    model_to.RNA_CRST_THETA7_TC =  model_from.RNA_CRST_THETA7_TC ;
    model_to.RNA_CRST_THETA8_A =  model_from.RNA_CRST_THETA8_A ;
    model_to.RNA_CRST_THETA8_B =  model_from.RNA_CRST_THETA8_B ;
    model_to.RNA_CRST_THETA8_T0 =  model_from.RNA_CRST_THETA8_T0 ;
    model_to.RNA_CRST_THETA8_TS =  model_from.RNA_CRST_THETA8_TS ;
    model_to.RNA_CRST_THETA8_TC =  model_from.RNA_CRST_THETA8_TC ;
    model_to.RNA_CXST_R0 =  model_from.RNA_CXST_R0 ;
    model_to.RNA_CXST_RC =  model_from.RNA_CXST_RC ;
    model_to.RNA_CXST_K =  model_from.RNA_CXST_K ;
    model_to.RNA_CXST_BLOW =  model_from.RNA_CXST_BLOW ;
    model_to.RNA_CXST_RLOW =  model_from.RNA_CXST_RLOW ;
    model_to.RNA_CXST_RCLOW =  model_from.RNA_CXST_RCLOW ;
    model_to.RNA_CXST_BHIGH =  model_from.RNA_CXST_BHIGH ;
    model_to.RNA_CXST_RHIGH =  model_from.RNA_CXST_RHIGH ;
    model_to.RNA_CXST_RCHIGH =  model_from.RNA_CXST_RCHIGH ;
    model_to.RNA_CXST_THETA1_A =  model_from.RNA_CXST_THETA1_A ;
    model_to.RNA_CXST_THETA1_B =  model_from.RNA_CXST_THETA1_B ;
    model_to.RNA_CXST_THETA1_T0 =  model_from.RNA_CXST_THETA1_T0 ;
    model_to.RNA_CXST_THETA1_TS =  model_from.RNA_CXST_THETA1_TS ;
    model_to.RNA_CXST_THETA1_TC =  model_from.RNA_CXST_THETA1_TC ;
    model_to.RNA_CXST_THETA4_A =  model_from.RNA_CXST_THETA4_A ;
    model_to.RNA_CXST_THETA4_B =  model_from.RNA_CXST_THETA4_B ;
    model_to.RNA_CXST_THETA4_T0 =  model_from.RNA_CXST_THETA4_T0 ;
    model_to.RNA_CXST_THETA4_TS =  model_from.RNA_CXST_THETA4_TS ;
    model_to.RNA_CXST_THETA4_TC =  model_from.RNA_CXST_THETA4_TC ;
    model_to.RNA_CXST_THETA5_A =  model_from.RNA_CXST_THETA5_A ;
    model_to.RNA_CXST_THETA5_B =  model_from.RNA_CXST_THETA5_B ;
    model_to.RNA_CXST_THETA5_T0 =  model_from.RNA_CXST_THETA5_T0 ;
    model_to.RNA_CXST_THETA5_TS =  model_from.RNA_CXST_THETA5_TS ;
    model_to.RNA_CXST_THETA5_TC =  model_from.RNA_CXST_THETA5_TC ;
    model_to.RNA_CXST_THETA6_A =  model_from.RNA_CXST_THETA6_A ;
    model_to.RNA_CXST_THETA6_B =  model_from.RNA_CXST_THETA6_B ;
    model_to.RNA_CXST_THETA6_T0 =  model_from.RNA_CXST_THETA6_T0 ;
    model_to.RNA_CXST_THETA6_TS =  model_from.RNA_CXST_THETA6_TS ;
    model_to.RNA_CXST_THETA6_TC =  model_from.RNA_CXST_THETA6_TC ;
    model_to.RNA_CXST_PHI3_A =  model_from.RNA_CXST_PHI3_A ;
    model_to.RNA_CXST_PHI3_B =  model_from.RNA_CXST_PHI3_B ;
    model_to.RNA_CXST_PHI3_XC =  model_from.RNA_CXST_PHI3_XC ;
    model_to.RNA_CXST_PHI3_XS =  model_from.RNA_CXST_PHI3_XS ;
    model_to.RNA_CXST_PHI4_A =  model_from.RNA_CXST_PHI4_A ;
    model_to.RNA_CXST_PHI4_B =  model_from.RNA_CXST_PHI4_B ;
    model_to.RNA_CXST_PHI4_XC =  model_from.RNA_CXST_PHI4_XC ;
    model_to.RNA_CXST_PHI4_XS =  model_from.RNA_CXST_PHI4_XS ;

    model_to.p3_x = model_from.p3_x;
    model_to.p3_y = model_from.p3_y;
    model_to.p3_z = model_from.p3_z;

    model_to.p5_x = model_from.p5_x;
    model_to.p5_y = model_from.p5_y;
    model_to.p5_z = model_from.p5_z;

    model_to.RNA_POS_BACK_a1 = model_from.RNA_POS_BACK_a1;
    model_to.RNA_POS_BACK_a2 = model_from.RNA_POS_BACK_a2;
    model_to.RNA_POS_BACK_a3 = model_from.RNA_POS_BACK_a3;



}





template<typename number, typename number4>
CUDARNACTInteraction<number, number4>::CUDARNACTInteraction() {

    _grooving = false;
    _read_par = true;

    //Not copied over to device memory
    _spring_potential = NULL;
    _spring_eqdist = NULL;
    _affected_len = NULL;

    //Copied over to device memory
    _h_ang_params = NULL;
    _d_ang_params = NULL;

    _h_affected_indx = NULL;
    _d_affected_indx = NULL;

    _h_affected = NULL;
    _d_affected = NULL;

    _h_aff_eqdist = NULL;
    _d_aff_eqdist = NULL;

    _h_aff_gamma = NULL;
    _d_aff_gamma = NULL;

    _spring_param_size_number = 0;
    _ang_param_size = 0;
}

template<typename number, typename number4>
CUDARNACTInteraction<number, number4>::~CUDARNACTInteraction() {
    //Delete All pointers required for spring potential parameters
    if(_spring_potential != NULL) delete[] _spring_potential;
    if(_spring_eqdist != NULL) delete[] _spring_eqdist;

    if(_h_ang_params != NULL) delete[] _h_ang_params;
    if(_d_ang_params != NULL) CUDA_SAFE_CALL( cudaFree(_d_ang_params) );

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
void CUDARNACTInteraction<number, number4>::get_settings(input_file &inp) {
    _use_debye_huckel = false;
    _mismatch_repulsion = false;
    std::string inter_type;
    if (!getInputString(&inp, "parfile", this->_parameterfile, 0) == KEY_FOUND){
        throw oxDNAException("Key 'parfile' not found. Necessary for Protein sims.");

    char s[5] = "none";
    if(strcmp(this->_parameterfile, s) == 0) _read_par = false;

    }
    if (!getInputString(&inp, "topology", this->_topology_filename, 0) == KEY_FOUND){
        throw oxDNAException("Key 'topology_file' not found.");
    }

    //Function pointers constructed in the implicit RNACT Constructor Call

    _use_debye_huckel = true;
    // copy-pasted from the DNA2Interaction constructor
    _debye_huckel_half_charged_ends = true;
    this->_grooving = true;
    // end copy from DNA2Interaction

    // copied from DNA2Interaction::get_settings() (CPU), the least bad way of doing things
    getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
    getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

    // lambda-factor (the dh length at T = 300K, I = 1.0)
    _debye_huckel_lambdafactor =  0.3667258;
    getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

    // the prefactor to the Debye-Huckel term
    _debye_huckel_prefactor = 0.0858;
    getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
    // End copy from DNA2Interaction

    int mismatches = 0;
    if(getInputBoolAsInt(&inp, "mismatch_repulsion", &mismatches, 0) == KEY_FOUND) {
        this->_mismatch_repulsion = (bool) mismatches;
    }
    if(this->_mismatch_repulsion)
    {
        float temp;
        if(getInputFloat(&inp, "mismatch_repulsion_strength", &temp, 0) == KEY_FOUND) {
            this->_RNA_HYDR_MIS = temp;
        }
        else
        {
            this->_RNA_HYDR_MIS = 1;
        }

    }
    this->RNACTInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDARNACTInteraction<number, number4>::cuda_init(number box_side, int N) {
    CUDABaseInteraction<number, number4>::cuda_init(box_side, N);

//    Addition of Reading Parameter File -> Moved from get_settings due to needing to fill variables that are filled in the CPU version of DNACTInteraction::read_topology
    fstream top;
    int tmp1, tmp2;
    top.open(this->_topology_filename, ios::in);
    if (top.is_open()){
        top >> tmp1 >> tmp2 >> this->nrna >> this->npro >> this->nrnas;
        top >> this->_firststrand;
        top.close();
    } else {
        throw oxDNAException("Could not open Topology File");
    }


    if(this->_firststrand < 0) offset = 0;
    else if(this->_firststrand > 0) offset = this->nrna;
    else throw oxDNAException("No Strand should have an ID of 0");

    if (_read_par) {
        //Initalizing Some Host and Device Arrays for Spring Parameters
        _spring_param_size_number = sizeof(number) * (this->npro*this->npro);
        _ang_param_size = sizeof(number) * (this->npro*4);

        _spring_potential = new number[this->npro*this->npro]();
        _spring_eqdist = new number[this->npro*this->npro]();

        //Initializing Host and Device Arrays for Angular Parameters
        _h_ang_params = new number[this->npro*4]();
        CUDA_SAFE_CALL( cudaMalloc(&_d_ang_params, _ang_param_size));

        char potswitch = 'x';
        number potential = 0.f, dist = 0.f;
        for(int i = 0; i< (this->npro*this->npro); i++){
            _spring_eqdist[i] = dist;
            _spring_potential[i] = potential;
            if(i < (this->npro *4 )) _h_ang_params[i] = dist;
        }

        //Checkers as Lambdas
        auto valid_angles = [](double a, double b, double c, double d)
        {
            double anglemin = min({a, b, c, d});
            double anglemax = max({a, b, c, d});
            if (anglemin < -1.0 || anglemax > 1.0){
                throw oxDNAException("Cos of Angle in Parameter File not in Valid bounds");
            }
        };

        auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
            if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
            if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
            if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
            if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
            if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
        };

        //Reading Parameter File
        int key1, key2 = 0;
        number a0, b0, c0, d0;
        string carbons;
        fstream parameters;
        parameters.open(this->_parameterfile, ios::in);
        getline (parameters,carbons);

        //total connections
        int spring_connection_num = 0;

        //allocate and declare affected_len vector
        _affected_len = new int[this->npro]();
        for(int i = 0; i < this->npro; i++) _affected_len[i] = 0;

        //Read Parameter File
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
                //potswitch is currently unused but may be later

                if (key2 - key1 == 1){
                    //Angular Parameters
                    parameters >> a0 >> b0 >> c0 >> d0;
                    valid_angles(a0, b0, c0, d0);
                    _h_ang_params[key1*4] = a0;
                    _h_ang_params[key1*4+1] = b0;
                    _h_ang_params[key1*4+2] = c0;
                    _h_ang_params[key1*4+3] = d0;

                    _spring_potential[key1*this->npro + key2] = potential;
                    _spring_eqdist[key1*this->npro + key2] = dist;

                    _spring_potential[key2*this->npro + key1] = potential;
                    _spring_eqdist[key2*this->npro + key1] = dist;

                } else {
                    _spring_potential[key1*this->npro + key2] = potential;
                    _spring_eqdist[key1*this->npro + key2] = dist;

                    _spring_potential[key2*this->npro + key1] = potential;
                    _spring_eqdist[key2*this->npro + key1] = dist;
                }
            }
            parameters.close();
        } else {
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

        //Parameters for Bending/Torsional, _h_ang_params is filled in parameter file reading
        CUDA_SAFE_CALL(cudaMemcpy(_d_ang_params, _h_ang_params, _ang_param_size, cudaMemcpyHostToDevice));

        //Memory Used by Parameters
        float param_memory_mb = (spring_connection_num * 2 * sizeof(int) + 2 * spring_connection_num * 2 * sizeof(number)
                + (this->npro + 1) * sizeof(int) + 4 * this->npro * sizeof(number))/SQR(1024);
        OX_LOG(Logger::LOG_INFO, "Spring Parameters Size: %.2f MB", param_memory_mb);

    } else OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");

    // Copied from CUDADNAINTERACTION
    RNAInteraction<number>::init();

    float f_copy = 1.0;//this->_hb_multiplier;
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

    //mismatch repulsion modification
    if(this->_mismatch_repulsion)
    {
        float tempmis = -1.0 * this->_RNA_HYDR_MIS / this->model->RNA_HYDR_EPS;
        this->F1_EPS[RNA_HYDR_F1][0][0] *= tempmis;
        this->F1_SHIFT[RNA_HYDR_F1][0][0] *= tempmis;
    }

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



    CUDAModel cudamodel;
    copy_Model_to_CUDAModel_RNACT(*(this->model), cudamodel);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(rnamodel,&cudamodel,sizeof(CUDAModel))  );


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
        _debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(-4.*x*x*x * l * l * q );
        _debye_huckel_RC = x*(q*x + 3. * q* l )/(q * (x+l));


        number debyecut =  2. * sqrt(SQR(this->model->RNA_POS_BACK_a1)	+ SQR(this->model->RNA_POS_BACK_a2) + SQR(this->model->RNA_POS_BACK_a3)) + _debye_huckel_RC;

        // the cutoff radius for the potential should be the larger of rcut and debyecut
        if (debyecut > this->_rcut){
            this->_rcut = debyecut;
            this->_sqr_rcut = debyecut*debyecut;
        }
        // End copy from cpu interaction
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
        CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );
    }
    //Constants for RNA/Protein Interaction (Same as DNA ones)
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

    _kbend = this->_k_bend;
    _ktor = this->_k_tor;

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

    //Parameters for DNACT book keeping
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_nrna, &this->nrna, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_npro, &this->npro, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_offset, &this->offset, sizeof(int)) );

    //kb and kt Parameters
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_kb, &_kbend, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_kt, &_ktor, sizeof(float)) );

}

template<typename number, typename number4>
void CUDARNACTInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
        if (_v_lists->use_edge()) {
            rnact_forces_edge_nonbonded<number, number4>
                    << < (_v_lists->_N_edges - 1) / (this->_launch_cfg.threads_per_block) + 1,
                    this->_launch_cfg.threads_per_block >> >
            (d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges, d_bonds, this->_average, this->_use_debye_huckel, this->_mismatch_repulsion, d_box);

            this->_sum_edge_forces_torques(d_forces, d_torques);

            // potential for removal here
            cudaThreadSynchronize();
            CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

            rnact_forces_edge_bonded<number, number4>
                    << < this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block >> >
            (d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_average, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box, _d_aff_eqdist, _d_aff_gamma, _d_ang_params, _d_affected_indx, _d_affected);

        } else throw oxDNAException("Edge Approach is only implemented for DNACT Interaction using CUDA approach. Please add use_edge = 1 to your input file.");

	} else throw oxDNAException("Must Use with Lists to run simulation");
}

template class CUDARNACTInteraction<float, float4>;
template class CUDARNACTInteraction<double, LR_double4>;
