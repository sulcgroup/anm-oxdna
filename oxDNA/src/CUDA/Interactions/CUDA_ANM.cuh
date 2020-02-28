/* System constants */
/*
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];

__constant__ float MD_hb_multi[1];
__constant__ float MD_F1_A[2];
__constant__ float MD_F1_RC[2];
__constant__ float MD_F1_R0[2];
__constant__ float MD_F1_BLOW[2];
__constant__ float MD_F1_BHIGH[2];
__constant__ float MD_F1_RLOW[2];
__constant__ float MD_F1_RHIGH[2];
__constant__ float MD_F1_RCLOW[2];
__constant__ float MD_F1_RCHIGH[2];
// 50 = 2 * 5 * 5
__constant__ float MD_F1_EPS[50];
__constant__ float MD_F1_SHIFT[50];

__constant__ float MD_F2_K[2];
__constant__ float MD_F2_RC[2];
__constant__ float MD_F2_R0[2];
__constant__ float MD_F2_BLOW[2];
__constant__ float MD_F2_RLOW[2];
__constant__ float MD_F2_RCLOW[2];
__constant__ float MD_F2_BHIGH[2];
__constant__ float MD_F2_RCHIGH[2];
__constant__ float MD_F2_RHIGH[2];

__constant__ float MD_F5_PHI_A[4];
__constant__ float MD_F5_PHI_B[4];
__constant__ float MD_F5_PHI_XC[4];
__constant__ float MD_F5_PHI_XS[4];

__constant__ float MD_dh_RC[1];
__constant__ float MD_dh_RHIGH[1];
__constant__ float MD_dh_prefactor[1];
__constant__ float MD_dh_B[1];
__constant__ float MD_dh_minus_kappa[1];
__constant__ bool MD_dh_half_charged_ends[1];
*/

//added constants for DNANM
__constant__ float MD_pro_backbone_sigma;
__constant__ float MD_pro_backbone_rstar;
__constant__ float MD_pro_backbone_b;
__constant__ float MD_pro_backbone_rc;
__constant__ float MD_pro_backbone_stiffness;

__constant__ float MD_pro_base_sigma;
__constant__ float MD_pro_base_rstar;
__constant__ float MD_pro_base_b;
__constant__ float MD_pro_base_rc;
__constant__ float MD_pro_base_stiffness;

__constant__ float MD_pro_sigma;
__constant__ float MD_pro_rstar;
__constant__ float MD_pro_b;
__constant__ float MD_pro_rc;

__constant__ int _npro;
__constant__ int _ndna;
__constant__ int _offset;

#include "../cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__global__ void dnanm_forces_edge_nonbonded(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, edge_bond *edge_list, int n_edges, LR_bonds *bonds, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, CUDABox<number, number4> *box) {
    if (IND >= n_edges) return;

    number4 dF = make_number4 < number, number4>(0, 0, 0, 0);
    number4 dT = make_number4 < number, number4>(0, 0, 0, 0);

    edge_bond b = edge_list[IND];

    // get info for particles
    number4 ppos = poss[b.from];
    number4 qpos = poss[b.to];
    int pbtype = get_particle_btype < number, number4>(ppos);
    int qbtype = get_particle_btype < number, number4>(qpos);

    // particle axes according to Allen's paper
    number4 a1, a2, a3;
    get_vectors_from_quat<number, number4>(orientations[b.from], a1, a2, a3);
    // get info for particle 2
    number4 b1, b2, b3;
    get_vectors_from_quat<number, number4>(orientations[b.to], b1, b2, b3);

    if (pbtype < 0 && qbtype < 0) {
        //Protein-Protein Excluded Volume **NONBONDED
        //copy over dnanm constants? done!
        number4 r = box->minimum_image(ppos, qpos);
        _excluded_volume(r, dF, MD_pro_sigma, MD_pro_rstar, MD_pro_b, MD_pro_rc);
        dF.w *= (number) 0.5f;
        //Add force to p index
        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
        // Should this be a different index?
        //int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //Add force to q index
        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        //int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);

    } else if (pbtype >= 0 && qbtype < 0) {
        //Protein-DNA Excluded Volume **NONBONDED
        number4 ppos_back = POS_BACK * a1;
        number4 ppos_base = POS_BASE * a1;
        number4 rback = box->minimum_image(ppos_back, qpos);
        number4 rbase = box->minimum_image(ppos_base, qpos);
        number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);
        _excluded_volume(rback, Ftmp, MD_pro_backbone_sigma, MD_pro_backbone_rstar, MD_pro_backbone_b, MD_pro_backbone_rc);
        dT += _cross<number, number4>(ppos_back, Ftmp);
        dF += Ftmp;

        _excluded_volume(rbase, Ftmp, MD_pro_base_sigma, MD_pro_base_rstar, MD_pro_base_b, MD_pro_base_rc);
        dT += _cross<number, number4>(ppos_base, Ftmp); //same here
        dF += Ftmp;

        dF.w *= (number) 0.5f;
        dT.w *= (number) 0.5f;
        // Add force AND TORQUE to p index
        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
        //int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[from_index]), dT);

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //Add force to q index
        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        //int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);

    } else if(pbtype < 0 && qbtype >= 0) {
        number4 qpos_back = POS_BACK * b1;
        number4 qpos_base = POS_BASE * b1;
        number4 rback = box->minimum_image(ppos, qpos_back);
        number4 rbase = box->minimum_image(ppos, qpos_base);
        number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);

        _excluded_volume(rback, Ftmp, MD_pro_backbone_sigma, MD_pro_backbone_rstar, MD_pro_backbone_b, MD_pro_backbone_rc);
        dT += _cross<number, number4>(qpos_back, Ftmp); //check the shit out of this
        dF += Ftmp;

        //Ftmp is set to 0 in excluded volume function
        _excluded_volume(rbase, Ftmp, MD_pro_base_sigma, MD_pro_base_rstar, MD_pro_base_b, MD_pro_base_rc);
        dT += _cross<number, number4>(qpos_base, Ftmp);
        dF += Ftmp;

        // NEED TO HANDLE TORQUE **SHOULD THIS BE NEGATIVE??
        dF.w *= (number) 0.5f;
        dT.w *= (number) 0.5f;
        // Add Force to p index
        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
        //int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);

        //NO TORQUE ON PROTEIN PARTICLES
        // Allen Eq. 6 pag 3: DO I NEED THIS WITH HOW TORQUE IS CALCULATED HERE
        //number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
        //number4 crx = _cross < number, number4> (dr, dF);
        //dT.x = -dT.x + crx.x;
        //dT.y = -dT.y + crx.y;
        //dT.z = -dT.z + crx.z;

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //Add Force AND TORQUE to q index
        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        //int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);
        //if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
        //    LR_atomicAddXYZ(&(torques[to_index]), dT);

    } else if(pbtype >= 0 && qbtype >= 0){
        LR_bonds pbonds = bonds[b.from];
        LR_bonds qbonds = bonds[b.to];
        _particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving,
                                                        use_debye_huckel, use_oxDNA2_coaxial_stacking, pbonds, qbonds,
                                                        b.from, b.to, box);

        dF.w *= (number) 0.5f;
        dT.w *= (number) 0.5f;

        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
        //int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[from_index]), dT);

        // Allen Eq. 6 pag 3:
        number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
        number4 crx = _cross < number, number4> (dr, dF);
        dT.x = -dT.x + crx.x;
        dT.y = -dT.y + crx.y;
        dT.z = -dT.z + crx.z;

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        //int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[to_index]), dT);
    }
}

// bonded interactions for edge-based approach
template <typename number, typename number4>
__global__ void dnanm_forces_edge_bonded(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, LR_bonds *bonds, bool grooving, bool use_oxDNA2_FENE, bool use_mbf, number mbf_xmax, number mbf_finf, CUDABox<number, number4> *box, number *_d_spring_eqdist, number *_d_spring_potential) {
    if(IND >= MD_N[0]) return;

    number4 F0, T0;

    F0.x = forces[IND].x;
    F0.y = forces[IND].y;
    F0.z = forces[IND].z;
    F0.w = forces[IND].w;
    T0.x = torques[IND].x;
    T0.y = torques[IND].y;
    T0.z = torques[IND].z;
    T0.w = torques[IND].w;

    number4 dF = make_number4<number, number4>(0, 0, 0, 0);
    number4 dT = make_number4<number, number4>(0, 0, 0, 0);
    number4 ppos = poss[IND];
    // get btype of p
    int pbtype = get_particle_btype <number, number4>(ppos);
    if (pbtype >= 0){
        LR_bonds bs = bonds[IND];
        // particle axes according to Allen's paper

        number4 a1, a2, a3;
        get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

        if(bs.n3 != P_INVALID) {
            number4 qpos = poss[bs.n3];

            number4 b1, b2, b3;
            get_vectors_from_quat<number, number4>(orientations[bs.n3], b1, b2, b3);

            _bonded_part<number, number4, true>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving,
                                                use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
        }
        if(bs.n5 != P_INVALID) {
            number4 qpos = poss[bs.n5];

            number4 b1, b2, b3;
            get_vectors_from_quat<number, number4>(orientations[bs.n5], b1, b2, b3);

            _bonded_part<number, number4, false>(qpos, b1, b2, b3, ppos, a1, a2, a3, dF, dT, grooving,
                                                 use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
        }

        // the real energy per particle is half of the one computed (because we count each interaction twice)
        dF.w *= (number) 0.5f;
        dT.w *= (number) 0.5f;

        forces[IND] = (dF + F0);
        torques[IND] = (dT + T0);

        torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, torques[IND]);
    } else{
//        printf("p %d btype %d\n", IND, pbtype);
        number4 ftotal = make_number4<number, number4>(0,0,0,0);
        for(int i = _npro*(IND - _offset); i < _npro*(IND - _offset)+_npro; i++){
            int qindex = i - _npro*(IND - _offset) +_offset;
            number eqdist = _d_spring_eqdist[i];
            if(eqdist > (number) 0){
                number4 qpos = poss[qindex];
                number4 r = make_number4<number, number4>(qpos.x - ppos.x, qpos.y - ppos.y, qpos.z - ppos.z, (number) 0);
                number d = _module<number, number4>(r);
                number gamma = _d_spring_potential[i];

                number fmod = (-1.0f * gamma) * (d - eqdist) / d;

                dF = fmod*r;
                dF.w = -0.5f * gamma * powf(d-eqdist, 2);

                ftotal -= dF;
    //                if(IND == 217 && qindex == 550) {
                //printf("p %d q %d g %.2f d %.6f ro %.2f df.x %.8f, df.y %.8f, df.z %.8f p.x %.8f p.y %.8f p.z %.8f q.x %.8f q.y %.8f q.z %.8f\n",
                  //     IND, qindex, gamma, d, eqdist, dF.x, dF.y, dF.z, forces[IND].x, forces[IND].y, forces[IND].z,
                    //   forces[qindex].x,forces[qindex].y, forces[qindex].z);
            };
        }
        LR_atomicAddXYZ(&(forces[IND]), ftotal);
        //printf("IND %d, f.x %.6f, f.y %.6f, f.z %.6f\n", IND, forces[IND].x, forces[IND].y, forces[IND].z);
    }
}

