/* System constants */
//added constants for DNANM
__constant__ float MD_pro_backbone_sigma;
__constant__ float MD_pro_backbone_rstar;
__constant__ float MD_pro_backbone_b;
__constant__ float MD_pro_backbone_rc;

__constant__ float MD_pro_base_sigma;
__constant__ float MD_pro_base_rstar;
__constant__ float MD_pro_base_b;
__constant__ float MD_pro_base_rc;

__constant__ float MD_pro_sigma;
__constant__ float MD_pro_rstar;
__constant__ float MD_pro_b;
__constant__ float MD_pro_rc;

__constant__ int _npro;
__constant__ int _ndna;
__constant__ int _offset;

#include "../cuda_utils/CUDA_lr_common.cuh"

template<typename number, typename number4>
__forceinline__ __device__ void _excluded_volume_quart(const number4 &r, number4 &F, number sigma, number rstar, number b, number rc) {
    number rsqr = CUDA_DOT(r, r);

    F.x = F.y = F.z = F.w = (number) 0.f;
    if(rsqr < SQR(rc)) {
        if(rsqr > SQR(rstar)) {
            number rmod = sqrt(rsqr);
            number rrc = rmod - rc;
            number fmod = 4.f * EXCL_EPS * b * CUB(rrc) / rmod;
            F.x = r.x * fmod;
            F.y = r.y * fmod;
            F.z = r.z * fmod;
            F.w = EXCL_EPS * b * SQR(SQR(rrc));
        }
        else {
            number lj_part = CUB(SQR(sigma)/rsqr);
            number fmod = 24.f * EXCL_EPS * (lj_part - 2.f*SQR(lj_part)) / rsqr;
            F.x = r.x * fmod;
            F.y = r.y * fmod;
            F.z = r.z * fmod;
            F.w = 4.f * EXCL_EPS * (SQR(lj_part) - lj_part);
        }
    }
}

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
        //Protein-Protein Excluded Volume
        number4 r = box->minimum_image(ppos, qpos);
        _excluded_volume_quart(r, dF, MD_pro_sigma, MD_pro_rstar, MD_pro_b, MD_pro_rc);

        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from; //pindex
        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to; //qindex

        //Add force to p index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //Add force to q index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);

    } else if (pbtype >= 0 && qbtype < 0) {
        //Protein-DNA Excluded Volume **NONBONDED
        number4 ppos_back;
        if(grooving) ppos_back = POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2;
        else ppos_back = POS_BACK * a1;
        number4 ppos_base = POS_BASE * a1;

        number4 r = box->minimum_image(ppos, qpos);
        number4 rback = r - ppos_back;
        number4 rbase = r - ppos_base;

        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;

        number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);
        _excluded_volume_quart(rback, Ftmp, MD_pro_backbone_sigma, MD_pro_backbone_rstar, MD_pro_backbone_b, MD_pro_backbone_rc);
        dT += _cross<number, number4>(ppos_back, Ftmp);
        dF += Ftmp;

        _excluded_volume_quart(rbase, Ftmp, MD_pro_base_sigma, MD_pro_base_rstar, MD_pro_base_b, MD_pro_base_rc);
        dT += _cross<number, number4>(ppos_base, Ftmp);
        dF += Ftmp;

        // Add force and torque to p index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[from_index]), dT);

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //Add force to q index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f){
            LR_atomicAddXYZ(&(forces[to_index]), dF);
//            printf("p %d q %d f.x %.4f f.y %.4f f.z %.4f\n", from_index, to_index, -dF.x, -dF.y, -dF.z);
        }
    } else if(pbtype < 0 && qbtype >= 0) {
        number4 qpos_back;
        if(grooving) qpos_back = POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2;
        else qpos_back = POS_BACK * b1;
        number4 qpos_base = POS_BASE * b1;
        number4 qpos_stack = POS_STACK * b1;

        int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
        int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;

        //vector r is reversed from the above case
        number4 r = box->minimum_image(qpos, ppos);
        number4 rback = r - qpos_back;
        number4 rbase = r - qpos_base;

        number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);
        _excluded_volume_quart(rback, Ftmp, MD_pro_backbone_sigma, MD_pro_backbone_rstar, MD_pro_backbone_b, MD_pro_backbone_rc);
        dT += _cross<number, number4>(qpos_back, Ftmp);
        dF += Ftmp;

        _excluded_volume_quart(rbase, Ftmp, MD_pro_base_sigma, MD_pro_base_rstar, MD_pro_base_b, MD_pro_base_rc);
        dT += _cross<number, number4>(qpos_base, Ftmp);
        dF += Ftmp;

        // Add Force and torque to q index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[to_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[to_index]), dT);

        dF.x = -dF.x;
        dF.y = -dF.y;
        dF.z = -dF.z;

        //add forces to p index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f) {
            LR_atomicAddXYZ(&(forces[from_index]), dF);
//            printf("p %d q %d f.x %.4f f.y %.4f f.z %.4f\n", from_index, to_index, -dF.x, -dF.y, -dF.z);
        }
    } else if(pbtype >= 0 && qbtype >= 0){
        LR_bonds pbonds = bonds[b.from];
        LR_bonds qbonds = bonds[b.to];
        _particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving,
                                                        use_debye_huckel, use_oxDNA2_coaxial_stacking, pbonds, qbonds,
                                                        b.from, b.to, box);


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
__global__ void dnanm_forces_edge_bonded(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, LR_bonds *bonds, bool grooving, bool use_oxDNA2_FENE, bool use_mbf, number mbf_xmax, number mbf_finf, CUDABox<number, number4> *box, number *_d_aff_eqdist, number *_d_aff_gamma, int *_d_affected_indx, int *_d_affected) {
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

        forces[IND] = (dF + F0);
        torques[IND] = (dT + T0);

        torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, torques[IND]);
    } else{

        number4 ftotal = make_number4<number, number4>(0,0,0,0);
        int pindex = IND - _offset;

        //Spring
        //bounds for spring loop determined per particle
        int &lb = _d_affected_indx[pindex];
        int &ub = _d_affected_indx[pindex+1];
        //printf("IND %d lb %d ub %d\n", IND, lb, ub);

        for(int i = lb; i < ub; i++) {
            int j = _d_affected[i];
            number4 &qpos = poss[j + _offset];
            number4 r = box->minimum_image(ppos, qpos);
            number d = _module<number, number4>(r);

            number &gamma = _d_aff_gamma[i];
            number &eqdist = _d_aff_eqdist[i];

            number fmod = (-1.0f * gamma) * (d - eqdist) / d;

            dF = fmod * r;
            dF.w = 0.5 * (0.5f * gamma * SQR(d - eqdist));

            ftotal -= dF;
        }
        LR_atomicAddXYZ(&(forces[IND]), ftotal);
        //printf("IND %d, f.x %.6f, f.y %.6f, f.z %.6f\n", IND, forces[IND].x, forces[IND].y, forces[IND].z);
    }
}

