/* System constants */
//added constants for DNACT
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

__constant__ float _kb;
__constant__ float _kt;

__constant__ int _npro;
__constant__ int _ndna;
__constant__ int _offset;

#include "../cuda_utils/CUDA_lr_common.cuh"
template<typename number, typename number4>
__forceinline__ __device__ void _ang_pot(number4 &dF, number4 &dT, const number4 &a1, const number4 &a3, const number4 &b1, const number4 &b3, const number4 &r, const number4 &ang_params, number (&kbkt)[2], int neighbor) {
    dF.x = dF.y = dF.z = dF.w = (number) 0.;
    dT.x = dT.y = dT.z = dT.w = (number) 0.;
    if(neighbor != 3 && neighbor != 5) return;
    number4 rij = make_number4<number, number4>(r.x/_module<number, number4>(r), r.y/_module<number>(r), r.z/_module<number>(r), 0.);
    number4 rji = (double) -1. * rij;

    number o1 = (CUDA_DOT(rij, a1)) - ang_params.x;
    number o2 = (CUDA_DOT(rji, b1)) - ang_params.y;
    number o3 = (CUDA_DOT(a1, b1)) - ang_params.z;
    number o4 = (CUDA_DOT(a3, b3)) - ang_params.w;
    number& kb = kbkt[0];
    number& kt = kbkt[1];
    number energy = kb / 2 * (SQR(o1) + SQR(o2)) + kt / 2 * (SQR(o3) + SQR(o4));

    if (neighbor == 5) {
        number4 F = -1*(((_cross<number, number4>(rji, _cross<number, number4>(rij, a1))) * kb * o1) -
                (_cross<number, number4>(rji, _cross<number, number4>(rij, b1))) * kb * o2) /
              _module<number, number4>(r);

        dF -= F;

        dF.w = 0.5*energy;

        number4 dTa = (_cross<number, number4>(rij, a1) * o1 *kb) - ((_cross<number, number4>(a1, b1) * o3 * kt) +
                                                                      (_cross<number, number4>(a3, b3) * o4 * kt));

        dT += dTa;

        //Debugging
        /*if(IND==104){
            printf("Neigh 5 p=%d\n", IND);
            printf("p %d q %d E= %.5f Angles o1 %.7f o2 %.7f o3 %.7f o4 %.7f\n", IND, IND +1, energy, o1, o2, o3, o4);
            printf("rij %.7f %.7f %.7f rji %.7f %.7f %.7f\n", rij.x, rij.y, rij.z, rji.x, rji.y, rji.z);
            printf("3 a1 %.7f %.7f %.7f b1 %.7f %.7f %.7f a3 %.7f %.7f %.7f b3 %.7f %.7f %.7f\n", a1.x, a1.y, a1.z,
                   b1.x, b1.y, b1.z, a3.x, a3.y, a3.z, b3.x, b3.y, b3.z);
            printf("force %.5f %.5f %.5f torque %.5f %.5f %.5f\n", dF.x, dF.y, dF.z, dT.x, dT.y, dT.z);
        }*/
    }

    if (neighbor == 3) {
        number4 F = -1*(((_cross<number, number4>(rji, _cross<number, number4>(rij, a1))) * kb * o1) -
                        (_cross<number, number4>(rji, _cross<number, number4>(rij, b1))) * kb * o2) /
                    _module<number, number4>(r);

        dF += F;

        dF.w = 0.5*energy;

        number4 dTb = (_cross<number, number4>(rji, b1) * o2 *kb) + ((_cross<number, number4>(a1, b1) * o3 * kt) +
                                                                      (_cross<number, number4>(a3, b3) * o4 * kt));

        dT += dTb;
        //Debugging
        /*if(IND==104){
            printf("Neigh 3 p=%d\n", IND);
            printf("p %d q %d E= %.5f Angles o1 %.7f o2 %.7f o3 %.7f o4 %.7f\n", IND, IND - 1, energy, o1, o2, o3, o4);
            printf("rij %.7f %.7f %.7f rji %.7f %.7f %.7f\n", rij.x, rij.y, rij.z, rji.x, rji.y, rji.z);
            printf("3 a1 %.7f %.7f %.7f b1 %.7f %.7f %.7f a3 %.7f %.7f %.7f b3 %.7f %.7f %.7f\n", a1.x, a1.y, a1.z,
                   b1.x, b1.y, b1.z, a3.x, a3.y, a3.z, b3.x, b3.y, b3.z);
            printf("force %.5f %.5f %.5f torque %.5f %.5f %.5f\n", dF.x, dF.y, dF.z, dT.x, dT.y, dT.z);
        }*/
    }

}

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
__global__ void dnact_forces_edge_nonbonded(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, edge_bond *edge_list, int n_edges, LR_bonds *bonds, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, CUDABox<number, number4> *box) {
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

        dT.x = -dT.x;
        dT.y = -dT.y;
        dT.z = -dT.z;

        //Add force to q index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f){
            LR_atomicAddXYZ(&(forces[to_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[to_index]), dT);
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

        dT.x = -dT.x;
        dT.y = -dT.y;
        dT.z = -dT.z;

        //add forces to p index
        if ((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (number) 0.f)
            LR_atomicAddXYZ(&(forces[from_index]), dF);
        if ((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (number) 0.f)
            LR_atomicAddXYZ(&(torques[from_index]), dT);

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
__global__ void dnact_forces_edge_bonded(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, LR_bonds *bonds, bool grooving, bool use_oxDNA2_FENE, bool use_mbf, number mbf_xmax, number mbf_finf, CUDABox<number, number4> *box, number *_d_aff_eqdist, number *_d_aff_gamma, number *_d_ang_params, number *_d_ang_kbkt, int *_d_affected_indx, int *_d_affected) {
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
        //Nucleotide
        LR_bonds bs = bonds[IND];
        // particle axes according to Allen's paper

        number4 a1, a2, a3;
        get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

        if(bs.n3 != P_INVALID) {
            number4 qpos = poss[bs.n3];
            number4 r = box->minimum_image(ppos, qpos);

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
    } else {
        //Protein Particles
        LR_bonds bs = bonds[IND];
//        printf("n3 %d p %d n5 %d\n", bs.n3, IND, bs.n5);
        int pindex = IND - _offset;

        int n3, n5 = -1;
        if (bs.n3 != -1) n3 = bs.n3 - _offset;
        if (bs.n5 != -1) n5 = bs.n5 - _offset;

        number4 p1, p2, p3;
        get_vectors_from_quat<number, number4>(orientations[IND], p1, p2, p3);
        //printf("p %d p1 %.4f %.4f %.4f p2 % .4f %.4f %.4f p3 %.4f %.4f %.4f\n", IND, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);

        number4 ftotal = make_number4<number, number4>(0, 0, 0, 0);
        number4 ttotal = make_number4<number, number4>(0, 0, 0, 0);

        number n3kbkt[2] = {_kb, _kt};
        number n5kbkt[2] = {_kb, _kt};

        if (_d_ang_kbkt != NULL) {
            n3kbkt[0] = _d_ang_kbkt[n3];
            n3kbkt[1] = _d_ang_kbkt[n3 + 1];
            n5kbkt[0] = _d_ang_kbkt[pindex];
            n5kbkt[1] = _d_ang_kbkt[pindex + 1];
        }

        //Previous neighbor
        if (n3 >= 0){
            number4 n3pos = poss[IND-1];
            number4 r = box->minimum_image(n3pos, ppos);
            number4 angle_o = make_number4<number, number4>(_d_ang_params[n3*4], _d_ang_params[n3*4+1], _d_ang_params[n3*4+2], _d_ang_params[n3*4+3]);
//            printf("p %d q %d Angle Params a0 %.4f b0 %.4f c0 %.4f d0 %.4f\n", pindex, IND-1, angle_o.x, angle_o.y, angle_o.z, angle_o.w);
            // In this case our particle is particle j of the pairwise interaction
            number4 &b1 = p1;
            number4 &b3 = p3;

            number4 a1, a2, a3;
            get_vectors_from_quat<number, number4>(orientations[IND-1], a1, a2, a3);
            //printf("p %d q %d a1 %.4f %.4f %.4f b1 %.4f %.4f %.4f\n", IND-1, IND, a1.x, a1.y, a1.z, b1.x, b1.y, b1.z);
            _ang_pot<number, number4>(dF, dT, a1, a3, b1, b3, r, angle_o, n3kbkt, 3);

            //Add contribution from Angular Potential
            //Energy is halved in each _ang_pot call
            //Torque is stored in dT
            ftotal += dF;
            ttotal += dT;

            //Excluded Volume Bonded n3 Neighbor (Since bonded neighs not included in edge list)
            _excluded_volume_quart(r, dF, MD_pro_sigma, MD_pro_rstar, MD_pro_b, MD_pro_rc);
            ftotal -= dF;
        }

//        if(IND == 2) printf("CUDA q2 Angular F %.4f %.4f %.4f T %.4f %.4f %.4f %.4f\n", dF.x, dF.y, dF.z, dT.x, dT.y, dT.z, dT.w);
        //Next neighbor
        if (n5 > 0) {
            number4 n5pos = poss[IND+1];
            number4 r = box->minimum_image(ppos, n5pos);
            number4 angle_o = make_number4<number, number4>(_d_ang_params[pindex*4], _d_ang_params[pindex*4+1], _d_ang_params[pindex*4+2], _d_ang_params[pindex*4+3]);
            // In this case our particle is particle i of the pairwise interaction
            number4 &a1 = p1;
            number4 &a3 = p3;

            number4 b1, b2, b3;
            get_vectors_from_quat<number, number4>(orientations[IND+1], b1, b2, b3);
            //printf("p %d q %d a1 %.4f %.4f %.4f b1 %.4f %.4f %.4f\n", IND, IND+1, a1.x, a1.y, a1.z, b1.x, b1.y, b1.z);
            _ang_pot<number, number4>(dF, dT, a1, a3, b1, b3, r, angle_o, n5kbkt, 5);

            //Add contribution from Angular Potential
            //Energy is halved in each _ang_pot call
            //Torque is stored in dT
            ftotal += dF;
            ttotal += dT;

            //Excluded Volume Bonded Neighbor n5 (Since bonded neighs not included in edge list)
            _excluded_volume_quart(r, dF, MD_pro_sigma, MD_pro_rstar, MD_pro_b, MD_pro_rc);
            ftotal += dF;
        }

        //Debugging
        //if(IND == 104) printf("CUDA p42 Angular F %.7f %.7f %.7f T %.7f %.7f %.7f %.7f\n", ftotal.x, ftotal.y, ftotal.z, ttotal.x, ttotal.y, ttotal.z, ttotal.w);

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

        //Add totals to particle lists
        LR_atomicAdd(&(forces[IND]), ftotal);
        LR_atomicAdd(&(torques[IND]), ttotal);

        //Necessary for Quaternion Calculations
        torques[IND] = _vectors_transpose_number4_product(p1, p2, p3, torques[IND]);

        //Debugging
        //if(IND == 104) printf("CUDA p2 Angular+Spring F %.7f %.7f %.7f T %.7f %.7f %.7f\n", forces[IND].x, forces[IND].y, forces[IND].z, torques[IND].x, torques[IND].y, torques[IND].z);
    }
}

