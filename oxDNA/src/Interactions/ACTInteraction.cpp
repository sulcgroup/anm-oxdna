/*
 * ACTInteraction.cpp
 *
 *  Created on: 16/July/2020
 *      Author: jonah
 */

#include "ACTInteraction.h"
#include "../Particles/ACTParticle.h"
#include <sstream>
#include <fstream>
#include <unistd.h>


template<typename number>
ACTInteraction<number>::ACTInteraction(){
    this->_int_map[ANGULAR_POTENTIAL] = (number (ACInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &ACTInteraction<number>::_ang_pot;
    this->_int_map[this->SPRING_POTENTIAL] = (number (ACInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &ACTInteraction<number>::_spring;
    this->_int_map[this->EXC_VOL] = (number (ACInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &ACTInteraction<number>::_exc_volume;
}

template<typename number>
ACTInteraction<number>::~ACTInteraction() = default;

template<typename number>
void ACTInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char parameterfile[500];
	float kb_tmp, kt_tmp;
	getInputString(&inp, "parfile", parameterfile, 1);
	getInputFloat(&inp, "bending_k", &kb_tmp, 0);
	_k_bend = (number) kb_tmp;

	getInputFloat(&inp, "torsion_k", &kt_tmp, 0);
	_k_tor = (number) kt_tmp;

	if(_k_bend < 0 || _k_tor < 0) throw oxDNAException("Invalid Bending or Torsion Constants");

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
	int key1, key2;
	char potswitch;
	double potential, dist;
	double a0, b0, c0, d0;
	string carbons;
	fstream parameters;
	parameters.open(parameterfile, ios::in);
	getline (parameters,carbons);
//	printf("Carbons %s", carbons.c_str());
	int N = stoi(carbons);
	if (parameters.is_open())
	{
		while (parameters.good())
		{
			parameters >> key1 >> key2 >> dist >> potswitch >> potential;
			valid_spring_params(N, key1, key2, dist, potswitch, potential);
			if (key2 - key1 == 1){
			    //Angular Parameters
			    parameters >> a0 >> b0 >> c0 >> d0;
			    valid_angles(a0, b0, c0, d0);
			    vector<double> angles {a0, b0, c0, d0};
                _ang_vals[key1] = angles;

                pair <int, int> lkeys (key1, key2);
                pair <char, double> pot (potswitch, potential);

                this->_rknot[lkeys] = dist;
                this->_potential[lkeys] = pot;
			} else {
                pair <int, int> lkeys (key1, key2);
                pair <char, double> pot (potswitch, potential);
                this->_rknot[lkeys] = dist;
                this->_potential[lkeys] = pot;
            }
		}
	parameters.close();
	}
	else
	{
		throw oxDNAException("ParameterFile Could Not Be Opened");
	}
}

template<typename number>
void ACTInteraction<number>::init() {
    //initialize Excluded Volume Parameters
    ACInteraction<number>::init();
}

template<typename number>
void ACTInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new ACTParticle<number>();
}

template<typename number>
void ACTInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->strand_id = -1;
	}

	int my_N, my_N_strands;

	char line[5120];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if (!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting",
				this->_topology_filename);

	topology.getline(line, 5120);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	char aminoacid[256];
	int strand, i = 0;
	while (topology.good()) {
		topology.getline(line, 2040);
		if (strlen(line) == 0 || line[0] == '#')
			continue;
		if (i == N)
			throw oxDNAException(
					"Too many particles found in the topology file (should be %d), aborting",
					N);

		int nside, cside;
		std::stringstream ss(line);
		ss >> strand >> aminoacid >> nside >> cside;

		int x;
		std::set<int> myneighs;
		while(ss.good())
		{
			ss >> x;
			if(x < 0 || x >= N)
			{
				throw oxDNAException("Line %d of the topology file has an invalid syntax, neighbor has invalid id",
									i + 2);
			}

			myneighs.insert(x);
		}
		ACTParticle<number> *p = dynamic_cast< ACTParticle<number>  * > (particles[i]);

		if(strlen(aminoacid) == 1) {
		    p->btype = Utils::decode_aa(aminoacid[0]);
		}


		if (nside < 0)
			p->n3 = P_VIRTUAL;
		else
            p->add_bonded_neighbor(dynamic_cast<ACTParticle<number>  *> (particles[nside]) );
		if (cside < 0)
			p->n5 = P_VIRTUAL;
		else
			p->add_bonded_neighbor(dynamic_cast<ACTParticle<number>  *> (particles[cside]) );

		for(std::set<int>::iterator k = myneighs.begin(); k != myneighs.end(); ++k )
		{
			if(p->index < *k)
			{
			    p->add_bonded_neighbor(dynamic_cast<ACTParticle<number>  *> (particles[*k]) );
			}
		}

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = abs(strand) - 1;
		p->index = i;
		i++;

	}

	if (i < my_N)
		throw oxDNAException(
				"Not enough particles found in the topology file (should be %d). Aborting",
				my_N);

	topology.close();

	if (my_N != N)
		throw oxDNAException(
				"Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");

	*N_strands = my_N_strands;


}

template<typename number>
number ACTInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_nonbonded(p, q, r, update_forces);
	energy += pair_interaction_bonded(p, q, r, update_forces);

	return energy;
}

template<typename number>
number ACTInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

    ACTParticle<number> *cp = dynamic_cast< ACTParticle<number> * > (p);
    if ((*cp).ACTParticle<number>::is_bonded(q)){
        number energy = this->_spring(p,q,r,update_forces);
        energy += this->_exc_volume(p,q,r,update_forces);
        if (q->index - p->index == 1) energy += _ang_pot(p, q, r, update_forces);
        return energy;
    } else {
        return 0.f;
    }
}

template<typename number>
number ACTInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

    ACTParticle<number> *cp = dynamic_cast< ACTParticle<number> * > (p);
    if ((*cp).ACTParticle<number>::is_bonded(q)){
        return 0.f;
    } else {
        number energy = this->_exc_volume(p, q, r, update_forces);
        return energy;
    }
}

template<typename number>
void ACTInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
}

template<typename number>
number ACTInteraction<number>::_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
    number eqdist;
    pair <int,int> keys (std::min(p->index, q->index), std::max(p->index, q->index));
    eqdist = this->_rknot[keys];
    //Harmonic Spring Potential
    number _k = this->_potential[keys].second; //stiffness of the spring
    number rnorm = r->norm();
    number rinsta = sqrt(rnorm);
    number energy = 0.5 * _k * SQR(rinsta-eqdist);

    if (update_forces) {
        LR_vector<number> force(*r );
        force *= (-1.0f * _k ) * (rinsta-eqdist)/rinsta;

        p->force -= force;
        q->force += force;
    }

    return energy;
}

template<typename number>
number ACTInteraction<number>::_ang_pot(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
    // Get Angular Parameters
    vector<double> &ang_params = _ang_vals[p->index];
    double &a0 = ang_params[0];
    double &b0 = ang_params[1];
    double &c0 = ang_params[2];
    double &d0 = ang_params[3];

    LR_vector<number> rij_unit = *r;
    rij_unit.normalize();
    LR_vector<number> rji_unit = rij_unit*-1.f;

    LR_vector<number> &a1 = p->orientationT.v1;
    LR_vector<number> &b1 = q->orientationT.v1;
    LR_vector<number> &a3 = p->orientationT.v3;
    LR_vector<number> &b3 = q->orientationT.v3;
//    printf("i %d j %d a1 %.4f %.4f %.4f\n", p->index, q->index, a1.x, a1.y, a1.z);
//    printf("i %d j %d a3 %.4f %.4f %.4f\n", p->index, q->index, a3.x, a3.y, a3.z);
//    printf("i %d j %d b1 %.4f %.4f %.4f\n", p->index, q->index, b1.x, b1.y, b1.z);
//    printf("i %d j %d b3 %.4f %.4f %.4f\n", p->index, q->index, b3.x, b3.y, b3.z);


    double o1 = rij_unit*a1-a0;
    double o2 = rji_unit*b1-b0;
    double o3 = a1*b1-c0;
    double o4 = a3*b3-d0;

//    double o1 = rij_unit*a1;
//    double o2 = rji_unit*b1;
//    double o3 = a1*b1;
//    double o4 = a3*b3;
//    printf("i %d j %d dANg = %.5f, %.5f, %.5f, %.5f \n", p->index, q->index, o1, o2, o3, o4);


    //Torsion and Bending
    number energy = _k_bend/2 * (SQR(o1) + SQR(o2)) + _k_tor/2 * (SQR(o3) + SQR(o4));
//    printf("Angular i %d j %d E=%.5f \n", p->index, q->index, energy);

    if(update_forces){

        LR_vector<number> force = -(rji_unit.cross(rij_unit.cross(a1))*_k_bend*o1 - rji_unit.cross(rij_unit.cross(b1))*_k_bend*o2)/r->module();

        p->force -= force;
        q->force += force;

        LR_vector<number> ta = rij_unit.cross(a1) * o1 * _k_bend;
        LR_vector<number> tb = rji_unit.cross(b1) * o2 * _k_bend;

        p->torque += p->orientationT*ta;
        q->torque += q->orientationT*tb;

        LR_vector<number> torsion = a1.cross(b1) * o3 * _k_tor + a3.cross(b3) * o4 * _k_tor;

        p->torque += p->orientationT*-torsion;
        q->torque += q->orientationT*torsion;
    }

    return energy;
}

////WORKING FULL VERSION -debug statements
//template<typename number>
//number ACTInteraction<number>::_ang_pot(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
//    // Get Angular Parameters
//    vector<double> &ang_params = _ang_vals[p->index];
//    double &a0 = ang_params[0];
//    double &b0 = ang_params[1];
//    double &c0 = ang_params[2];
//    double &d0 = ang_params[3];
//
//    LR_vector<number> rij_unit = *r;
//    rij_unit.normalize();
//    LR_vector<number> rji_unit = rij_unit * -1.f;
//
//    LR_vector<number> &a1 = p->orientationT.v1;
//    LR_vector<number> &b1 = q->orientationT.v1;
//    LR_vector<number> &a3 = p->orientationT.v3;
//    LR_vector<number> &b3 = q->orientationT.v3;
//
//
//    double o1 = rij_unit * a1 - a0;
//    double o2 = rji_unit * b1 - b0;
//    double o3 = a1 * b1 - c0;
//    double o4 = a3 * b3 - d0;
//
//    //Torsion and Bending
//    number energy = _k_bend / 2 * (SQR(o1) + SQR(o2)) + _k_tor / 2 * (SQR(o3) + SQR(o4));
//
//    if (update_forces) {
//
//        LR_vector<number> force = -(rji_unit.cross(rij_unit.cross(a1)) * _k_bend * o1 -
//                                    rji_unit.cross(rij_unit.cross(b1)) * _k_bend * o2) / r->module();
//
//        p->force -= force;
//        q->force += force;
//
//        LR_vector<number> ta = rij_unit.cross(a1) * o1 * _k_bend;
//        LR_vector<number> tb = rji_unit.cross(b1) * o2 * _k_bend;
//
//        p->torque += p->orientationT * ta;
//        q->torque += q->orientationT * tb;
//
//        LR_vector<number> torsion = a1.cross(b1) * o3 * _k_tor + a3.cross(b3) * o4 * _k_tor;
//
//        p->torque += p->orientationT * -torsion;
//        q->torque += q->orientationT * torsion;
//    }
//
//    return energy;
//}









template class ACTInteraction<float>;
template class ACTInteraction<double>;
