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
	this->_int_map[this->SPRING_POTENTIAL] = &ACTInteraction<number>::_spring;
	this->_int_map[this->EXC_VOL] = &ACTInteraction<number>::_exc_volume;
    this->_int_map[ANGULAR_POTENTIAL] = &ACTInteraction<number>::_exc_volume; //change function to Angular Potential
}

template<typename number>
ACTInteraction<number>::~ACTInteraction() = default;

template<typename number>
void ACTInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char parameterfile[500];
	getInputString(&inp, "PARFILE", parameterfile, 0);

	//Addition of Reading Parameter File
	int key1, key2;
	char potswitch;
	double potential, dist;
	double a0, b0, c0;
	string carbons;
	fstream parameters;
	parameters.open(parameterfile, ios::in);
	getline (parameters,carbons);
	if (parameters.is_open())
	{
		while (parameters.good())
		{
			parameters >> key1 >> key2 >> dist >> potswitch >> potential;
			if (key2 - key1 == 1){
			    parameters >> a0 >> b0 >> c0;
			    tuple <double, double, double> angles = make_tuple(a0, b0, c0);
                pair <int, int> lkeys (key1, key2);
                pair <char, double> pot (potswitch, potential);
                _ang_vals[key1] = angles;
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
//    _sigma = 0.117f;
//    _rstar= 0.087f;
//    _b = 671492.f;
//    _rc = 0.100161f;

//    _sigma = 0.10f;
//    _rstar= 0.08f;
//    _b = 163352410.227f;
//    _rc = 0.102644f;

    //quartic parameters
//    _sigma = 0.3f;
//    _rstar= 0.23f;
//    _b = 3614904.f;
//    _rc = 0.298003f;

    //large possible parameters
//    _sigma = 0.55f;
//    _rstar = 0.47f;
//    _b = 80892.1;
//    _rc = 0.588787;

    //oof not those lets try these
    this->_sigma = 0.35f;
    this->_rstar = 0.349f;
    this->_b = 306484596.421f;
    this->_rc = 0.352894;

    //Angular Constants
    _sigma_bond_sqr = 0.5;
    _sigma_tor_sqr = 0.5;
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

	char line[2048];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if (!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting",
				this->_topology_filename);

	topology.getline(line, 2040);

	sscanf(line, "%d %d\n", &my_N, &my_N_strands);

	char aminoacid[2040];
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
number ACTInteraction<number>::_ang_pot(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
    // Get Angular Parameters
    tuple<double, double, double> ang_params = _ang_vals[p->index];
    double a0 = get<0>(ang_params);
    double b0 = get<1>(ang_params);
    double c0 = get<2>(ang_params);

    LR_vector<number> &r_unit = *r;
    r_unit.normalize();
    LR_vector<number> &a1 = p->orientationT.v1;
    LR_vector<number> &b1 = q->orientationT.v1;

    double o1 = r_unit*a1-a0;
    double o2 = -r_unit*b1-b0;
    double o3 = a1*b1-c0;


    number energy = exp(SQR(o1)/(2*_sigma_bond_sqr) + SQR(o2)/(2*_sigma_bond_sqr) + SQR(o3)/(2*_sigma_tor_sqr));

    if(update_forces){
        LR_vector<number> force = -energy * ((a1*(o1/_sigma_bond_sqr)) - b1*(o2/_sigma_bond_sqr));
        p->force -= force;
        q->force += force;

        LR_vector<number> torq_piece = (o3/_sigma_tor_sqr)*a1.cross(b1);
        p->torque += energy*(o1/_sigma_bond_sqr * r_unit.cross(a1) + torq_piece);
        q->torque += energy*(o2/_sigma_bond_sqr * r_unit.cross(b1) - torq_piece);
    }

    return energy;
}



template class ACTInteraction<float>;
template class ACTInteraction<double>;