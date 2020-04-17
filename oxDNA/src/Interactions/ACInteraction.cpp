/*
 * ACInteraction.cpp
 *
 *  Created on: 10/feb/2018
 *      Author: jonah
 */

#include "ACInteraction.h"
#include "../Particles/ACParticle.h"
#include <sstream>
#include <fstream>
#include <unistd.h>


// TODO: Get the files loaded with their strings saved in the init


template<typename number>
ACInteraction<number>::ACInteraction(){
	this->_int_map[SPRING_POTENTIAL] = &ACInteraction<number>::_spring;
	this->_int_map[EXC_VOL] = &ACInteraction<number>::_exc_volume;

}

template<typename number>
ACInteraction<number>::~ACInteraction() {

}

template<typename number>
void ACInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char parameterfile[500];
	getInputString(&inp, "PARFILE", parameterfile, 0);

	//Addition of Reading Parameter File for ACInteraction Only!
	int key1, key2;
	char potswitch;
	double potential, dist;
	string carbons;
	fstream parameters;
	parameters.open(parameterfile, ios::in);
	getline (parameters,carbons);
	if (parameters.is_open())
	{
		while (parameters.good())
		{
			parameters >> key1 >> key2 >> dist >> potswitch >> potential;
            pair <int, int> lkeys (key1, key2);
            pair <char, double> pot (potswitch, potential);
			_rknot[lkeys] = dist;
			_potential[lkeys] = pot;
		}
	parameters.close();
	}
	else
	{
		throw oxDNAException("ParameterFile Could Not Be Opened");
	}
}

template<typename number>
void ACInteraction<number>::init() {
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
    _sigma = 0.35f;
    _rstar = 0.349f;
    _b = 306484596.421f;
    _rc = 0.352894;
}

template<typename number>
void ACInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new ACParticle<number>();
}

template<typename number>
void ACInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
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
		ACParticle<number> *p = dynamic_cast< ACParticle<number>  * > (particles[i]);

		if(strlen(aminoacid) == 1) {
		    p->btype = Utils::decode_aa(aminoacid[0]);
		}


		if (nside < 0)
			p->n3 = P_VIRTUAL;
		else
            p->add_bonded_neighbor(dynamic_cast<ACParticle<number>  *> (particles[nside]) );
		if (cside < 0)
			p->n5 = P_VIRTUAL;
		else
			p->add_bonded_neighbor(dynamic_cast<ACParticle<number>  *> (particles[cside]) );

		for(std::set<int>::iterator k = myneighs.begin(); k != myneighs.end(); ++k )
		{
			if(p->index < *k)
			{
			    p->add_bonded_neighbor(dynamic_cast<ACParticle<number>  *> (particles[*k]) );
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
number ACInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_nonbonded(p, q, r, update_forces);
	energy += pair_interaction_bonded(p, q, r, update_forces);

	return energy;
}

template<typename number>
number ACInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

    ACParticle<number> *cp = dynamic_cast< ACParticle<number> * > (p);
    if ((*cp).ACParticle<number>::is_bonded(q)){
        number energy = _spring(p,q,r,update_forces);
        energy += _exc_volume(p,q,r,update_forces);
        return energy;
    } else {
        return 0.f;
    }
}

template<typename number>
number ACInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

    ACParticle<number> *cp = dynamic_cast< ACParticle<number> * > (p);
    if ((*cp).ACParticle<number>::is_bonded(q)){
        return 0.f;
    } else {
        number energy = this->_exc_volume(p, q, r, update_forces);
        return energy;
    }
}

template<typename number>
void ACInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {
}

template class ACInteraction<float>;
template class ACInteraction<double>;
