/*
 * GCInteraction.cpp
 *
 *  Created on: 10/feb/2018
 *      Author: jonah
 */

#include "GCInteraction.h"
#include "../Particles/GCParticle.h"
#include <sstream>

template<typename number>
GCInteraction<number>::GCInteraction() : BaseInteraction<number, GCInteraction<number> >() {
	this->_int_map[SPRING_POTENTIAL] = &GCInteraction<number>::_spring;
	// this->_int_map[EXC_VOL] = &GCInteraction<number>::_exc_volume;
}

template<typename number>
GCInteraction<number>::~GCInteraction() {

}

template<typename number>
void GCInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
}

template<typename number>
void GCInteraction<number>::init() {

//TODO: Figure out these values
	_r = 0.07865696f;
	_k = 0.75f;
	_sigma = 0.06755f;
	_rstar= 0.0787f;
	_b = -155.35f;
	_rcut = 0.1573f;
}

template<typename number>
void GCInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new GCParticle<number>();
}

template<typename number>
void GCInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
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

		//int res = sscanf(line, "%d %s %d %d", &strand, base, &tmpn3, &tmpn5);
		std::stringstream ss(line);
		ss >> strand >> aminoacid >> nside >> cside;
		if(!ss.good())
		{
			throw oxDNAException(
								"Line %d of the topology file has an invalid syntax",
								i + 2);
		}

		int x;
		std::set<int> myneighs;
		ss >> x;
		while(ss.good())
		{
			if(x < 0 || x >= N)
			{
				throw oxDNAException(
									"Line %d of the topology file has an invalid syntax, neighbor has invalid id",
									i + 2);
			}
			myneighs.insert(x);
			ss >> x;
		}


		GCParticle<number> *p = dynamic_cast< GCParticle<number>  * > (particles[i]);

		if (nside < 0)
			p->n3 = P_VIRTUAL;
		else
			p->n3 = particles[nside];
		if (cside < 0)
			p->n5 = P_VIRTUAL;
		else
			p->n5 = particles[cside];

		for(std::set<int>::iterator k = myneighs.begin(); k != myneighs.end(); ++k )
		{
			if(p->index < *k)
			{
			  p->add_bonded_neighbor(dynamic_cast<GCParticle<number>  *> (particles[*k]) );
			}
		}

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = strand - 1;


		//TODO WHAT IS IT DOING WITH THE ATOI stuff????
		//THIS is from DNAInteraction.cpp... How much of it do I really need?
		// the base can be either a char or an integer
		if (strlen(aminoacid) == 1) {
			p->type = Utils::decode_aminoacid(aminoacid[0]);
			p->btype = Utils::decode_aminoacid(aminoacid[0]);
		} else {
			if (atoi(aminoacid) > 0)
				p->type = atoi(aminoacid) % 4;
			else
				p->type = 3 - ((3 - atoi(aminoacid)) % 4);
			p->btype = atoi(aminoacid);
		}

		// until here replacem after here is fine


		if (p->type == P_INVALID)
			throw oxDNAException(
					"Particle #%d in strand #%d contains a non valid aminoacid '%c'. Aborting",
					i, strand, aminoacid);

		p->index = i;
		i++;

		// here we fill the affected vector
		if (p->n3 != P_VIRTUAL)
			p->affected.push_back(ParticlePair<number>(p->n3, p));
		if (p->n5 != P_VIRTUAL)
			p->affected.push_back(ParticlePair<number>(p, p->n5));
	}
	// TODO: Is this the right N?
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
number GCInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_nonbonded(p, q, r, update_forces);
	energy += pair_interaction_bonded(p, q, r, update_forces);

	return energy;
}

template<typename number>
number GCInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			computed_r = this->_box->min_image(p->pos, q->pos);
			r = &computed_r;
	}

	return (number) this->_spring(p,q,r,update_forces);
}

template<typename number>
number GCInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	/*LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _exc_volume(p, q, r, update_forces);*/
	return (number) 0.f;
}

template<typename number>
void GCInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class GCInteraction<float>;
template class GCInteraction<double>;
