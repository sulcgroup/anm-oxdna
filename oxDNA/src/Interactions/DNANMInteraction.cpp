/*
 * DNANMInteraction.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 *  Inherits from the DNA2Interaction Class
 *  Uses DNA2 Model and ACProtein Model
 */


#include "DNANMInteraction.h"
#include "DNA2Interaction.h"
#include "ACInteraction.h"
#include <sstream>
#include <fstream>

#include "../Particles/DNANucleotide.h"
#include "../Particles/ACParticle.h"

template<typename number>
void DNANMInteraction<number>::get_settings(input_file &inp){
	this->DNA2Interaction<number>::get_settings(inp);
	char parameterfile[500];
	getInputString(&inp, "PARFILE", parameterfile, 0);

	//Addition of Reading Parameter File
	int key1;
	int key2;
	char potswitch;
	double potential;
	pair <int, int> lkeys;
	pair <char, double> pot;
	double dist;
	string carbons;
	fstream parameters;
	parameters.open(parameterfile, ios::in);
	getline (parameters,carbons);
	if (parameters.is_open())
	{
		while (parameters.good())
		{
			parameters >> key1 >> key2 >> dist >> potswitch >> potential;
			lkeys.first=key1;
			lkeys.second=key2;
			pot.first=potswitch;
			pot.second=potential;
			_rknot[lkeys] = dist;
			_potential[lkeys]=pot;
		}
	parameters.close();
	}
	else
	{
		throw oxDNAException("ParameterFile Could Not Be Opened");
	}
}

template<typename number>
void DNANMInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N){
	this->DNA2Interaction<number>::check_input_sanity(particles,N);
}


template<typename number>
void DNANMInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < npro; i++) particles[i] = new ACParticle<number>();
	for(int i = npro; i < ndna; i++) particles[i] = new DNANucleotide<number>(_grooving);
}

template<typename number>
void DNANMInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;
	/*
	allocate_particles(particles, N); //How/Can I Use This?
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->strand_id = -1;
	}
	*/
	int my_N, my_N_strands;

	char line[2048];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if (!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting",
				this->_topology_filename);

	topology.getline(line, 2040);
	sscanf(line, "%d %d %d %d\n", &my_N, &my_N_strands,&ndna,&npro);

	char aminoacid[256];
	char base[256];
	int strand, i = 0;
	while (topology.good()) {
		topology.getline(line, 2040);
		if (strlen(line) == 0 || line[0] == '#')
			continue;
		if (i == N)
			throw oxDNAException(
					"Too many particles found in the topology file (should be %d), aborting",
					N);
		BaseParticle<number> *p = particles[i];
		std::stringstream ss(line);
		ss >> strand;

		if (strand<0){
			int nside, cside;
			ss >> aminoacid >> nside >> cside;
			if(strlen(base) == 1) {
				p->btype = Utils::decode_aa(base[0]);
			}
			int x;
			std::set<int> myneighs;
			while(ss.good())
			{
				ss >> x;
				if(x < 0 || x >= N)
				{
					throw oxDNAException(
					"Line %d of the topology file has an invalid syntax, neighbor has invalid id",i + 2);
				}
				myneighs.insert(x);
			}
			ACParticle<number> *p = dynamic_cast< ACParticle<number>  * > (particles[i]);

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
				  p->add_bonded_neighbor(dynamic_cast<ACParticle<number>  *> (particles[*k]) );
				}
			}
			p->strand_id = strand - 1;
			p->index = i;
			i++;

			// here we fill the affected vector
			if (p->n3 != P_VIRTUAL)
				p->affected.push_back(ParticlePair<number>(p->n3, p));
			if (p->n5 != P_VIRTUAL)
				p->affected.push_back(ParticlePair<number>(p, p->n5));

		} else {
			int tmpn3, tmpn5;
			ss>>base>>tmpn3>>tmpn5;

			BaseParticle<number> *p = particles[i];

			if(tmpn3 < 0) p->n3 = P_VIRTUAL;
			else p->n3 = particles[tmpn3];
			if(tmpn5 < 0) p->n5 = P_VIRTUAL;
			else p->n5 = particles[tmpn5];

			// store the strand id
			// for a design inconsistency, in the topology file
			// strand ids start from 1, not from 0
			p->strand_id = strand - 1;

			// the base can be either a char or an integer
			if(strlen(base) == 1) {
				p->type = Utils::decode_base(base[0]);
				p->btype = Utils::decode_base(base[0]);
			}
			else {
				if(atoi (base) > 0) p->type = atoi (base) % 4;
				else p->type = 3 - ((3 - atoi(base)) % 4);
				p->btype = atoi(base);
			}

			if(p->type == P_INVALID) throw oxDNAException ("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);
			p->index = i;
			i++;

			// here we fill the affected vector
			if (p->n3 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p->n3, p));
			if (p->n5 != P_VIRTUAL) p->affected.push_back(ParticlePair<number>(p, p->n5));
		}

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

};

template<typename number>
number DNANMInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (p->btype >= 0 && q->btype >=0){ //DNA-DNA Interaction
		number energy= this->DNA2Interaction<number>::pair_interaction_bonded(p,q,r,update_forces);
		return energy;
	} else if ((p->btype >= 0 && q->btype <0) ||(p->btype <0 && q->btype >=0)){ //DNA,Protein Interaction
		// Not Yet Supported
		number energy=0.0;
		return energy;
	} else if ((p->btype <0 && q->btype <0)){ //Protein,Protein Interaction
		number energy= ACInteraction<number>::pair_interaction_bonded(p,q,r,update_forces);
		return energy;
	} else{
		number energy=0.0;
		return energy;
	}
}

template<typename number>
number DNANMInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (p->btype >= 0 && q->btype >=0){ //DNA-DNA Interaction
		number energy= this->DNA2Interaction<number>::pair_interaction_nonbonded(p,q,r,update_forces);
		return energy;
	} else if ((p->btype >= 0 && q->btype <0) ||(p->btype <0 && q->btype >=0)){ //DNA,Protein Interaction
		number energy= _protein_dna_exc_volume(p,q,r,update_forces);
		return energy;
	} else if ((p->btype <0 && q->btype <0)){ //Protein,Protein Interaction
		number energy= ACInteraction<number>::pair_interaction_nonbonded(p,q,r,update_forces);
		return energy;
	} else{
		number energy=0.0;
		return energy;
	}
}



template<typename number>
number DNANMInteraction<number>::_protein_dna_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)
{
	 //printf("Calling crodwer DNA on %d %d \n",p->index,q->index);
	 BaseParticle<number> *protein;
	 BaseParticle<number> *nuc;

	 LR_vector<number> force(0, 0, 0);
     LR_vector<number> rcenter = *r;

	 if(p->btype >= 0 && q->btype < 0)
     {
		 //rcenter = -rcenter;
		 protein = q;
		 nuc = p;
	 }
	 else if (p->btype < 0 && q->btype >= 0 )
	 {
		 rcenter = -rcenter;
		 protein = p;
		 nuc = q;
	 }
	 else
	 {
		 return 0.f;
	 }

	 LR_vector<number> r_to_back = rcenter  - nuc->int_centers[DNANucleotide<number>::BACK];
	 LR_vector<number> r_to_base = rcenter  - nuc->int_centers[DNANucleotide<number>::BASE];

	 LR_vector<number> torquenuc(0,0,0);

	 number energy = this->_protein_dna_repulsive_lj(r_to_back, force, update_forces);
	 if (update_forces) {
		    torquenuc  -= nuc->int_centers[DNANucleotide<number>::BACK].cross(force);
	 		nuc->force -= force;
		 	protein->force += force;
	 }


	 energy += this->_protein_dna_repulsive_lj(r_to_base, force, update_forces);

	 if(update_forces) {

		    torquenuc  -= nuc->int_centers[DNANucleotide<number>::BASE].cross(force);
		    nuc->torque += nuc->orientationT * torquenuc;

		    //crowder->torque -= crowder->orientationT * torquenuc;
	 		nuc->force -= force;
	 		protein->force += force;
	 }


	 return energy;
}

template<typename number>
number DNANMInteraction<number>::_protein_dna_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(_pro_dna_rcut)) {
		if(rnorm > SQR(_pro_dna_rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - _pro_dna_rcut;
			energy = _pro_dna_stiffness * _pro_dna_b * SQR(rrc);
			if(update_forces) force = -r * (2 * _pro_dna_stiffness * _pro_dna_b * rrc / rmod);
		}
		else {
			number tmp = SQR(_pro_dna_sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * this->_pro_dna_stiffness * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * _pro_dna_stiffness * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}





template<typename number>
void DNANMInteraction<number>::init() {
	this->DNA2Interaction<number>::init();
//TODO: Figure out these values
	_pro_dna_sigma = 0.0786f;
	_pro_dna_rstar= 0.0746f;
	_pro_dna_b = 72471.9f;
	_pro_dna_rcut = 0.0798717f;
}

template class DNAInteraction<float>;
template class DNAInteraction<double>;
