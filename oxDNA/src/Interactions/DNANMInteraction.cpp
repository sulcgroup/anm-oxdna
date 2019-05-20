/*
 * DNANMInteraction.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: jonah
 *  Inherits from the DNA2Interaction Class
 *  Uses DNA2 Model and ACProtein Model
 *  ACProtein Functions are implemented here as to avoid a multitude of multiple inheritance
 *  problems of which I am unsure are truly solvable
 */


#include "DNANMInteraction.h"
#include "DNA2Interaction.h"
#include "ACInteraction.h"
#include <sstream>
#include <fstream>

#include "../Particles/DNANucleotide.h"
#include "../Particles/ACParticle.h"

template<typename number>
DNANMInteraction<number>::DNANMInteraction() : DNA2Interaction<number>() { // @suppress("Class members should be properly initialized")
	//typedef std::map<int, number (child::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)> interaction_map;

	//Protein Methods Function Pointers
	this->_int_map[SPRING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNANMInteraction<number>::_protein_spring;
	this->_int_map[EXC_VOL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNANMInteraction<number>::_protein_exc_volume;

	//DNA Methods Function Pointers
	this->_int_map[this->BACKBONE] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNANMInteraction<number>::_backbone;
	this->_int_map[this->COAXIAL_STACKING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_coaxial_stacking;
	this->_int_map[this->CROSS_STACKING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_cross_stacking;
	this->_int_map[this->BONDED_EXCLUDED_VOLUME] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_bonded_excluded_volume;
	this->_int_map[this->STACKING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_stacking;
	this->_int_map[this->HYDROGEN_BONDING] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_hydrogen_bonding;
	this->_int_map[this->NONBONDED_EXCLUDED_VOLUME] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_nonbonded_excluded_volume;
	this->_int_map[this->DEBYE_HUCKEL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces))  &DNANMInteraction<number>::_debye_huckel;
}

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
	//this->DNA2Interaction<number>::check_input_sanity(particles,N);
	//Need to make own function that checks the input sanity
}

template<typename number>
void DNANMInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	if (ndna==0 || ndnas==0) throw oxDNAException("Missing DNA particles");
	for(int i = 0; i < ndna; i++) particles[i] = new DNANucleotide<number>(this->_grooving);
	for(int i = ndna; i < N; i++) particles[i] = new ACParticle<number>();
}

template<typename number>
void DNANMInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;
	int my_N, my_N_strands;

	char line[2048];
	std::ifstream topology;
	topology.open(this->_topology_filename, ios::in);

	if (!topology.good())
		throw oxDNAException("Can't read topology file '%s'. Aborting",
				this->_topology_filename);

	topology.getline(line, 2040);
	sscanf(line, "%d %d %d %d %d\n", &my_N, &my_N_strands, &ndna, &npairs, &ndnas);

	allocate_particles(particles, N); //How/Can I Use This?
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = 0;
	   particles[i]->strand_id = 0;
	}

	char aminoacid[256];
	char base[256];
	int strand, i = 0;
	while (topology.good()) {
		topology.getline(line, 2040);
		if (strlen(line) == 0 || line[0] == '#')
			continue;
		if (i == N+npairs)
			throw oxDNAException(
					"Too many particles found in the topology file (should be %d), aborting",
					N);
		//BaseParticle<number> *p = particles[i];
		std::stringstream ss(line);
		ss >> strand;

		if (strand<0){
			int nside, cside;
			ss >> aminoacid >> nside >> cside;

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

			if(strlen(base) == 1) {
				p->btype = Utils::decode_aa(base[0]);
			}
			if (nside < 0)
				p->n3 = P_VIRTUAL;
			else
				p->n3 = particles[nside];
			if (cside < 0)
				p->n5 = P_VIRTUAL;
			else
				p->n5 = particles[cside];

			p->strand_id = abs(strand)+ndnas-1;
			p->index = i;
			for(std::set<int>::iterator k = myneighs.begin(); k != myneighs.end(); ++k )
			{
				if(p->index < *k)
				{
				  p->add_bonded_neighbor(dynamic_cast<ACParticle<number>  *> (particles[*k]) );
				}
			}

			i++;
			// here we fill the affected vector
			if (p->n3 != P_VIRTUAL)
				p->affected.push_back(ParticlePair<number>(p->n3, p));
			if (p->n5 != P_VIRTUAL)
				p->affected.push_back(ParticlePair<number>(p, p->n5));

		} if(strand>0) {
			int tmpn3, tmpn5;
			ss>>base>>tmpn3>>tmpn5;

			DNANucleotide<number> *p = dynamic_cast< DNANucleotide<number>  * > (particles[i]);

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
		} if (strand == 0){
			int pro, dna;
			ss>>base>>dna>>pro; //Currently base is completely arbitrary but could be used in the future
			std::pair <int,int> pdpair;
			pdpair.first=dna;
			pdpair.second=pro;
			_ippairs.push_back(pdpair);
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
number DNANMInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces){
	number energy = pair_interaction_nonbonded(p, q, r, update_forces);
	energy += pair_interaction_bonded(p, q, r, update_forces);

	return energy;
}

template<typename number>
number DNANMInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	if (p->btype >= 0 && q->btype >=0){ //DNA-DNA Interaction
		number energy= this->DNA2Interaction<number>::pair_interaction_bonded(p,q,r,update_forces);
		return energy;
	} else if ((p->btype >= 0 && q->btype <0) ||(p->btype <0 && q->btype >=0)){ //DNA,Protein Interaction
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}

		std::pair<int,int> lookup;
		if (p->index>q->index){ //the higher index will always be the protein particle
			lookup.first=p->index;
			lookup.second=q->index;
		} else{
			lookup.second=p->index;
			lookup.first=q->index;
		}

		if ( std::find(this->_ippairs.begin(), this->_ippairs.end(), lookup) != this->_ippairs.end() ){
			number energy=_protein_spring(p,q,r,update_forces);
			return energy;
		}else{
			number energy=0.0;
			return energy;
		}
	} else if ((p->btype <0 && q->btype <0)){ //Protein,Protein Interaction
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}
		number energy= this->_protein_spring(p,q,r,update_forces);
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
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}
		number energy= this->_protein_dna_exc_volume(p,q,r,update_forces);
		return energy;
	} else if ((p->btype <0 && q->btype <0)){ //Protein,Protein Interaction
		LR_vector<number> computed_r(0, 0, 0);
		if(r == NULL) {
			if (q != P_VIRTUAL && p != P_VIRTUAL) {
				computed_r = q->pos - p->pos;
				r = &computed_r;
			}
		}
		number energy= this->_protein_exc_volume(p,q,r,update_forces);
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

	 number energy = this->_protein_dna_repulsive_lj(r_to_back, force, update_forces, _pro_backbone_sigma, _pro_backbone_b, _pro_backbone_rstar,_pro_backbone_rcut,_pro_backbone_stiffness);
	 if (update_forces) {
		    torquenuc  -= nuc->int_centers[DNANucleotide<number>::BACK].cross(force);
	 		nuc->force -= force;
		 	protein->force += force;
	 }


	 energy += this->_protein_dna_repulsive_lj(r_to_base, force, update_forces, _pro_base_sigma, _pro_base_b, _pro_base_rstar, _pro_base_rcut, _pro_base_stiffness);

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
number DNANMInteraction<number>::_protein_dna_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rcut)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rcut;
			energy = stiffness * b * SQR(rrc);
			if(update_forces) force = -r * (2 * stiffness * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * stiffness * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * stiffness * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

template<typename number>
void DNANMInteraction<number>::init() {
	this->DNA2Interaction<number>::init();
//TODO: Figure out these values
    //Backbone-Protein Excluded Volume Parameters
	_pro_backbone_sigma = 1.05;
	_pro_backbone_rstar= 0.97285f;
	_pro_backbone_b = 483.718f;
	_pro_backbone_rcut = 1.05998f;
    _pro_backbone_stiffness = 1.0f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.68f;
	_pro_base_rstar= 0.62452f;
	_pro_base_b = 1256.04f;
	_pro_base_rcut = 0.683987f;
    _pro_base_stiffness = 1.0f;
    //Protein-Protein Excluded Volume Parameters
	_pro_sigma = 0.3480514f;
	_pro_rstar= 0.348f;
	_pro_b = 336800.0f;
	_pro_rcut = 0.348103f;
	ndna=0;
	npairs=0;
	ndnas=0;
}

//Functions from ACInteraction.h
//Stolen due to inheritance issues

template<typename number>
number DNANMInteraction<number>::_protein_repulsive_lj(const LR_vector<number> &r, LR_vector<number> &force, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(_pro_rcut)) {
		if(rnorm > SQR(_pro_rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - _pro_rcut;
			energy = EXCL_EPS * _pro_b * SQR(rrc);
			if(update_forces) force = -r * (2 * EXCL_EPS * _pro_b * rrc/ rmod);
		}
		else {
			number tmp = SQR(_pro_sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r* (24 * EXCL_EPS * (lj_part - 2*SQR(lj_part))/rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

template<typename number>
number DNANMInteraction<number>::_protein_exc_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (p->index != q->index){
		LR_vector<number> force(0,0,0);

		number energy =  DNANMInteraction<number>::_protein_repulsive_lj(*r, force, update_forces);

		if(update_forces)
		{
			p->force -= force;
			q->force += force;
		}

		return energy;
	} else {
		return (number) 0.f;
	}
}

template<typename number>
number DNANMInteraction<number>::_protein_spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	pair <int,int> keys;
	number eqdist;
	char interactiontype;
	if (p->index != q->index)
	{
		if (p->index > q->index)
		{
			keys.first=q->index;
			keys.second=p->index;
		} else {
			keys.first=p->index;
			keys.second=q->index;
		}
		eqdist = _rknot[keys];
		interactiontype = _potential[keys].first;
		if (eqdist != 0.0){ //only returns number if eqdist is in .par file
			switch (interactiontype){
				case 's':
					{
						//Harmonic Spring Potential
						if ((eqdist < 0.0) || (eqdist > 2.0))  //ensures r0 is less than 7 Angstrom cutoff and nonnegative
						{
							if (keys.first+1 != keys.second){
								throw oxDNAException("No rknot or invalid rknot value for particle %d and %d rknot was %f", q->index, p->index, eqdist);
							}
						}
						number _k = _potential[keys].second; //stiffness of the spring
						if ((_k == 0) || (_k < 0)){
							throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
						}
						number rnorm = r->norm();
						number rinsta = sqrt(rnorm);
						number energy = 0.5 * _k * SQR(rinsta-eqdist);

						if (update_forces)
						{
							LR_vector<number> force(*r );
							force *= (-1.0f * _k ) * (rinsta-eqdist)/rinsta;
							p->force -= force;
							q->force += force;
						//printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
						//printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
						}
						return energy;
					} break;
				case 'i':
					{
						//Every Possible Pair of Particles Needs to be Calculated
						number _k = _potential[keys].second; //stiffness of the spring
						if ((_k == 0) || (_k < 0)){
							throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
						}
						number rnorm = r->norm();
						number rinsta = sqrt(rnorm);
						number energy = 0.5 * _k * SQR(rinsta-eqdist)*(1/eqdist);

						if (update_forces)
						{
							LR_vector<number> force(*r );
							force *= (-1.0f * _k ) * ((rinsta-eqdist)/rinsta) * (1/eqdist);
							p->force -= force;
							q->force += force;
						//printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
						//printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
						}
						return energy;
					} break;
				case 'e':
					{
						//Every Possible Pair of Particles Needs to be Calculated
						number _k = _potential[keys].second; //stiffness of the spring
						if ((_k == 0) || (_k < 0)){
							throw oxDNAException("No Spring Constant or invalid Spring Constant for particle %d and %d spring constant was %f", p->index, q->index, _k);
						}
						number rnorm = r->norm();
						number rinsta = sqrt(rnorm);
						number energy = 0.5 * _k * SQR(rinsta-eqdist)*(1/(eqdist*eqdist));

						if (update_forces)
						{
							LR_vector<number> force(*r );
							force *= (-1.0f * _k ) * ((rinsta-eqdist)/rinsta) * (1/(eqdist*eqdist));
							p->force -= force;
							q->force += force;
						//printf("@@@: particle %d and %d rinsta=%f , eqdist=%f, r-r0 = %f, prefactor = %f, force = %f,%f,%f, ener=%f \n",p->index,q->index,rinsta,eqdist, rinsta-eqdist, (-1.0f * _k ) * (rinsta-eqdist)/rinsta, force.x,force.y,force.z,energy);
						//printf("@@@: %f %f \n",rinsta,(-1.0f * _k ) * (rinsta-eqdist)/rinsta);
						}
						return energy;
					} break;
				default:
					{
						throw oxDNAException("Interaction type specified in .par file is Invalid, particles %d and %d, switch %c",p->index,q->index,interactiontype);
					}
			}
		} else {
			return (number) 0.f; //returns 0 if no rknot value in parameter value aka they aren't bonded
		}
	} else {
		return (number) 0.f; //returns 0 if particle pair consists of particle and itself
	}

}

template<typename number>
DNANMInteraction<number>::~DNANMInteraction() {

}





