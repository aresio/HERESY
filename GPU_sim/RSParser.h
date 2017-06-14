#ifndef __RSP__
#define __RSP__

#include <string>
#include <vector>
#include <map>

const unsigned int FINE_REAGENTI = UINT_MAX;
const unsigned int FINE_INIBITORI = UINT_MAX-1;
const unsigned int FINE_REAZIONE = UINT_MAX-2;

class ReactionSystemsParser {

public:
	bool OpenFile(std::string path, std::string path_m0, std::string context);
	std::vector<unsigned int> vettore_regole;
	std::vector<unsigned int> vettore_offset;
	std::vector<char> vettore_stati;
	std::vector<char> vector_context;
	unsigned int iterations;

	ReactionSystemsParser() { this->iterations =0; };
	
	std::map<std::string, unsigned int> insieme_specie; 
	std::map<unsigned int, std::string> rev_insieme_specie;

	unsigned int get_number_of_species() { return this->insieme_specie.size(); };
	unsigned int get_number_of_reactions() { return this->vettore_offset.size(); };
	unsigned int get_iterations() { return this->iterations; };

};

#endif