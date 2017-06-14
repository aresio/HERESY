#include "RSParser.h"
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <regex>
#include <vector>
#include <istream>
#include <sstream>
#include <math.h>
#include <climits>

void tokenize( std::string s, std::vector<std::string>& v )  {
    std::istringstream buf(s.c_str());
    for(std::string token; getline(buf, token, '\t'); )
        v.push_back(token);
}


bool ReactionSystemsParser::OpenFile(std::string path, std::string path_m0, std::string path_context) {

	// printf(" * Opening Reaction System in file '%s'... ", path.c_str());

	unsigned int species_index  = 0;

	std::ifstream myfile(path.c_str());

	if ( myfile.is_open() ) {

		//printf("parsing for chemical species... ");

		// step 1: does line begin with substring "  Reaction "?
		std::string linea;
		std::string sublinea;
		std::smatch mr; 
		std::smatch mf; 
		std::regex e ( "\\\[.*?\\\]" );
		std::regex f ( "\\\".*?\\\"" );
		// std::regex e ("Reaction");
		

		while( myfile.good() ) {
			
			getline(myfile, linea);

			// gotcha: step 2, subdivide in [...] tokens
			if (linea.compare(0, 11, "  Reaction ")==0) {

				vettore_offset.push_back( vettore_regole.size() );
		
				// substep 1: parsa reagenti
				std::regex_search(linea, mr, e);
				// printf(" * Reactants: ");
				// printf("%s\n", mr[0].str().c_str());			

					// subsubstep 1: snocciola reagenti	
					sublinea = mr[0].str();					
					while( std::regex_search(sublinea, mf, f) ) {
						// printf("%s\n", mf[0].str().c_str());			

						if (insieme_specie.count(mf[0].str())==0) {
							insieme_specie[mf[0].str()] = species_index;
							rev_insieme_specie[species_index]=mf[0].str();
							species_index++;
						}
						vettore_regole.push_back( insieme_specie[mf[0].str()] );

						sublinea = mf.suffix().str();
					}

				linea = mr.suffix().str();

				vettore_regole.push_back( FINE_REAGENTI );
				
				// substep 2: parsa inibitori
				std::regex_search(linea, mr, e);
				// printf(" * Inhibitors: ");
				// printf("%s\n", mr[0].str().c_str());	

					// subsubstep 2: snocciola reagenti	
					sublinea = mr[0].str();					
					while( std::regex_search(sublinea, mf, f) ) {
						// printf("%s\n", mf[0].str().c_str());			

						if (insieme_specie.count(mf[0].str())==0) {
							insieme_specie[mf[0].str()] = species_index;
							rev_insieme_specie[species_index]=mf[0].str();
							species_index++;
						}
						vettore_regole.push_back( insieme_specie[mf[0].str()] );

						sublinea = mf.suffix().str();
					}

				linea = mr.suffix().str();

				vettore_regole.push_back( FINE_INIBITORI );

				// substep 3: parsa prodotti
				std::regex_search(linea, mr, e);
				//printf(" * Products: ");
				// printf("%s\n", mr[0].str().c_str());			

					// subsubstep 3: snocciola reagenti	
					sublinea = mr[0].str();					
					while( std::regex_search(sublinea, mf, f) ) {
						// printf("%s\n", mf[0].str().c_str());			

						if (insieme_specie.count(mf[0].str())==0) {
							insieme_specie[mf[0].str()] = species_index;
							rev_insieme_specie[species_index]=mf[0].str();
							species_index++;
						}
						vettore_regole.push_back( insieme_specie[mf[0].str()] );

						sublinea = mf.suffix().str();
					}
				
				vettore_regole.push_back( FINE_REAZIONE );
			
			// printf("\n");
			}


		}

		// printf("OK!\n");


	} else {
		perror ("ERROR reading input file");
		return false;
	}

	myfile.close();

	for (unsigned int i=0; i<this->get_number_of_species(); i++) this->vettore_stati.push_back(0);
	for (unsigned int i=0; i<this->get_number_of_species(); i++) this->vettore_stati.push_back(0);

	if (path_m0!="") {

		myfile.open( path_m0 );
		if (myfile.is_open()) {

			std::string linea;

			while(myfile.good()) {
				getline(myfile, linea);
				linea = "\""+linea.substr(0, linea.size())+"\"";
				// printf( "%u\n", this->insieme_specie[linea] );
				this->vettore_stati[ this->insieme_specie[linea] ] = 1;
			}
		} else {
			perror ("ERROR reading input m0 file");
			return false;
		}
		myfile.close();

	}

	
	unsigned int detected_iterations = 0;
	myfile.open( path_context );
	if (myfile.is_open()) {

		std::string linea;
				
		while(myfile.good()) {

			std::vector<char> empty;
			empty.resize(this->insieme_specie.size(), 0);

			getline(myfile, linea);			
			std::vector<std::string> v;
			tokenize(linea, v);
			if (v.size()==0) continue;
			detected_iterations++;
			
			if (v[0].compare(".")!=0) {				
				for (unsigned int i=0; i<v.size(); i++) {
					unsigned int indice = this->insieme_specie[v[i].c_str()];
					empty[indice]=1;
					// printf("<WTF> %s, %d.\n", v[i].c_str(), i);					
				}
			}
			for (unsigned int i=0; i<empty.size(); i++){
				vector_context.push_back(empty[i]);
			}
		}
	} else {
		perror ("ERROR reading input m0 file");
		return false;
	}
	myfile.close();
	// printf(" * Detected %d iterations from the context file.\n", detected_iterations);
	this->iterations = detected_iterations;

	return true;

}
