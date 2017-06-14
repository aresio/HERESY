
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <curand_kernel.h>
#include <limits.h>
#include <cuda_runtime.h>

#include "RSParser.h"
#include "optionparser.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>

__constant__ unsigned int DEV_CONST_REACTIONS[8000];
__constant__ unsigned int DEV_CONST_OFFSET[8000];

template <class T>
inline T minimo( T v1, T v2 ) {
	return
		v1 > v2? v2 : v1;		
}

template <typename T, typename T2>
inline T minimo( T v1, T2 v2 ) {
	return
		v1 > v2? v2 : v1;		
}

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Required1(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }
  
  static option::ArgStatus Required2(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Required3(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Numeric1(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }

    static option::ArgStatus Numeric2(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }

};

void dump_grezzo_reazioni(ReactionSystemsParser* rsp) {

	for (unsigned int r=0; r<rsp->get_number_of_reactions(); r++ ) {
		printf(" * Reazione %u: ", r);
		unsigned int pos = rsp->vettore_offset[r];
		while( rsp->vettore_regole[pos]!=FINE_REAZIONE ) {
			printf("%u, ", rsp->vettore_regole[pos]);
			pos++;
		}
		printf("\n");
	}
	printf("\n");

}

void dump_grezzo_stato_iniziale(ReactionSystemsParser* rsp) {
	printf(" * Dumping initial state of the RS:\n");
	for (unsigned int s=0; s<rsp->get_number_of_species(); s++) {
			printf("%u\t", rsp->vettore_stati[s]);
		}
	printf("\n");
}

void start_profiling(cudaEvent_t* start, cudaEvent_t* stop) {

	/// TIMER 1	
	cudaEventCreate(start);  
	cudaEventCreate(stop);
	cudaEventRecord(*start, 0);

}

float stop_profiling(cudaEvent_t* start, cudaEvent_t* stop) {

	cudaEventRecord( *stop, 0 );
	cudaEventSynchronize( *stop );
	float tempo;
	cudaEventElapsedTime( &tempo, *start, *stop );
	tempo /= 1000;
	// printf("Tempo di esecuzione: %f s\n", tempo);
	return tempo;

}

template <bool use_context>
__global__ void Context( char* Stato, unsigned int num_stato, unsigned int num_sostanze, char* dev_context, unsigned int step ) {

	unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;

	if ( blockIdx.x == gridDim.x-1 ) {
		if (tid>num_sostanze-1) return;
	}

	char stato_corrente = num_stato;
	char stato_prossimo = num_stato ^ 1;

	// reset the next state for results
	Stato[ tid + stato_prossimo*num_sostanze ] = 0;

	// overwrite the CURRENT state with context
	if (use_context) {
		char val = dev_context[tid+step*num_sostanze];
		if (val==1) Stato[ tid + stato_corrente*num_sostanze ] = val;
	} 

//	printf("TID %d current state %d.\n", tid,   Stato[ tid + stato_corrente*num_sostanze ]);
//	printf("TID %d next state %d.\n", tid,		Stato[ tid + stato_prossimo*num_sostanze ]);

}

template<bool use_const>
__global__ void Simulate_Lightweight( 
	const unsigned int* Reazioni, char* Stato, const unsigned int* offset, 
	const unsigned int num_stato, const unsigned int reazioni, const unsigned int num_sostanze) {

	unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;

	if ( blockIdx.x == gridDim.x-1 ) {
		if (tid>reazioni-1) return;
	}

	unsigned int pos;
	if (use_const) 
		pos = DEV_CONST_OFFSET[tid];
	else
		pos = offset[tid];

	char stato_corrente = num_stato ^ 1;
	char stato_prossimo = num_stato ;

	bool reagente = Stato[ stato_corrente*num_sostanze + DEV_CONST_REACTIONS[pos] ];
	bool inibitore = Stato[ stato_corrente*num_sostanze + DEV_CONST_REACTIONS[pos+2] ];			

	if ( reagente && (!inibitore) ) {
		Stato[ stato_prossimo*num_sostanze + DEV_CONST_REACTIONS[pos+4] ] = 1;
	} else {
		Stato[ stato_prossimo*num_sostanze + DEV_CONST_REACTIONS[pos+4] ] = 0;
	}
	
}

template<bool use_const>
__global__ void Simulate( 
	const unsigned int* Reazioni, char* Stato, const unsigned int* offset, 
	const unsigned int num_stato, const unsigned int reazioni, const unsigned int num_sostanze) {

	unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;

	if ( blockIdx.x == gridDim.x-1 ) {
		if (tid>reazioni-1) return;
	}

	unsigned int pos;
	if (use_const) 
		pos = DEV_CONST_OFFSET[tid];
	else
		pos = offset[tid];

	char stato_corrente = num_stato ^ 1;
	char stato_prossimo = num_stato ;

	bool risultato = true;


	// processiamo i reagenti
	if (use_const) {
		while(DEV_CONST_REACTIONS[pos]!=FINE_REAGENTI) {		
			risultato &= Stato[ stato_corrente*num_sostanze + DEV_CONST_REACTIONS[pos] ];			
			pos++;
		}
	} else {
		while(Reazioni[pos]!=FINE_REAGENTI) {		
			risultato &= Stato[ stato_corrente*num_sostanze + Reazioni[pos] ];
			pos++;
		}
	}

//	printf("[Reactants] TID: %d result %d.\n", tid, risultato);
	
	pos++;

	// processiamo gli inibitori
	if (use_const) {
		while(DEV_CONST_REACTIONS[pos]!=FINE_INIBITORI) {
			risultato &= !Stato[ stato_corrente*num_sostanze + DEV_CONST_REACTIONS[pos] ];
			pos++;
		}
	} else {
		while(Reazioni[pos]!=FINE_INIBITORI) {
			risultato &= !Stato[ stato_corrente*num_sostanze + Reazioni[pos] ];
			pos++;
		}
	}

//	printf("[Inhibitors] TID: %d result %d.\n", tid, risultato);

	pos++;

	// calcoliamo i risultati
	if (risultato) {
		if (use_const) {
			while(DEV_CONST_REACTIONS[pos]!=FINE_REAZIONE) {
				Stato[ stato_prossimo*num_sostanze +  DEV_CONST_REACTIONS[pos] ] = 1;
				pos++;				
			}
		} else {
			while(Reazioni[pos]!=FINE_REAZIONE) {
				Stato[ stato_prossimo*num_sostanze +  Reazioni[pos] ] = 1;				
				// printf("[Update (res=1)] TID: %d result %d.\n", tid, Stato[ stato_prossimo*num_sostanze +  Reazioni[pos] ]);
				pos++;
			}
		}
	} else {
		// printf("[Update (res=0)] TID: %d result %d.\n", tid, Stato[ stato_prossimo*num_sostanze +  Reazioni[pos] ]);
	}
}


__global__ void SaveTrace(char* state, char* trace, const long unsigned int step, const long unsigned int species, const unsigned int numstato) {

	unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid>species-1) return;	

	trace[ species*step + tid ] = state[ numstato*species + tid ] ; 
	// printf("[Trace] TID %d state %d.\n", tid, trace[species*step + tid ]);

}


unsigned int numBlocchi( unsigned int NUM_THREADS, unsigned int TPB ) {

	return   NUM_THREADS / TPB + 1;
	
}


void calculateGroupsAndStates(const size_t est, const size_t freem, const unsigned int species, unsigned long int* groups, unsigned long int* statesPerGroup) {
	unsigned long int half = floorl(freem/2);
	*statesPerGroup = floorl( half/((sizeof(char))*species) );
	*groups = floorl( est / (*statesPerGroup) ) + 1;
}

enum  optionIndex { HELP, REQUIRED1, REQUIRED2, REQUIRED3, REQUIRED4, NUMERIC1, NUMERIC2, CONSMEMORY, VERBOSE, LIGHTWEIGHT };
const option::Descriptor usage[] = {
	{ HELP,    0,"", "help",    Arg::None,    "  \t--help  \tPrint usage and exit." },
	{ REQUIRED1, 0,"r","rules",Arg::Required1,"  -r <arg>, \t--rules=<arg>  \tInput file specifying the rules." },
	{ REQUIRED2, 0,"i","initial",Arg::Required2,"  -i <arg>, \t--initial=<arg>  \tInput file specifying the initial state of the system." },
	{ REQUIRED3, 0,"o","output",Arg::Required2,"  -o <arg>, \t--output=<arg>  \tInput file specifying the output file of the trace." },
	{ REQUIRED4, 0,"c","context",Arg::Required2,"  -c <arg>, \t--context=<arg>  \tInput file specifying the context of the reaction system." },
	{ NUMERIC1, 0,"s","steps", Arg::Numeric1, "  -s <num>, \t--steps=<num>  \tRequires a number as argument." },
	{ NUMERIC2, 0,"b","blocks", Arg::Numeric2, "  -b <num>, \t--blocks=<num>  \tRequires a number as argument." },
	{ CONSMEMORY, 0,"n","no_constant", Arg::None, "  -n, \t--no_constant  \tDisables the constant memory." },
	{ VERBOSE, 0,"v","verbose", Arg::None, "  -v, \t--verbose\tEnables verbose mode." },
	{ LIGHTWEIGHT, 0,"l","lightweight", Arg::None, "  -l, \t--lightweight\tEnables lightweight kernel for normal systems." },
	{ 0, 0, 0, 0, 0, 0 } };


int main(int argc, char**argv )
{
	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);

	#ifdef __GNUC__
		// GCC supports C99 VLAs for C++ with proper constructor calls.
		option::Option options[stats.options_max], buffer[stats.buffer_max];
	#else
		// use calloc() to allocate 0-initialized memory. It's not the same
		// as properly constructed elements, but good enough. Obviously in an
		// ordinary C++ program you'd use new[], but this file demonstrates that
		// TLMC++OP can be used without any dependency on the C++ standard library.
		option::Option* options = (option::Option*)calloc(stats.options_max, sizeof(option::Option));
		option::Option* buffer  = (option::Option*)calloc(stats.buffer_max,  sizeof(option::Option));
	#endif

	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())    return 1;
	if (options[HELP] || argc == 0)
	{
		int columns = getenv("COLUMNS")? atoi(getenv("COLUMNS")) : 80;
		option::printUsage(fwrite, stdout, usage, columns);
		return 0;
	}

	std::string input_path("");
	std::string state_path("");
	std::string output_path("");
	std::string context_path("");
	unsigned long int MAX_PASSI = 0;
	unsigned int DIM_BLOCCO = 32;
	bool use_constant = false;
	bool use_context  = false;
	bool verbose = false;
	bool force_disable_constant_memory  = false;
	bool lightweight = false;
	bool output_to_console = true;

	for (int i = 0; i < parse.optionsCount(); ++i)
	{
		option::Option& opt = buffer[i];
		// fprintf(stdout, "Argument #%d is ", i);
		switch (opt.index())
		{
		  case HELP:
			// not possible, because handled further above and exits the program
		  case REQUIRED1:
			input_path = opt.arg;
			break;
		  case REQUIRED2:
			state_path = opt.arg;
			break;
		  case REQUIRED3:			
			output_path = opt.arg;
			break;
		  case REQUIRED4:			
			context_path = opt.arg;
			break;
		  case NUMERIC1:
			MAX_PASSI = atoi(opt.arg);
			break;		  
		  case NUMERIC2:			
			DIM_BLOCCO = atoi(opt.arg);
			break;		  
		  case CONSMEMORY:
			force_disable_constant_memory = true;
			break;
		  case VERBOSE:			
			verbose = true;
			break;
		  case LIGHTWEIGHT:			
			lightweight = true;
			break;
		 /* case UNKNOWN:
			// not possible because Arg::Unknown returns ARG_ILLEGAL
			// which aborts the parse with an error
			break; */
		}
  }

	if (verbose) {
		fprintf(stdout, " * Verbose mode enabled.\n");
		fprintf(stdout, " * RS rules loaded from file: '%s'\n", input_path);
		fprintf(stdout, " * Output file: '%s'\n", output_path);
		fprintf(stdout, " * Context file: '%s'\n", context_path);
		// fprintf(stdout, " * Initial state loaded from file:  '%s'\n", state_path);
		if (force_disable_constant_memory) fprintf(stdout, " * Constant memory disabled.\n");
		if (lightweight) fprintf(stdout, " * Lightweight kernel enabled.\n");
		// fprintf(stdout, " * --blocks with argument '%s'\n", DIM_BLOCCO);
		// fprintf(stdout, " * Simulation steps: %s\n", MAX_PASSI);		
	}




	ReactionSystemsParser rsp;
	// rsp.OpenFile(argv[1], argv[2]);	
	rsp.OpenFile(input_path, state_path, context_path);

	if (rsp.vector_context.size()>0) use_context = true;
	if (use_context && verbose) printf(" * Using context.\n");

	/* arguments: 
		- rules file
		- initial state file
		- running steps
		- block dimension
	*/

	// query device (0) for properties
	cudaDeviceProp devProp;
    cudaGetDeviceProperties(&devProp, 0);

	// force number of steps?
	if (MAX_PASSI==0) {
		MAX_PASSI = rsp.get_iterations();				
	} else {
		MAX_PASSI = minimo(MAX_PASSI, rsp.get_iterations());		
	}
	if (verbose) fprintf(stdout, " * Simulation will perform %d iterations.\n", MAX_PASSI);

	unsigned long int estimatedGM;
	// unsigned long int bytes_per_stream;
	estimatedGM = MAX_PASSI * (unsigned long int)(rsp.get_number_of_species()) * (unsigned long int)(sizeof(char)) ;

	size_t freem, total;
	cudaMemGetInfo(&freem, &total);  

	if (verbose) {
		printf(" * Total global memory: %lu\n", total);
		printf(" * Available global memory: %lu\n", freem);
		printf(" * Estimation of the global memory needed to store the dynamics of the RS: %lu\n", estimatedGM);
		fprintf(stdout, " * Number of detected chemical species: %d.\n", rsp.get_number_of_species() );
		fprintf(stdout, " * Number of detected reactions: %d.\n", rsp.get_number_of_reactions() );
	}

	unsigned long int groups = 0;
	unsigned long int statesPerGroup = 0;

	if ( freem < estimatedGM )  {
		//printf("WARNING: output file is larger than GPU's global memory: using multiple streams.\n");		
		/*
			Dividiamo l'output in due gruppi: stream 0 e stream 1 (lettura durante simulazione).
			I gruppi sono spaccati in più sottoinsiemi sequenziali.
			Dobbiamo determinare quante iterazioni copre ogni gruppo.
		*/		
		calculateGroupsAndStates(estimatedGM, freem, rsp.get_number_of_species(), &groups, &statesPerGroup);		
	} else {
		groups = 1;
		statesPerGroup = MAX_PASSI;
	}

	if (verbose) fprintf(stdout, " * Using %lu groups with %lu states\n", groups, statesPerGroup );


	// TODO
	groups = 1;
	statesPerGroup = MAX_PASSI;
	if (verbose) fprintf(stdout, "WARNING: multiple groups and streams disabled.\n");


	char* dev_results[2];	
	char* host_results[2];	

	if (groups==1) {
		cudaMalloc(&dev_results[0], sizeof(char)*rsp.get_number_of_species() );
		host_results[0] = (char*) malloc ( sizeof(char) * rsp.get_number_of_species() );		
	} else {
		cudaMalloc(&dev_results[0], sizeof(char)*statesPerGroup );
		cudaMalloc(&dev_results[1], sizeof(char)*statesPerGroup );
		host_results[0] = (char*) malloc ( sizeof(char)*statesPerGroup );
		host_results[1] = (char*) malloc ( sizeof(char)*statesPerGroup );
	}
				
	if (verbose) printf(" * Requested %u blocks\n", numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO) );

	char* dev_stato;
	unsigned int* dev_regole;
	unsigned int* dev_offset;
	char* read_back = (char*) malloc ( sizeof(char)*rsp.get_number_of_species()*2 );

	cudaMalloc( &dev_stato, sizeof(char) * rsp.get_number_of_species() * 2 );	
	cudaMalloc( &dev_regole, sizeof(unsigned int) * rsp.vettore_regole.size() );
	cudaMalloc( &dev_offset, sizeof(unsigned int) * rsp.vettore_offset.size() );

	if (verbose) fprintf(stdout, " * Rules vector requires %d bytes.\n", sizeof(unsigned int) * rsp.vettore_regole.size() );
	if (verbose) fprintf(stdout, " * Offsets vector requires %d bytes.\n", sizeof(unsigned int) * rsp.vettore_offset.size() );
	if ( sizeof(unsigned int) * rsp.vettore_regole.size() < 8000) {
		if (verbose) fprintf(stdout, " * It could possible to exploit the constant memory to store the model.\n");
		unsigned int model_size = sizeof(unsigned int) * rsp.vettore_regole.size() ;
		cudaMemcpyToSymbol(DEV_CONST_REACTIONS, &(rsp.vettore_regole[0]), sizeof(unsigned int) * rsp.vettore_regole.size()); 
		cudaMemcpyToSymbol(DEV_CONST_OFFSET, &(rsp.vettore_offset[0]), sizeof(unsigned int) * rsp.vettore_offset.size()); 
		use_constant = true;
		if (force_disable_constant_memory) {
			use_constant=false; // override
			if (verbose) fprintf(stdout, "WARNING: constant memory disabled.\n");
		}
	}
	
	// memory context
	char* dev_context;
	// unsigned int* dev_context_offset;
	cudaMalloc( &dev_context, sizeof(char)*rsp.vector_context.size() );
	// cudaMalloc( &dev_context_offset, sizeof(unsigned int)*rsp.vector_context_offset.size() );
	cudaMemcpy( dev_context, &(rsp.vector_context[0]), sizeof(char)*rsp.vector_context.size() , cudaMemcpyHostToDevice );
	// cudaMemcpy( dev_context_offset, &(rsp.vector_context_offset[0]), sizeof(unsigned int)*rsp.vector_context_offset.size() , cudaMemcpyHostToDevice );
	
	cudaMemcpy( dev_regole, &(rsp.vettore_regole[0]), sizeof( unsigned int )*rsp.vettore_regole.size(), cudaMemcpyHostToDevice );
	cudaMemcpy( dev_offset, &(rsp.vettore_offset[0]), sizeof( unsigned int )*rsp.vettore_offset.size(), cudaMemcpyHostToDevice );	
	cudaMemcpy( dev_stato,  &(rsp.vettore_stati[0]),  sizeof( char ) * rsp.get_number_of_species()*2 , cudaMemcpyHostToDevice );

	// memoria per la traccia
	char* host_trace = (char*) malloc ( sizeof(char) * rsp.get_number_of_species() * MAX_PASSI );
	char* dev_trace;
	memset(host_trace, 0, sizeof(char) * rsp.get_number_of_species() * MAX_PASSI);
	cudaMalloc( &dev_trace, sizeof(char) * rsp.get_number_of_species() * MAX_PASSI );
	cudaMemcpy( dev_trace, host_trace, sizeof(char) * rsp.get_number_of_species() * MAX_PASSI, cudaMemcpyHostToDevice );
	
	// Profiling 
	cudaEvent_t start,  stop;
	start_profiling(&start, &stop);

	//printf(" * Initial state loaded from input files:\n");
	if (verbose)  {
		for (unsigned int i=0; i<rsp.get_number_of_species()*2; i++) {
			printf("Species %d: %d\t", i, rsp.vettore_stati.at(i));
		}
		printf("\n\n");
	}	

	unsigned int num_stato = 0;

	/*
	if (use_context) 
		Context<true><<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>(dev_stato, num_stato^1, rsp.get_number_of_species(), dev_context, 0);
	else
		Context<false><<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>(dev_stato, num_stato^1, rsp.get_number_of_species(), dev_context, 0);

	// write on trace
	SaveTrace<<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>( dev_stato, dev_trace, 0, rsp.get_number_of_species(), num_stato );

	*/

	for (unsigned int i=0; i<MAX_PASSI; i++) {
		/*
		printf(" * Switch state: %d.\n", num_stato);
		

		
				// dump a video dei risultati		
		cudaMemcpy( read_back, dev_stato, sizeof(char)*rsp.get_number_of_species()*2, cudaMemcpyDeviceToHost );
		cudaThreadSynchronize();
		for (unsigned int s=0; s<rsp.get_number_of_species()*2; s++) {
			printf("[Readback] Species %d state %d.\n", s, read_back[s+(num_stato^1)*rsp.get_number_of_species()]);
		}
		printf("\n");

		*/
		if (use_context) 
			Context<true><<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>(dev_stato, num_stato^1, rsp.get_number_of_species(), dev_context, i);
		else
			Context<false><<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>(dev_stato, num_stato^1, rsp.get_number_of_species(), dev_context, i);

		// Simulate<<< numBlocchi(rsp.get_number_of_reactions(), DIM_BLOCCO) , DIM_BLOCCO, 0, streams[selStream] >>>( dev_regole, dev_stato, dev_offset, num_stato, rsp.get_number_of_reactions(), rsp.get_number_of_species() );
		if (lightweight) {
			Simulate_Lightweight<true><<< numBlocchi(rsp.get_number_of_reactions(), DIM_BLOCCO) , DIM_BLOCCO >>>( dev_regole, dev_stato, dev_offset, num_stato, rsp.get_number_of_reactions(), rsp.get_number_of_species() );
		} else {
			if (use_constant)
				Simulate<true><<< numBlocchi(rsp.get_number_of_reactions(), DIM_BLOCCO) , DIM_BLOCCO >>>( dev_regole, dev_stato, dev_offset, num_stato, rsp.get_number_of_reactions(), rsp.get_number_of_species() );
			else
				Simulate<false><<< numBlocchi(rsp.get_number_of_reactions(), DIM_BLOCCO) , DIM_BLOCCO >>>( dev_regole, dev_stato, dev_offset, num_stato, rsp.get_number_of_reactions(), rsp.get_number_of_species() );
		}
		// write on trace
		SaveTrace<<< numBlocchi(rsp.get_number_of_species(), DIM_BLOCCO), DIM_BLOCCO>>>( dev_stato, dev_trace, i, rsp.get_number_of_species(), num_stato );

		cudaThreadSynchronize();

		// dump a video dei risultati		
		/*
		cudaMemcpy( read_back, dev_stato, sizeof(char)*rsp.get_number_of_species()*2, cudaMemcpyDeviceToHost );
		cudaThreadSynchronize();
		for (unsigned int s=0; s<rsp.get_number_of_species()*2; s++) {
			printf("[Readback] Species %d state %d.\n", s, read_back[s+(num_stato^1)*rsp.get_number_of_species()]);
		}
		printf("\n");
		*/
		
		
		num_stato ^= 1;			

		/*
		if (++selCount == statesPerGroup-1 ) {
			selCount = 0;

			unsigned long int copy_bytes = sizeof(char) * minimo( MAX_PASSI-(blockCount*statesPerGroup), statesPerGroup );			

			// memcpy asincrona
			cudaMemcpyAsync( host_results[selCount], dev_results[selCount], copy_bytes, cudaMemcpyDeviceToHost, streams[selStream]);
			
			selStream ^= 1;
			blockCount ++;
		}
		*/
		
	}

	stop_profiling(&start, &stop);

	cudaMemcpy(host_trace, dev_trace, sizeof(char) * rsp.get_number_of_species() * MAX_PASSI, cudaMemcpyDeviceToHost);

	if (!output_to_console) {
		std::ofstream output_file(output_path.c_str());
		if (output_file.is_open()) {

			output_file << "0" << "\t" ;

			for (unsigned int s=0; s<rsp.get_number_of_species(); s++) {
				if (rsp.vector_context[ s ]) {
					output_file << rsp.rev_insieme_specie[s]  << "\t";
				}
			}
			output_file << std::endl;

			for (unsigned int step=0; step<MAX_PASSI; step++) {
				output_file << step+1 << "\t";
				for (unsigned int species =0; species < rsp.get_number_of_species(); species++) {
					if ( host_trace[ rsp.get_number_of_species()*step + species ] == 1 ) {
						output_file << rsp.rev_insieme_specie[species] << "\t";
					}				
				}
				output_file << "\n";
			}
			output_file.close();
		} else {
			perror("ERROR: cannot save output dynamics to file.\n");
			exit(-1);
		}
	}  else { // output to console

		std::cout << "0 ";
		for (unsigned int s=0; s<rsp.get_number_of_species(); s++) {			
			if (rsp.vector_context[ s ]) {
				std::cout << rsp.rev_insieme_specie[s]  << " ";
			}
		}
		std::cout << std::endl;

		for (unsigned int step=0; step<MAX_PASSI; step++) {
			std::cout << step+1 << " ";
			for (unsigned int species =0; species < rsp.get_number_of_species(); species++) {
				if ( host_trace[ rsp.get_number_of_species()*step + species ] == 1 ) {
					std::cout << rsp.rev_insieme_specie[species] << " ";
				}				
			}
			std::cout << "\n";
		}

	}

	// system("pause");
    return 0;
}
