/*
 ============================================================================
 Name        : GI-learning.cpp
 Author      : Pietro Cottone
 Version     :
 Copyright   : You can do whatever you want with this code, if you cite me ;)
 Description : Hello World in C++,
 ============================================================================
 */
 
//MIO Usare -std=c++11 per il compilatore versione 11. Altrimenti cose come auto non è possibile usarle
//MIO Inoltre qui questo aspetto è già configurato nel makefile(quindi non bisogna fare niente
#include <ctime>
#include <chrono>
#include <iostream>
#include <string>
#include "dfa.h"
#include "edsm.h"
#include "messages.h"
#include <boost/filesystem.hpp>

#include "lstar.h"
#include "observationpack.h"
#include "esperimenti.h"                                    //QUESTO POI ANDRA' ELIMINATO

#define SAMPLE_DIR "examples" 								// training set: put your training data here

#define EDSM_FILE "examples.txt" 							// training set: put your training data here
#define EDSM_FILE_BIG "examples_big.txt" 				// training set: put your training data here

#define LSTAR_FILE "lstar.txt"										// file for lstar
#define LSTAR_FILE_BIG "lstar_big.txt" 						// file for lstar

#define OBSERVATIONPACK_FILE "observationpack.txt"               //file for observation pack
#define OBSERVATIONPACK_FILE_BIG "observationpack_big.txt"       //file for observation pack

#define DOT_DIR "results"											// dot file with inferred DFA
#define DOT_FILE_EDSM "inferredDFA_edsm.dot"
#define DOT_FILE_LSTAR "inferredDFA_lstar.dot"
#define DOT_FILE_OBSERVATIONPACK "inferredDFA_observationpack.dot"
//#define EDSM_RES_PATH "DFA_dot" 							// by-products of inference

#define COUNTEREXAMPLE_FUNCTION "rivest-shapire"           //default counterexample analysis method for observation pack

#define MAX_BUFFER_SIZE 256
#define BASE_PATH_EXE 0											// the working dir is where the executable is

#define MAX_ARGC 5
#define MIN_ARGC 4                                             //In addition to the executable you have to pass at least 3 strings

using namespace std;

namespace fs=boost::filesystem;

void parse_input_args(int, char**, string *, string *, string *, string *);
void parse_input_args_aux(string * , const char* , const char* , const char*);
string check_res_dir(const string);

/**
 *
 * @param argc
 * @param argv
 * @return
 */


int main(int argc, char* argv[]){
	//dfa
	//dfa EDSM_dfa = NULL;
	gi::edsm* edsm_exe = NULL;

	clock_t tStart;

	//file
	string edsm_example_file="";
	string lstar_dfa_file="";
    string observationpack_dfa_file="";
    string obpack_counterexample_function="";
    
	// working dir
	string base_path;
	string res_path;

	//parse input arguments
	parse_input_args(argc, argv, &edsm_example_file, &lstar_dfa_file , &observationpack_dfa_file, &obpack_counterexample_function);
	// Setting the locale, in order to choose the correct language
	//string curr_os_locale = setlocale(LC_CTYPE, "");
	//cout<<"Current locale: "<<curr_os_locale<<endl;

	if(BASE_PATH_EXE){ //MIO E' SETTATO A ZERO QUINDI VA NELL'ELSE
		// complete path of the executable
		fs::path selfpath=fs::system_complete(argv[0]); //MIO selfpath E' UNA VARIABILE DI TIPO fs:path
		// path of the folder where the executable is: it is assumed as the working directory
		base_path = selfpath.parent_path().c_str();
	}else{
		// the working dir is the current dir
		base_path = fs::current_path().c_str();

	}
	
	cout<<"Working dir: "<< base_path <<endl;
    
    res_path = check_res_dir(base_path);
	//ORIGINALE res_path = res_path + fs::path::preferred_separator;//DOMANDA MIO QUESTO NON E' SUPERFLUO PROVA AD ELIMINARLO (gia' in check_res_dir aggiungo alla fine il preferred_separator)
	
	gi::esperimenti exper; //QUESTO POI ANDRA' ELIMINATO
	
    //ESEGUI EDSM SOLO SE edsm_example_file DIVERSO DA NOTHING. CIOE' SI DA LA POSSIBILITA' ALL'UTENTE INSERENDO NOTHING DI SALTARE L'ESECUZIONE DI EDSM
    if( edsm_example_file.compare("nothing") ) //It's possible compare a sts::string with a litteral string (const *char)
    {
	    // --- EDSM algorithm ---
		cout << endl<< "************************";
		cout << endl<< "********  EDSM *********" << endl;
		cout << "************************" << endl;

		//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;
	
	
		// print current example file
		cout<<"Example file: "<< edsm_example_file <<endl; //MIO edsm_example_file E' IL FILE PASSATO DA RIGA DI COMANDO
	
		edsm_exe = new gi::edsm(edsm_example_file.c_str());
	
		// START TIME
	    std::chrono::time_point<std::chrono::system_clock> start, end;
	    start = std::chrono::system_clock::now();
	
		gi::dfa* EDSM_dfa = edsm_exe->run(res_path);
	    // STOP TIME
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
		// Print execution time
		cout<<"Time taken to EDSM: " <<elapsed_seconds.count()*1.0 << endl;
	
	
		// Print transition table of the inferred automaton
		EDSM_dfa->print_dfa_ttable("- Edsm dfa -");
	
		// Create dot figure for the inferred automaton
		// DOMANDA MIO MA IN QUESTO MODO IL FILE DOT NON VIENE MESSO nella cartella results ANZICHE' NELLA CARTELLA INDICATA DA res_path
		// DOMANDA MIO CHE CONTIENE APPESO A results ANCHE DATA E ORA
		//ORIGINALE EDSM_dfa.print_dfa_dot_mapped_alphabet("EDSM", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_EDSM).c_str());
	    EDSM_dfa->print_dfa_dot_mapped_alphabet("EDSM", (res_path + DOT_FILE_EDSM).c_str());
		//MIO LA RIGA SOPRA L'HO AGGIUNTA IO
		// free allocated memory
		if(edsm_exe!=NULL)
			delete edsm_exe;
			
		if(EDSM_dfa!=NULL)
			delete EDSM_dfa;	
	}
    
    gi::dfa* L_dfa;
    if( lstar_dfa_file.compare("nothing") )
    {
		// --- LSTAR algorithm ---
		cout << endl<< "*************************";
		cout << endl<< "********  LSTAR  ********" << endl;
		cout <<  "*************************" << endl;
	
	    
		gi::lstar* lstar_exe = new gi::lstar(gi::dfa::read_dfa_file(lstar_dfa_file));//MIO Legge un dfa dal file lstar_dfa_file. Il processo va fatto dinamicamente perchè il nome del file da cui creo il DFA è inserito a runtime (e non solo)
	
		// START TIME
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	
		L_dfa = lstar_exe->run(false, "");//MIO Qui viene chiamato il costruttore di copia (quando si fa un'assegnazione durante l'inizializzazione viene chiamato. 
	// STOP TIME
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
		exper.writeTime("../esperimenti/experiments.txt",elapsed_seconds.count()*1.0 );
	
		cout<<"Time taken to LSTAR: "<< elapsed_seconds.count()*1.0 << endl;
	    // DOMANDA MIO MA IN QUESTO MODO IL FILE DOT NON VIENE MESSO nella cartella results ANZICHE' NELLA CARTELLA INDICATA DA res_path
		// DOMANDA MIO CHE CONTIENE APPESO A results ANCHE DATA E ORA
		//ORIGINALE L_dfa.print_dfa_dot_mapped_alphabet("LSTAR", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_LSTAR).c_str());
	    //L_dfa->print_dfa_dot_mapped_alphabet("LSTAR", (res_path + DOT_FILE_LSTAR).c_str());
	    if(lstar_exe!=NULL)
		    delete lstar_exe;
		
	}
	
	gi::dfa* OP_dfa;
	if( observationpack_dfa_file.compare("nothing") )
    {
	    // --- OBSERVATION PACK algorithm ---
	    cout << endl<< "************************************";
		cout << endl<< "********  OBSERVATION PACK  ********" << endl;
		cout <<        "************************************" << endl;
		
		gi::observationpack* observationpack_exe = new gi::observationpack(gi::dfa::read_dfa_file(observationpack_dfa_file) , obpack_counterexample_function);
		
		// START TIME
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
		
		OP_dfa = observationpack_exe->run_observationpack(false, "");//MIO Qui viene fatta una copia membro a membro. Quindi per i membri dinamici vengono solo copiati i puntatori ma in questo caso va bene.NO NO viene chiamato il costrutto di copia (quando si inizializza una variabile nella dichiarazione vine chiamato il costruttore di copia)
        // STOP TIME
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
        exper.writeTime("../esperimenti/experiments.txt",elapsed_seconds.count()*1.0 );
		cout<<"Time taken to OBSERVATION PACK: "<< elapsed_seconds.count()*1.0  << endl;
		 //OP_dfa->print_dfa_dot_mapped_alphabet("OBSERVATIONPACK", (res_path + DOT_FILE_OBSERVATIONPACK ).c_str());
		 if(observationpack_exe!=NULL)
		    delete observationpack_exe;
		 
		 
		   	
	}        	
	
		
		if(!obpack_counterexample_function.compare("rivest-shapire"))
		{
	//QUESTO POI ANDRA' ELIMINATO    
		vector<SYMBOL> esito = OP_dfa->equivalence_query(L_dfa);
		if(esito.size() ==0 && (OP_dfa->get_num_states()==L_dfa->get_num_states()))
		cout<<endl<<endl<<"SONO EQUIVALENTI"<<endl;
		else
		cout<<endl<<endl<<"NON SONO EQUIVALENTI"<<endl;
		//QUESTO POI ANDRA' ELIMINATO  	
	}
	
	if(L_dfa!=NULL)
		delete L_dfa;
    if(OP_dfa != NULL)
		delete OP_dfa;
}







//MIO SALVA IN BS ED IN LS ED IN OP I NOMI DEI FILE (EVENTUALMENTE I PERCORSI) USATI PER L'INFERENZA
//SALVA IN CDF LA COUNTEREXAMPLE DECOMPOSITIONFUNCTION
void parse_input_args(int argc, char* argv[], string *bs, string *ls, string *op, string *cdf){
	if(argc>MAX_ARGC || argc<MIN_ARGC){
		cerr<<MSG_WRONG_ARGC<<endl;

		exit(EXIT_FAILURE);
	}

    //MIO CIOE' POSSO PASSARE COME PRIMO PARAMETRO (DOPO IL NOME DEL FILE BINARIO)ANZICHE' IL NOME DEL FILE LITTLE OPPURE BIG PER DIRE CHE IL TRAINING SET SARA' NEL FILE
    //MIO example oppure examples_big.txt   OPPURE POSSO PASSARE DIRETTAMENTE IL NOME DEL FILE COMPLETO (COMPRENSIVO DEL PERCORSO)
	 //MIO LO STESSO DI SOPRA. SE PASSO COME SECONDO PARAMETRO LITTLE L'ALGORITMO LSTAR ASSUMERA' COME DFA TARGET QUELLO
        //CONTENUTO NEL FILE lstar.txt ALTRIMENTI SE BIG QUELLO NEL FILE  lstar_big.txt ALTRIMENTI (NE BIG NE LITTLE) USA
        //DIRETTAMENTE IL NOME DEL FILE COMPLETO (COMPRENSIVO DEL PERCORSO)    
    parse_input_args_aux(bs , argv[1] , EDSM_FILE , EDSM_FILE_BIG);
    parse_input_args_aux(ls , argv[2],  LSTAR_FILE , LSTAR_FILE_BIG);
    parse_input_args_aux(op , argv[3] , OBSERVATIONPACK_FILE , OBSERVATIONPACK_FILE_BIG);
		
	
	if( ! op->compare("nothing") ) //if op è uguale a nothing
	{
	    if(argc == 5)
	        cerr<<MSG_WARNING_ARGV<<": "<<argv[4]<<endl;
	}
	else //op != "nothing"
	    if(argc == 4)
	        (*cdf) = (*cdf) +  COUNTEREXAMPLE_FUNCTION; //Use rivest-shapire
	    else //argc==5
	        (*cdf) = (*cdf) +  argv[4]; 
	
}

void parse_input_args_aux(string *file , const char* file_entered , const char* file_little, const char* file_big)
{
	if(!strcmp("little", file_entered)){ //MIO Nel caso di argv[1], significa che il secondo parametro quello appena dopo il nome dell'eseguibile, è uguale little
			(*file) = (*file) + file_little;
		}else if (!strcmp("big", file_entered)){
			(*file) = (*file) + file_big;
		}else{
			//cerr<<MSG_WRONG_ARGV<< argv[1] <<endl;
			(*file) = (*file) + file_entered;
		}
} 
//MIO SE NELLA DIRECTORY DI LAVORO (base_path)  ESISTE GIA' UNA CARTELLA results I RISULTATI DELL'ESECUZIONE CORRENTE
//MIO VENGONO INSERITI NELLA CARTELLA results_dataeoracorrente
//INOLTRE QUI VIENE CREATA LA DIRECTORY IN QUESTIONE
string check_res_dir(const string  base_path){
	string res_path;
	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	res_path = base_path + fs::path::preferred_separator + DOT_DIR;

	// move res_dir if it exists
	if(fs::exists(res_path)){
		char append_time[MAX_BUFFER_SIZE];
		time_t rawtime = std::time(NULL);

		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );


		strftime(append_time, MAX_BUFFER_SIZE, "%Y_%m_%d_%H_%M_%S", timeinfo);

        res_path = res_path + "_" + append_time;
		//ORIGINALE fs::rename(res_path, res_path + "_" + append_time);
		}

	// create res_dir
	fs::create_directory( res_path );
	
	res_path = res_path + fs::path::preferred_separator;

	return res_path;
}
