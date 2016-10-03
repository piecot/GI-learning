/*
 ============================================================================
 Name        : GI-learning.cpp
 Author      : Pietro Cottone, Marco Ortolani and Gabriele Pergola
 Version     :
 Copyright   : You can do whatever you want with this code, if you cite me ;)
 Description : 
 ============================================================================
 */

#include <ctime>
#include <iostream>
#include <string>
#include "dfa.h"
#include "bluefringe.h"
#include "rpni.h"
#include "edsm.h"
#include "lstar.h"
#include "blueStar.h"
#include "idsm.h"
#include "messages.h"
#include <boost/filesystem.hpp>


#define SAMPLE_DIR "examples" 								// training set: put your training data here

#define EDSM_FILE "examples.txt" 							// training set: put your training data here
#define EDSM_FILE_BIG "examples_big.txt" 				// training set: put your training data here

#define LSTAR_FILE "lstar.txt"										// file for lstar
#define LSTAR_FILE_BIG "lstar_big.txt" 						// file for lstar

#define DOT_DIR "results"											// dot file with inferred DFA
#define DOT_FILE_IDSM "inferredDFA_IDSM.dot"
#define DOT_FILE_BLUESTAR "inferredDFA_bluestar.dot"
#define DOT_FILE_EDSM "inferredDFA_edsm.dot"
#define DOT_FILE_RPNI "inferredDFA_rpni.dot"
#define DOT_FILE_LSTAR "inferredDFA_lstar.dot"
//#define EDSM_RES_PATH "DFA_dot" 							// by-products of inference

#define MAX_BUFFER_SIZE 256
#define BASE_PATH_EXE 0											// the working dir is where the executable is

#define MAX_ARGC 3
#define MIN_ARGC 1

using namespace std;

namespace fs=boost::filesystem;

void parse_input_args(int, char**, string *, string *);
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
	gi::rpni* 		rpni_exe = NULL;
	gi::edsm* 		edsm_exe = NULL;
	gi::blueStar*	bluestar_exe = NULL;
	gi::idsm* 		idsm_exe = NULL;


	clock_t tStart;

	//file
	string example_file="";
	string lstar_dfa_file="";

	// working dir
	string base_path;
	string res_path;

	//parse input arguments
	parse_input_args(argc, argv, &example_file, &lstar_dfa_file);

	// Setting the locale, in order to choose the correct language
	//string curr_os_locale = setlocale(LC_CTYPE, "");
	//cout<<"Current locale: "<<curr_os_locale<<endl;

	if(BASE_PATH_EXE){
		// complete path of the executable
		fs::path selfpath=fs::system_complete(argv[0]);

		// path of the folder where the executable is: it is assumed as the working directory
		base_path = selfpath.parent_path().c_str();
	}else{
		// the working dir is the current dir
		base_path = fs::current_path().c_str();
	}



	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	res_path = check_res_dir(base_path);
	res_path = res_path + fs::path::preferred_separator;


	cout<<"Working dir: "<< base_path <<endl;


int a=0;
if(a != 0)
{

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- IDSM algorithm ---
	cout << endl<< "************************";
	cout << endl<< "********  IDSM *********" << endl;
	cout << "************************" << endl;

	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;


	// print current example file
	cout<<"Example file: "<< example_file <<endl;

	idsm_exe = new gi::idsm(example_file.c_str(), 0.05, 1000.0);


	gi::dfa* IDSM_dfa = idsm_exe->run(res_path);


	// Print transition table of the inferred automaton
	IDSM_dfa->print_dfa_ttable("- IDSM dfa -");


	// Create dot figure for the inferred automaton
	IDSM_dfa->print_dfa_dot_mapped_alphabet("IDSM", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_IDSM).c_str());

	//Error rate
	cout << "Error-rate: "<< idsm_exe->get_error_rate_final_dfa() << endl;

	// free allocated memory
	if(idsm_exe!=NULL)
		delete idsm_exe;




	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- Blue* algorithm ---
	cout << endl<< "************************";
	cout << endl<< "********  BLUESTAR *********" << endl;
	cout << "************************" << endl;

	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;


	// print current example file
	cout<<"Example file: "<< example_file <<endl;

	bluestar_exe = new gi::blueStar(example_file.c_str(), 0.05, 1000.0);


	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	//res_path = check_res_dir(base_path);
	//res_path = res_path + fs::path::preferred_separator;


	gi::dfa* BLUESTAR_dfa = bluestar_exe->run(res_path);


	// Print transition table of the inferred automaton
	BLUESTAR_dfa->print_dfa_ttable("- BlueStar dfa -");

	// Create dot figure for the inferred automaton
	BLUESTAR_dfa->print_dfa_dot_mapped_alphabet("BlueStar", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_BLUESTAR).c_str());

	//Error rate
	cout << "Error-rate: "<< bluestar_exe->get_error_rate_final_dfa() << endl;

	// free allocated memory
	if(bluestar_exe!=NULL)
		delete bluestar_exe;

}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- LSTAR algorithm ---
	cout << endl<< "*************************";
	cout << endl<< "********  LSTAR  ********" << endl;
	cout <<  "*************************" << endl;

	gi::lstar* lstar_exe = new gi::lstar(gi::dfa::read_dfa_file(lstar_dfa_file));

	// start timer to compute execution time
	tStart = clock();

	gi::dfa* L_dfa = lstar_exe->run(false, "");

	// Stat lstar
	//n_memb_query_lstar = L_dfa->get_n_memb_query();

	cout<<"Time taken to LSTAR: "<< (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;

	L_dfa->print_dfa_dot_mapped_alphabet("LSTAR", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_LSTAR).c_str());


	delete lstar_exe;


	exit(EXIT_SUCCESS);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- RPNI algorithm ---
	cout << endl<< "************************";
	cout << endl<< "********  RPNI *********" << endl;
	cout << "************************" << endl;

	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;


	// print current example file
	cout<<"Example file: "<< example_file <<endl;

	rpni_exe = new gi::rpni(example_file.c_str());


	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	//res_path = check_res_dir(base_path);

	res_path = res_path + fs::path::preferred_separator;

	// start timer to compute execution time
	tStart = clock();

	gi::dfa* RPNI_dfa = NULL;
	double millisec=  rpni_exe->run_elapsed_time(res_path,  &RPNI_dfa );
	//gi::dfa* RPNI_dfa = rpni_exe->run(res_path);

	// Print execution time
	cout<<"Time taken to RPNI: " <<(double)(clock() - tStart)/CLOCKS_PER_SEC << endl;


	// Print transition table of the inferred automaton
	RPNI_dfa->print_dfa_ttable("- Rpni dfa -");

	// Create dot figure for the inferred automaton
	RPNI_dfa->print_dfa_dot_mapped_alphabet("RPNI", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_RPNI).c_str());

	cout << "DFA di RPNI stampato in: "<<base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_RPNI <<endl;

	// free allocated memory
	if(rpni_exe!=NULL)
		delete rpni_exe;



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- EDSM algorithm ---
	cout << endl<< "************************";
	cout << endl<< "********  EDSM *********" << endl;
	cout << "************************" << endl;

	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;


	// print current example file
	cout<<"Example file: "<< example_file <<endl;

	edsm_exe = new gi::edsm(example_file.c_str());


	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	//res_path = check_res_dir(base_path);

	res_path = res_path + fs::path::preferred_separator;

	// start timer to compute execution time
	tStart = clock();

	gi::dfa* EDSM_dfa = edsm_exe->run(res_path);

	// Print execution time
	cout<<"Time taken to EDSM: " <<(double)(clock() - tStart)/CLOCKS_PER_SEC << endl;


	// Print transition table of the inferred automaton
	EDSM_dfa->print_dfa_ttable("- Edsm dfa -");

	// Create dot figure for the inferred automaton
	EDSM_dfa->print_dfa_dot_mapped_alphabet("EDSM", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_EDSM).c_str());

	cout << "DFA di EDSM stampato in: "<<base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_EDSM <<endl;

	// free allocated memory
	if(edsm_exe!=NULL)
		delete edsm_exe;


	if(RPNI_dfa != NULL)
		delete RPNI_dfa;
	if(EDSM_dfa != NULL)
		delete EDSM_dfa;
//	if(BLUESTAR_dfa != NULL)
//		delete BLUESTAR_dfa;
	if(L_dfa != NULL)
		delete L_dfa;
}








void parse_input_args(int argc, char* argv[], string *bs, string *ls){
	if(argc>MAX_ARGC || argc<MIN_ARGC){
		cerr<<MSG_WRONG_ARGC<<endl;

		exit(EXIT_FAILURE);
	}

	if(argc>=2){
		if(!strcmp("little", argv[1])){
			(*bs) = (*bs) + EDSM_FILE;
		}else if (!strcmp("big", argv[1])){
			(*bs) = (*bs) + EDSM_FILE_BIG;
		}else{
			//cerr<<MSG_WRONG_ARGV<< argv[1] <<endl;
			(*bs) = argv[1];
		}

		if(3==argc){
			if(!strcmp("little", argv[1])){
				(*ls) = (*ls) + LSTAR_FILE;
			}else if (!strcmp("big", argv[1])){
				(*ls) = (*ls) + LSTAR_FILE_BIG;
			}else{
				(*ls) = (*ls) + argv[2];
			}
		}
	}else{
		(*bs) = (*bs) + EDSM_FILE;
		(*ls) = (*ls) + LSTAR_FILE;
	}


}


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

		fs::rename(res_path, res_path + "_" + append_time);
	}

	// create res_dir
	fs::create_directory( res_path );

	res_path = res_path + fs::path::preferred_separator;

	return res_path;
}





//	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/// W characterized set
//
//	gi::dfa dfa128 = gi::dfa::read_dfa_file("../examples/128.txt");
//	gi::dfa dfa144 = gi::dfa::read_dfa_file("../examples/144.txt");
//	gi::dfa dfa153 = gi::dfa::read_dfa_file("../examples/153.txt");
//
//	cout << "Size 128: "<<dfa128.get_num_states() << endl;
//	cout << "Size 144: "<<dfa144.get_num_states() << endl;
//	cout << "Size 153: "<<dfa153.get_num_states() << endl;
//
//
//
////	SYMBOL* table = targetdfa.table_filling();
////
////	vector<int>* eq_states = targetdfa.equivalent_states_list_from_table(table);
////
////	cout << "STATI EQUIVALENTI: "<<endl;
////	for(auto it=eq_states->begin(); it != eq_states->end(); ++it)
////		cout << "Stato "<< *it << endl;
////
////
////
////	vector<vector<SYMBOL>> cover_set = targetdfa.get_cover_set();
////
////	cout << "Cover set: "<<endl;
////	for(auto &it: cover_set)
////		cout << it << endl;
////
//
//
////	cout << "Characterization set: "<<endl;
////	vector<vector<SYMBOL>> char_set = targetdfa.get_characterization_set();
////
////	for(auto it=char_set.begin(); it != char_set.end(); ++it)
////		cout << (*it) << endl;
////
////
////
////	int num_states_target_dfa = 2;
////
////	vector<vector<SYMBOL>> aug_characterization_set = targetdfa.get_augmented_characterization_set(num_states_target_dfa);
////
////	cout << "AUGMENTED augmented_characterization_set: "<<endl;
////	for(auto &it: aug_characterization_set)
////		cout << it << endl;
////
////
////
//
//
//	cout << "W-method TEST SET: "<<endl;
//	vector<vector<SYMBOL>> test_set_128 = dfa128.get_w_method_test_set(dfa153.get_num_states());
//	vector<vector<SYMBOL>> test_set_144 = dfa144.get_w_method_test_set(dfa153.get_num_states());
//
//	vector<vector<SYMBOL>> test_set_153_128 = dfa153.get_w_method_test_set(dfa128.get_num_states());
//	vector<vector<SYMBOL>> test_set_153_144 = dfa153.get_w_method_test_set(dfa144.get_num_states());
//
//	//for(auto &it1 : test_set)
//	//	cout << it1 << endl;
//
//
//	cout << "Membership query: "<<endl;
//	int true_positive[2] = {0, 0};
//	int true_negative[2] = {0, 0};
//	int false_positive[2] = {0, 0};
//	int false_negative[2] = {0, 0};
//
//	double precision[2] = {0, 0};
//	double recall[2] = {0, 0};
//	double fmeasure[2] = {0, 0};
//
//	int index_user =0;
//	for(auto &it1 : test_set_144)
//	{
//		if(dfa144.membership_query(it1))
//		{
//			if(!dfa153.membership_query(it1))
//				false_positive[index_user]++;
//			else
//				true_positive[index_user]++;
//		}
//		else
//		{
//			if(dfa153.membership_query(it1))
//				false_negative[index_user]++;
//			else
//				true_negative[index_user]++;
//		}
//	}
//
//
//	index_user =1;
//	for(auto &it1 : test_set_128)
//	{
//		if(dfa128.membership_query(it1))
//		{
//			if( !dfa153.membership_query(it1))
//				false_positive[index_user]++;
//			else
//				true_positive[index_user]++;
//		}
//		else
//		{
//			if(dfa153.membership_query(it1))
//				false_negative[index_user]++;
//			else
//				true_negative[index_user]++;
//		}
//	}
//
//
//	for(int i=0; i<2; ++i)
//	{
//		if(i==1)
//			cout << "128" << endl;
//		else
//			cout << "144" << endl;
//		cout << "TP: " << true_positive[i] << endl;
//		cout << "TN: " << true_negative[i] << endl;
//		cout << "FP: " << false_positive[i] << endl;
//		cout << "FN: " << false_negative[i] << endl;
//		cout << endl;
//
//		precision[i] = (double) true_positive[i] / (double)(true_positive[i]+false_positive[i]);
//		recall[i] = (double) true_positive[i] / (double)(true_positive[i]+false_negative[i]);
//		fmeasure[i] = (double) (2.0 * precision[i] * recall[i]) / (double) (precision[i] + recall[i]) ;
//
//		cout << "Precision: "<<precision[i]<<endl;
//		cout << "Recall: "<<recall[i]<<endl;
//		cout << "Fmeasure: "<<fmeasure[i]<<endl;
//		cout << endl;
//		cout << endl;
//	}
////
////	cout << "Membership query: "<<endl;
////	int true_positive[2] = {0, 0};
////	int true_negative[2] = {0, 0};
////	int false_positive[2] = {0, 0};
////	int false_negative[2] = {0, 0};
////
////	double precision[2] = {0, 0};
////	double recall[2] = {0, 0};
////	double fmeasure[2] = {0, 0};
////
////	int index_user =0;
////	for(auto &it1 : test_set_153_144)
////	{
////		if(dfa153.membership_query(it1))
////		{
////			if(dfa144.membership_query(it1))
////				true_positive[index_user]++;
////			else
////				false_negative[index_user]++;
////		}
////		else
////		{
////			if(dfa144.membership_query(it1))
////				false_positive[index_user]++;
////			else
////				true_negative[index_user]++;
////		}
////	}
////
////
////	index_user =1;
////	for(auto &it1 : test_set_153_128)
////	{
////		if(dfa153.membership_query(it1))
////		{
////			if(dfa128.membership_query(it1))
////				true_positive[index_user]++;
////			else
////				false_negative[index_user]++;
////		}
////		else
////		{
////			if(dfa128.membership_query(it1))
////				false_positive[index_user]++;
////			else
////				true_negative[index_user]++;
////		}
////	}
////
////
////	for(int i=0; i<2; ++i)
////	{
////		if(i==1)
////			cout << "128" << endl;
////		else
////			cout << "144" << endl;
////
////		cout << "TP: " << true_positive[i] << endl;
////		cout << "TN: " << true_negative[i] << endl;
////		cout << "FP: " << false_positive[i] << endl;
////		cout << "FN: " << false_negative[i] << endl;
////		cout << endl;
////
////		precision[i] = (double) true_positive[i] / (double)(true_positive[i]+false_positive[i]);
////		recall[i] = (double) true_positive[i] / (double)(true_positive[i]+false_negative[i]);
////		fmeasure[i] = (double) (2.0 * precision[i] * recall[i]) / (double) (precision[i] + recall[i]) ;
////
////		cout << "Precision: "<<precision[i]<<endl;
////		cout << "Recall: "<<recall[i]<<endl;
////		cout << "Fmeasure: "<<fmeasure[i]<<endl;
////		cout << endl;
////		cout << endl;
////	}
//
//
//
//
////	map<int, vector<SYMBOL>> access_strings = targetdfa.get_access_strings();
////
////	cout << "Access strings: "<<endl;
////	for(auto it=access_strings.begin(); it != access_strings.end(); ++it)
////		cout << it->second << endl;
////
////	vector<string> mapped_access_strings =  targetdfa.get_access_strings_with_alphabet_symbols();

