/*
 * DFA definition.
 */

#include "dfa.h"
#include <utilities.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>	 //isstringstream
#include <map>
#include <list>
#include <cmath>
#include <algorithm>
#include <limits>


using namespace std;


#define BUFFER_SIZE 1024

#define INF numeric_limits<int>::max()						/*!< Usually adopted for a non initialized state. */

// For Table Filling Algorithm
#define DFA_TF_STATE_N numeric_limits<SYMBOL>::max()		/*!< Unknow status */
#define DFA_TF_STATE_X numeric_limits<SYMBOL>::max() -1  	/*!< Distinct state (accepting vs rejecting) */

// For W-Method
#define LIMIT_OF_TESTSET_W_METHOD 400000000


gi::dfa::dfa(){
	dim_alphabet = 0;
	num_states 	 = 0;
	start_state  = 0;

	ttable = NULL;
	alphabet=NULL;
}


gi::dfa::dfa(const int n_state, const int dim_alf, const string *alf, const int s_state){
	dim_alphabet = dim_alf;
	num_states 	 = n_state;
	start_state  = s_state;

	ttable = new int*[num_states];
	for(int i=0; i<num_states; ++i)
		ttable[i] = new int[dim_alphabet+1];				// "+1" for Type column

	for(int i=0; i<num_states; ++i)
		for(int j=0; j<dim_alphabet+1; ++j){
			if(j >= dim_alphabet){
				ttable[i][j] = 0;
				continue;
			}
			ttable[i][j]=ND;
		}


	// Alphabet symbols
	alphabet=NULL;
	set_alphabet(alf, dim_alf);
}



// Delegating constructors
gi::dfa::dfa(const int n_state, const int dim_alf, const string *alf)
:dfa(n_state, dim_alf, alf, 0){}


gi::dfa::dfa(const int n_state, const int dim_alf, const string *alf, const int s_state, const int** tt_copy )
:dfa(n_state, dim_alf, alf, 0){
	for(int i=0; i<n_state; ++i)
		for(int j=0; j<dim_alf+1; ++j)
			ttable[i][j] = tt_copy[i][j];//d1.get(i,j);
}



// Copy constructor
gi::dfa::dfa(const dfa &d1)
:dfa(d1.num_states, d1.dim_alphabet, d1.alphabet, d1.start_state, (const int**) d1.ttable){}



gi::dfa::~dfa(){
	if(ttable != NULL){
		for(auto i=0; i<num_states; ++i)
			if(ttable[i] != NULL)
				delete[] ttable[i];

		delete [] ttable;
	}

	if(alphabet != NULL){
		delete[] alphabet;
	}

}



// Return a union dfa of (this) and "dfa_hp"
gi::dfa* gi::dfa::unionDFA(dfa* dfa_hp)
{
	// Number of states of union dfa
	int count_state = dfa_hp->num_states + num_states;


	// Union DFA instance
	dfa* union_dfa = new dfa(count_state, dim_alphabet, alphabet);


	// Definition of union DFA:
	// target dfa's states are arranged first (they have smaller indexes), then hypothesis dfa's states are added
	for(int j=0; j<num_states; ++j)												// Target automata
		for(int k=0; k<dim_alphabet+1; ++k)
			union_dfa->ttable[j][k] = ttable[j][k];


	for(int j=0; j<dfa_hp->num_states; ++j){									// HP dfa
		for(int k=0; k<dim_alphabet+1; ++k){									// In union dfa, start state of HP dfa takes index "num_state" of target dfa
			if(k != dim_alphabet)
				union_dfa->ttable[num_states+j][k] = dfa_hp->ttable[j][k] + num_states;
			else
				union_dfa->ttable[num_states+j][k] = dfa_hp->ttable[j][k];
		}
	}


	return union_dfa;
}



SYMBOL* gi::dfa::table_filling(){
	// The Table considered is only the upper triangular matrix, that can be saved in a linear array of size n(n-1)/2
	// Conversion of index form matrix to array are:
	//
	// From linear index k, to (i,j) for tha matrix (i row, j column)
	// i = n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
	// j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
	//
	// From (i,j) to k
	// Order such that i<j, then:
	// k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
	//
	// check (http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix)


	// *** TABLE-FILLING ALGORITHM with witness ***
	int i,j,k;
	int n = num_states;
	int tf_l = (num_states*(num_states-1))/2;


	// Table with pair of distinct states
	SYMBOL* distinti = new SYMBOL[tf_l];


	// Acceptor and rejector states are surely distinc. It exploits this to make a first distinction.
	for(i=0; i<(num_states-1); ++i)
		for(j=i+1; j<num_states; ++j){
			k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			// Check if a state is an acceptor and another is a rejector.
			if(ttable[i][dim_alphabet] != ttable[j][dim_alphabet]){
				distinti[k] = DFA_TF_STATE_X;
			}else
				distinti[k]=  DFA_TF_STATE_N;
		}


	// Minimizing loop
	// Each itaration check if the table was modified.
	bool modificata = true;
	while(modificata)
	{
		modificata = false;
		int stato_arrivo_1;
		int stato_arrivo_2;

		for(i=0; i<(num_states-1); ++i){
			for(j=i+1; j<num_states; ++j){
				k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
				//V cout << endl << "SP1: "<<i << ", SP2: "<<j;
				if(distinti[k] == DFA_TF_STATE_N){

					for(SYMBOL w=0; w<dim_alphabet; ++w)
					{
						stato_arrivo_1 = ttable[i][w];
						stato_arrivo_2 = ttable[j][w];

						// Lo stato di arrivo letto nel dfa potrebbe essere una coppia del tipo (2,2)
						if(stato_arrivo_1 == stato_arrivo_2)
							continue;

						// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
						// quindi le coppie devono avere sempre i<j
						if(stato_arrivo_2 < stato_arrivo_1){
							int tmp = stato_arrivo_1;
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						// Se finisco su me stesso vado avanti
						if(stato_arrivo_1 == i && stato_arrivo_2 == j)
							continue;

						//V cout << endl <<"SA1: "<<stato_arrivo_1 << ", SA2: "<<stato_arrivo_2 << " --> " << (int)distinti[stato_arrivo_1][stato_arrivo_2] << endl;
						// Se la coppia di arrivo è distinta, lo è anche quella originaria
						int i1 = stato_arrivo_1, j1 = stato_arrivo_2;
						int k1 = (n*(n-1)/2) - (n-i1)*((n-i1)-1)/2 + j1 - i1 - 1;

						if(distinti[k1] != DFA_TF_STATE_N){
							//TODO: è corretto?
							//distinti[k1] = w + '0'; //scrive l'indice (in char) della lettera dell'alfabeto per cui i due stati differiscono
							//distinti[k] = tableFillingState::TF_STATE_O;
							distinti[k] = w;
							modificata = true;
							break;											// Necessario!
						}
					}
				}
			}
		}
	}


	// Stampo la tabella degli stati distinti
	//    '@'  	-->	 It is an equivalent pair of states because its cell in the table remained empty
	// 	  'X'  	-->  It is a distinct pair of states
	// Otherwise ->	 It prints the symbol distincting the pair of states starting from initial state
	#ifdef DEBUG_DFA
	cout << "--------------------------" << endl;
	cout << "Tabella dei distinti " << endl;
	cout << "--------------------------" << endl;
	for(i=0; i<(num_states-1); ++i){
		for(j=i+1; j<num_states; ++j){
			k=(n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			char toprint= (distinti[k]==DFA_TF_STATE_X)?'X':(distinti[k]==DFA_TF_STATE_N)?'@': (char)(distinti[k] + 48);
			cout << "("<< i << "," << j << "):" << toprint << "  ";
		}
		cout << endl;
	}
	cout << "--------------------------" << endl;
	#endif

	return distinti;
}


// Minimization using Table-filing algorithmm
gi::dfa* gi::dfa::minimize_TF()
{
	const int num_row = num_states;

	// Matrix distinct: table for record distinct states
	bool** distinct = new bool*[num_states-1];								// Tabella delle coppie di stati distinti
	for(int i=0; i<num_row-1; ++i)											// EX:per 6 stati:5 righe e 6 colonne (di queste 6 solo 5 usate effettivamente, ma serve)
		distinct[i] = new bool[num_states];

	// Initialization of distinct table
	for(int i=0; i<num_states-1; ++i)
		for(int j=0; j<num_states; ++j)
			distinct[i][j]=false;

	// Distinguo tra accettanti e non accettanti							// EX: 0 <= i <= 4, j=i+1 (1 <= j <= 5 per la prima iterazione)
	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j)
			if(ttable[i][dim_alphabet] != ttable[j][dim_alphabet])				// Se uno è accettante mentre l'altro no
				distinct[i][j] = true;

	// Loop  minimizzante
	bool modificata = true;
	while(modificata)
	{
		modificata = false;

		for(int i=0; i<(num_states-1); ++i){
			for(int j=i+1; j<num_states; ++j){
				if(!distinct[i][j]){

					for(int k=0; k<dim_alphabet; ++k)
					{
						int stato_arrivo_1 = ttable[i][k];
						int stato_arrivo_2 = ttable[j][k];

						if(stato_arrivo_1 == stato_arrivo_2)
							continue;

						if(stato_arrivo_2 < stato_arrivo_1){				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
							int tmp = stato_arrivo_1;
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						if(distinct[stato_arrivo_1][stato_arrivo_2]){		// Se la coppia di arrivo è distinta, lo è anche quella originaria
							distinct[i][j] = true;
							modificata = true;
						}
					}
				}
			}
		}
	}


	// ** Mi preparo a costruire l'automa minimizzato **
	// Creo un elenco per avere una lista delle coppie di stati equivalenti più direttamente accessibile,
	// 	se è ND allora non ha un equivalente.
	int stato_equivalente[num_states];
	for(int i=0; i<num_states; ++i)
		stato_equivalente[i] = ND;

	// Conto il numero di stati finali, ed assegno ad ogni stato il suo eventuale stato equivalente
	// (devo anche controllare nel caso ho (0,4) (0,5) che quando incontro (4,5) non lo conto come un ulteriore stato in meno)
	int stati_finali = num_states;
	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j)
			if(!distinct[i][j] && stato_equivalente[i] == ND && stato_equivalente[j] == ND){
				stato_equivalente[j] = i;
				stati_finali--;
			}

	// Stampo le informazioni riguardo agli stati equivalenti trovati o meno
	/*//V cout << "N di stati finali: " << stati_finali << endl;
	cout << "Equivalenze stati: " << endl;
	for(int i=0; i<num_state; ++i)
		if(stato_equivalente[i] != ND)
			cout << "S:"<<i<<" --> "<<stato_equivalente[i]<<endl;
		else
			cout << "S:"<<i<<endl
	*/

	// STAMPO il vecchio DFA
	#ifdef DEBUG_2
	this->print_dfa_ttable("DFA BEFORE MINIMIZATION");
	#endif

	// Istanzio il nuovo DFA minimizzato
	dfa* dfa_min = new dfa(stati_finali,  dim_alphabet, alphabet, 0);

	int** ttable_min = dfa_min->get_ttable();

	int count = 0;
	for(int i=0; i<num_states; ++i){
		if(stato_equivalente[i] == ND){
			for(int j=0; j<dim_alphabet+1; ++j)
				ttable_min[count][j]=ttable[i][j];
			count++;
		}
	}

	// Aggiorno le transizioni verso stati ormai scomparsi perchè sostituiti con gli equivalenti
	int equivalenze_finora= 0;
	for(int i=0; i<num_states; ++i){
		if(stato_equivalente[i] != ND)
		{
			for(int k=0; k<stati_finali; ++k)
				for(int t=0; t<dim_alphabet; ++t)
					// Sostituisco la transizione allo stato "i" con lo stato equivalente "stato_equivalente[i]"
					if(ttable_min[k][t] == i)
						ttable_min[k][t] = stato_equivalente[i];
		}
	}

	// Qui aggiorno la label, perché se ad esempio ho collassato 2 stati prima dello stato 6,
	//  adesso lo stato 6 si trova nella riga 4, però le transizioni sono rimaste
	//  verso 6 e devono essere sistemate
	for(int i=0; i<num_states; ++i)
	{
		if(stato_equivalente[i] != ND)
			equivalenze_finora++;

		if(equivalenze_finora != 0)
		{
			int nuova_label = i-equivalenze_finora;
			for(int k=0; k<stati_finali; ++k)
				for(int t=0; t<dim_alphabet; ++t)
					if(ttable_min[k][t] == i)
						ttable_min[k][t] = nuova_label;
		}
	}

	// STAMPO il dfa minimizzato
	#ifdef DEBUG_2
	dfa_min->print_dfa_ttable("MINIMIZED DFA");
	#endif

	//libero la memoria allocata
	if(distinct != NULL){
		for(int i=0; i<num_row -1; ++i){
			if(distinct[i] != NULL){ 
				delete[] distinct[i];
			}
		}
		delete [] distinct;
	}

	return dfa_min;
}


void gi::dfa::print_dfa_ttable(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "  ";
	for(int i=0; i<dim_alphabet; ++i)
		header = header + " | "+ intTostring(i);
	header = header + " - A";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)			// Valori delle transizioni, o ND o il valore
				cout << " N ";
			else if(j < dim_alphabet)
				cout << " "<< ttable[i][j] <<"  ";

			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
			{
				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
					cout << "  / ";
				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
					cout << " Ac ";
				else if(ttable[i][j] == DFA_STATE_REJECTING)
					cout << " Ri ";
				else
					cout << "  X ";
			}
		}
		cout << endl;
	}

	cout << "--------------------------" << endl;

}


void gi::dfa::print_dfa_ttable_mapped_alphabet(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
	// Uso l'alfabeto originale


	// Using Mapped Alphabet

	/*typedef	map<char, char>::const_iterator It;
	cout << "Mappatura inversa correta? "<<endl;
	for(It p1=inverse_mapped_alphabet.begin(); p1!=inverse_mapped_alphabet.end(); ++p1)
		cout << (*p1).first << "; " << (*p1).second << endl;
	cout << "FINE"<<endl;*/

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "  ";
	for(int i=0; i<dim_alphabet; ++i)
		header = header + " | "+ alphabet[i];
	header = header + " - A";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)			// Valori delle transizioni, o ND o il valore
				cout << " N ";
			else if(j < dim_alphabet)
				cout << " "<< alphabet[ttable[i][j]] <<"  ";

			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
			{
				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
					cout << "  / ";
				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
					cout << " Ac ";
				else if(ttable[i][j] == DFA_STATE_REJECTING)
					cout << " Ri ";
				else
					cout << "  X ";
			}
		}
		cout << endl;
	}

	cout << "--------------------------" << endl;

}


void gi::dfa::print_dfa_in_text_file(const string file_path)
{

	for(int i=0; i<num_states; ++i)
		for(int j=0; j<dim_alphabet+1; ++j)
			if(ttable[i][j] == ND){
				cout << "ERROR: DFA is not completely defined, there ND transition" << endl;
				exit(EXIT_FAILURE);	//TODO: exception
			}



	// **********************
	//   WRITE DFA in FILE
	// **********************

	// Opena file
	ofstream myfile;
	myfile.open(file_path.c_str());


	// Write alphabet size
	myfile << intTostring(dim_alphabet) << " ";

	// Write num of states
	myfile << intTostring(num_states) << " ";

	// Write dfa name
	myfile << "dfa" << "\n";

	// Write alphabet symbols
	for(int i=0; i<dim_alphabet; ++i)
		myfile << alphabet[i] << " ";
	myfile << "\n";


	// Write transition table
	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet+1; ++j)
		{
			myfile << "dfa[" <<intTostring(i)<<"][";
			if(j < dim_alphabet)
				myfile << alphabet[j] << "]="<< ttable[i][j] <<";\n";

			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
			{

				myfile << intTostring(dim_alphabet) << "]=";
				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
					myfile << "0";
				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
					myfile << "1";
				else if(ttable[i][j] == DFA_STATE_REJECTING)
					myfile << "0";

				myfile <<";\n";
			}
		}
	}


	myfile.close();
}




//void gi::dfa::print_dfa_with_color(string title)
//{
//	// STAMPA IL DFA
//	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
//
//	cout << endl<< "--------------------------" << endl;
//	cout << title << endl;
//	string header = "    ";
//		for(int i=0; i<dim_alphabet; ++i)
//			header = header + " | "+ intTostring(i);
//		header = header + " - A  - C";
//	cout << header << endl;
//
//	for(int i=0; i<num_state; ++i){
//		cout << "S"<<i<<"  ";
//
//		for(int j=0; j<dim_alphabet+1; ++j)
//		{
//
//			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
//				cout << " N ";
//
//			else if(j < dim_alphabet)
//				cout << " "<< ttable[i][j] <<" ";
//
//			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
//			{
//				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
//					cout << "  / ";
//				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
//					cout << " Ac ";
//				else if(ttable[i][j] == DFA_STATE_REJECTING)
//					cout << " Ri ";
//				else
//					cout << "  X ";
//			}
//		}
//
////		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
////			cout << " W ";
////		else if(this->is_inside_blue_states(i))
////			cout << " B";
////		else
////			cout << " R";
//
//		cout << endl;
//	}
//
//	cout << "--------------------------" << endl;
//}
//
//void gi::dfa::print_dfa_with_color_mapped_alphabet(string title)
//{
//	// STAMPA IL DFA
//	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
//
//	cout << endl<< "--------------------------" << endl;
//	cout << title << endl;
//	string header = "    ";
//		for(int i=0; i<dim_alphabet; ++i)
//			header = header + " | "+ alphabet[i];
//		header = header + " - A  - C";
//	cout << header << endl;
//
//	for(int i=0; i<num_state; ++i){
//		cout << "S"<<i<<"  ";
//
//		for(int j=0; j<dim_alphabet+1; ++j)
//		{
//
//			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
//				cout << " N ";
//			else if(j < dim_alphabet)
//				cout << " "<< alphabet[ttable[i][j]] <<" ";
//
//			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
//			{
//				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
//					cout << "  / ";
//				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
//					cout << " Ac ";
//				else if(ttable[i][j] == DFA_STATE_REJECTING)
//					cout << " Ri ";
//				else
//					cout << "  X ";
//			}
//		}
//
////		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
////			cout << " W ";
////		else if(this->is_inside_blue_states(i))
////			cout << " B";
////		else
////			cout << " R";
//
//		cout << endl;
//	}
//
//	cout << "--------------------------" << endl;
//
//}


void gi::dfa::set_ttable(int** ext_ttable){
	ttable = ext_ttable;
}

int** gi::dfa::get_ttable(){
	return ttable;
}

int gi::dfa::get_ttable(int i, int j){

	if(i<num_states && j<dim_alphabet+1)
		return ttable[i][j];
	else{
		cerr<<"dfa::get_ttable: out of bound"<<endl;
		exit(EXIT_FAILURE);
	}

}

void gi::dfa::set_ttable_entry(int i, int j, int v){
	if(i<num_states && j<dim_alphabet+1)
		ttable[i][j]=v;
	else{
			cerr<<"dfa::set_ttable: out of bound"<<endl;
			exit(EXIT_FAILURE);
		}
}


void	gi::dfa::set_acceptor_state(int state)
{
	ttable[state][dim_alphabet] = DFA_STATE_ACCEPTING;
}


void	gi::dfa::set_rejector_state(int state)
{
	ttable[state][dim_alphabet] = DFA_STATE_REJECTING;
}


int gi::dfa::get_dim_alphabet(){
	return dim_alphabet;
}

int gi::dfa::get_num_states() const{
	return num_states;
}

int gi::dfa::get_start_state(){
	return start_state;
}


void gi::dfa::print_dfa_dot(string title, const char *file_path)
{
	ofstream myfile;
	myfile.open(file_path);

	string header = "digraph "+title+" {\n";
	string start_state = "__start0 [label=\"\" shape=\"none\"];\n\n";

	start_state = start_state + "rankdir=LR;\nsize=\"8,5\";\n\n";

	//States
	string states = "";
	string shape = "";
	string style="";
	string color="";
	for(int i=0; i<num_states; ++i)
	{
		if(ttable[i][dim_alphabet] == DFA_STATE_UNREACHABLE)
			continue;

		if(ttable[i][dim_alphabet] == DFA_STATE_ACCEPTING){
			shape = "doublecircle";
			style = "rounded,filled";
		}
		else if(ttable[i][dim_alphabet] == DFA_STATE_REJECTING){
			shape = "circle";
			style = "filled";
		} else {
			shape = "circle";
			style = "filled";
		}

		color="white";

		states = states + "s"+intTostring(i)+" [style=\""+style+"\", color=\"black\", fillcolor=\""+color+"\" shape=\""+shape+"\", label=\""+intTostring(i)+"\"];\n";
	}

	// Transizioni
	string transitions = "";
	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet; ++j){
			int arrive_state = ttable[i][j];
			if(arrive_state == ND)
				continue;

			transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(arrive_state)+" [label=\""+	intTostring(j)+"\"];\n";
		}
	}

	string end = "__start0 -> 0;";
	string footer ="\n}";

	myfile << header << start_state << states << transitions /*<< end*/<<footer;

	myfile.close();
}



void gi::dfa::print_dfa_dot_mapped_alphabet(string title, const char *file_path)
{
	string state_name_prefix = "q";
	ofstream myfile;

	myfile.open(file_path);


	string header = "digraph "+title+" {\n";

	string start_state = "__start0 [style = invis, shape = none, label = \"\", width = 0, height = 0];\n\n";

	start_state = start_state + "rankdir=LR;\nsize=\"8,5\";\n\n";

	string start_arrow = "";
	start_arrow = "subgraph cluster_main { \n\tgraph [pad=\".75\", ranksep=\"0.15\", nodesep=\"0.15\"];\n\t style=invis; \n\t__start0 -> s0 [penwidth=2];\n}\n";

	//States
	string states = "";
	string shape = "";
	string style="";
	string color="";
	for(int i=0; i<num_states; ++i)
	{
		if(ttable[i][dim_alphabet] == DFA_STATE_UNREACHABLE)
			continue;

		if(ttable[i][dim_alphabet] == DFA_STATE_ACCEPTING){
			shape = "doublecircle";
			style = "rounded,filled";
		}
		else if(ttable[i][dim_alphabet] == DFA_STATE_REJECTING){
			shape = "circle";
			style = "filled";
		} else {
			shape = "circle";
			style = "filled";
		}

		color="white";

		states = states + "s"+intTostring(i)+" [style=\""+style+"\", color=\"black\", fillcolor=\""+color+"\" shape=\""+shape+"\", label=\""+state_name_prefix+intTostring(i)+"\"];\n";
	}


	// Transizioni
	string transitions = "";
	//map<string, string> label_for_transiction;					// La chiave individua una coppia di stati tra cui potrebbe esserci una transizione
																	// Il valore è la label da stampare, contenente tutti i simboli per cui avviene quella transizione

	vector< vector<string> > label_for_transiction(num_states, vector<string>(num_states));

	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet; ++j){

			int arrive_state = ttable[i][j];

			if(arrive_state == ND)
				continue;

			string transition_symbol = alphabet[j];

			if(label_for_transiction[i][arrive_state].length() == 0)
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + transition_symbol;
			else if(label_for_transiction[i][arrive_state].length() % 9 == 0)			// Inserisce ogni 7 simboli un ritorno a capo nella label
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + "\\n" + transition_symbol;
			else
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + "," +transition_symbol;


			// ORIGINALE un carattere - una transizione
			//transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(arrive_state)+" [label=\""+inverse_mapped_alphabet[alphabet_symbols[j]]+"\"];\n";
		}
	}


	for(int i=0; i<num_states; ++i)
		for(int j=0; j<num_states; ++j){
			if(label_for_transiction[i][j].compare(""))
				transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(j)+" [label=\""+label_for_transiction[i][j]+"\"];\n";
		}

	string end = "__start0 -> 0;";
	string footer ="\n}";

	myfile << header << start_state <<  states << start_arrow << transitions /*<< end*/<<footer;


	myfile.close();
}


void gi::dfa::set_num_state(int n){
	num_states = n;
}


void gi::dfa::set_alphabet(const string* alf, const int d_alf)
{
	// Erase the existing alphabet
	if(alphabet!=NULL) delete[] alphabet;

	// Set new size of alphabet
	dim_alphabet = d_alf;

	// Instance for new alphabet
	alphabet = new string[dim_alphabet];

	// Copy alphabet in input to alphabet for current dfa
	// inside mapped_alphabet for every symbol there is the associated index
	for(SYMBOL i=0; i<dim_alphabet; ++i){
		alphabet[i] = alf[i];
		mapped_alphabet[(string)alf[i]] = i;
	}

}


string* gi::dfa::get_alphabet(){
	return alphabet;
}


int gi::dfa::get_arrive_state(vector<SYMBOL> dfa_string)
{
	int state = 0;
	int next_state=ND;
	int lettera = 0;

	//typedef	vector<SYMBOL>::const_iterator It_posneg;

	//for(It_posneg i=dfa_string.begin(); i!=dfa_string.end(); ++i){
	for(auto i=dfa_string.begin(); i!=dfa_string.end(); ++i){
		next_state = ttable[state][(*i)];
		if(next_state == ND){
			state = ND;
			break;
		}
		state = next_state;
	}

	return state;		// Ritorna ND se la stringa non è compatibile, diversamente ritorna lo stato dove termina la stringa
}


gi::dfa gi::dfa::read_dfa_file(const string file_name)
{
	char nameDFA[BUFFER_SIZE];
	char line[BUFFER_SIZE];

	string *alphabet_file=NULL;

	int counter = 0;
	int num_total_line = 0;
	int cstato = 0;
	char calfabeto[10];
	int ctransizione = 0;
	int current_line = 0;
	int start_state = 0;

	string n;

	ifstream read;
	string template_line;

	dfa res;

	// Open connection to file
	read.open(file_name.c_str());

	if(read.is_open())
		cout << "File " << file_name << " is open."<< endl;
	else{
		cerr << "Error opening " << file_name << endl;
		exit(EXIT_FAILURE);
	}

	// initial state
	start_state = 0;

	// Read first line and set "num states", "dim alf" and "dfa name"
	read.getline(line,BUFFER_SIZE);
	counter = sscanf(line, "%d %d %s", &(res.dim_alphabet), &(res.num_states), nameDFA);

	// Check if the first line is complete
	if(counter != 3){
		cout << "Error in first line of file" << endl;

		exit(EXIT_FAILURE);
	}


	/////////////////////////  ALPHABET  //////////////////////////////////
	// read the alphabet
	alphabet_file = new string[res.dim_alphabet];
	read.getline(line,BUFFER_SIZE);

	istringstream iss(line);

	//cout << "DImensione alfabeto letta: "<< intTostring(res.dim_alphabet) << endl;
	counter=0;
	while (iss >> n){
		if(counter >= res.dim_alphabet)
			break;

		alphabet_file[counter++] = n;
	}



	// check read alphabet
	if(counter != res.dim_alphabet){
		cout << "Error in second line of file: issue with size of alphabet" << endl;

		if(alphabet_file)
			delete [] alphabet_file;

		exit(EXIT_FAILURE);
	}


	// Set alphabet and "mapped_alphabet" for the current dfa
	res.set_alphabet(alphabet_file, res.dim_alphabet);

	if(alphabet_file)
		delete [] alphabet_file;




	///////////// compute utility values//////////////////
	// expected lines of file
	num_total_line = (res.dim_alphabet+1)*(res.num_states);

	// template of the all lines of the file
	template_line = (string)nameDFA+"[%d][%[^]]] = %d;";
	/////////////////////////////////////////////////////



	// allocate memory for ttable
	res.ttable = new int*[res.num_states];

	// initialize ttable rows
	// "+2" is for algorithm like EDSM with states Type and Colour
	for(int i=0; i<res.num_states; ++i){
		res.ttable[i] = new int[res.dim_alphabet+1];
		for(int j=0; j<res.dim_alphabet+1; ++j)
			res.ttable[i][j]=ND;
	}


	// Parsing the file
	while(!read.eof())
	{
		read.getline(line,BUFFER_SIZE);

		string cline = line;

		// Handler for last line
		string trimmedline = trim(cline);
		if(trimmedline == "")										// Happen only in the last line
			continue;

		++current_line;

		// Integrity check
		if(current_line > num_total_line){

			cerr << "ERROR - Line number greater than max" << endl;
			exit(EXIT_FAILURE);	// TODO: exapection
		}

		// Read line and set transition
		counter = sscanf(line, template_line.c_str(), &cstato, &calfabeto, &ctransizione);
		if(counter != 3)
			cerr << "ERROR in current line"<<current_line<<endl;

		string transition_symbol = calfabeto;


		// It detects the row for type of state (accepting or rejecting)
		if(transition_symbol.compare(intTostring(res.dim_alphabet)) == 0)
			res.ttable[cstato][stringToint(transition_symbol)] = ctransizione;
		else
			res.ttable[cstato][res.mapped_alphabet[transition_symbol]] = ctransizione;

	}


	//res.print_dfa_ttable("TABELLA");

	// Close connection
	read.close();

	return res;
}


//void gi::dfa::read_example_file(string path_samples, vector<SYMBOL>* &positive, int* dim_positive, vector<SYMBOL>* &negative,  int *dim_negative, int* &wp = NULL, int* &wn = NULL)
//{
//	cout << "Reading strings from txt file: "<<endl;
//	int cp = 0;														// Numero di stringhe positive del linguaggio
//	int cn = 0;														//   -    -     -     negative  -      -
//	char ch;
//
//	cout << path_samples << endl;
//
//	fstream fin(path_samples.c_str(), fstream::in);
//
//	if(!fin){
//		cerr<<"An error occurred in opening file :("<<endl;
//		exit(EXIT_FAILURE);
//	}
//
//	while (fin >> noskipws >> ch) {									// Faccio in conteggio previo del numero di stringhe positive e negative presenti nel txt
//		if(ch == '\n')
//			continue;
//		else if(ch == '+')
//			cp++;
//		else if(ch == '-')
//			cn++;
//	}
//	(*dim_positive) = cp;
//	(*dim_negative) = cn;
//
//	positive = new vector<SYMBOL>[cp];
//	negative = new vector<SYMBOL>[cn];
//	wp	= new int[cp];
//	wn	= new int[cn];
//
//	for(int i=0; i<cp; ++i)
//		wp[i] = 0;
//	for(int i=0; i<cn; ++j)
//		wn[i] = 0;
//
//
//	cout << intTostring(cp) + " positivi"<<endl;
//	cout << intTostring(cn) + " negativi"<<endl;
//
//	int flagcp = 0;
//	int flagcn = 0;
//	bool casopositive = false;
//	bool casonegative = false;
//	bool primap = true;
//	bool priman = true;
//
//	ifstream infile(path_samples.c_str());
//
//	bool first = true;
//	bool second = false;
//	string line;
//
//	int local_dim_alphabet;
//
//	while (getline(infile, line))
//	{
//	    istringstream iss(line);
//	    int a;
//	    string n;
//
//
//	    // Read first line for dim alphabet
//	    if(first){
//	    	if (!(iss >> a)) { break; } // error
//	    	local_dim_alphabet = a;
//	    	//cout << "dimensione alfabeto " << a << endl;
//	    	first = false;
//	    	second = true;
//
//	    	continue;
//	    }
//
//	    if(second){
//
//	    	int counter=0;
//	    	while (iss >> n){
//	    		if(counter>=local_dim_alphabet)
//	    			break;
//
//	    		if(mapped_alphabet.find(n[0])!=mapped_alphabet.end()){
//	    			cerr<<"Error in reading example: alphabet different from that of dfa!"<<endl;
//
//	    			exit(EXIT_FAILURE);
//	    		}
//	    	}
//
//	    	// Alphabet
//	    	if(counter!= local_dim_alphabet){
//	    		cerr<<"Error in reading example: number of red alphabet symbols mismatches with the declared one!"<<endl;
//	    		cerr<<"Expected symbols: "<<dim_alphabet<<endl;
//	    		cerr<<"Red symbols: "<<counter<<endl;
//
//
//	    		exit(EXIT_FAILURE);
//	    	}else if(local_dim_alphabet != dim_alphabet){
//	    		cerr<<"Error in reading example: alphabet different from that of dfa!"<<endl;
//
//	    		exit(EXIT_FAILURE);
//	    	}
//
//	    	// alphabet ok ;)
//	    	second= false;
//	    }
//
//	    bool weight = true;
//
//	    // Read lines
//		while (iss >> n)
//		{
//			if( !n.compare("+") ){
//
//				weight = true;
//				casopositive = true;
//				casonegative = false;
//				if(primap){												// Se è il primo caso evito l'incremento
//					primap =false;
//					continue;
//				}
//				flagcp++;
//				continue;
//
//			}else if( !n.compare("-") ){
//				weight = true;
//				casonegative = true;
//				casopositive = false;
//				if(priman){												// Se è il primo caso evito l'incremento
//					priman =false;
//					continue;
//				}
//				flagcn++;
//				continue;
//			}
//
//
//			int tmp = mapped_alphabet[(char) n[0]];
//
//			if(weight){
//				weight = false;
//
//				if(casopositive)
//					wp[flagcp] = stringToint(n);
//				else if(casonegative)
//					wn[flagcn] = stringToint(n);
//
//			} else {
//
//				if(casopositive)
//							positive[flagcp].push_back(tmp);
//				else if(casonegative)
//							negative[flagcn].push_back(tmp);
//			}
//		}
//	}
//}


vector<SYMBOL> gi::dfa::equivalence_query(dfa* dfa_hp){
	vector<SYMBOL> witness;

	#ifdef DEBUG_2
	cout << endl << "--------------------------" << endl;
	cout << "EQUIVALENCE QUERY" << endl;
	cout << "--------------------------" << endl;
	#endif


	// Build union DFA of target dfa (thisone) and dfa_hp
	dfa* dfa_union = this->unionDFA(dfa_hp);
	#ifdef DEBUG_2
	dfa_union->print_dfa_ttable("DFA UNIONE");
	#endif


	// Table-filling algorithm on union dfa
	SYMBOL* distincts_table = dfa_union->table_filling();

	// Extract list of equivalent states from table of distinct states,
	// every vector contain a list of equivalent states for the state that correspond to the vector.
	vector<int>* equivalent_states_list = dfa_union->equivalent_states_list_from_table(distincts_table);

	// Verify if start states of dfas are equivalent:
	// Controllo se fra gli stati equivalenti allo stato 0 (stato iniziale dell'automa corrente)
	// c'è lo stato iniziale dell'automa ipotesi identificato dall'indice "num_state".
	// Nel caso fosse presente, non entro nell'if, e ritorno un vector vuoto come controesempio.
	// Se così non fosse (è il caso dell'"end()") genero un controesempio.
	if(equivalent_states_list[0].end() == find(equivalent_states_list[0].begin(), equivalent_states_list[0].end(), num_states) )
		witness = dfa_union->witness_from_table(distincts_table, num_states);


	// Free allocated memory
	if(distincts_table != NULL)
		delete[] distincts_table;

	delete [] equivalent_states_list;

	delete dfa_union;


	return witness;
}


vector<vector<SYMBOL>> gi::dfa::get_characterization_set()
{
	// Characterization set of examples for current DFA
	vector<vector<SYMBOL>> characterization_set(0, vector<SYMBOL>());


	// Table-filling algorithm over union dfa
	SYMBOL* distincts_table = this->table_filling();


	// Extract list of equivalent states from table of distinct states,
	// every vector contain a list of equivalent states for the state that correspond to the vector.
	//vector<int>* equivalent_states_list = this->equivalent_states_list_from_table(distincts_table);

	// Check if the curret automaton is minimal
	dfa* tmp_minimized_dfa = this->minimize_TF();
	if(tmp_minimized_dfa->get_num_states() != num_states){
		cout << "ERROR: Processed DFA is not minimial! Minimize it" << endl;
		delete tmp_minimized_dfa;
		exit(EXIT_FAILURE);
	}
	delete tmp_minimized_dfa;


	// Extract access strings for the DFA states
	map<int, vector<SYMBOL>> access_strings = this->get_access_strings();


	// track pairs of states already checked
	map<int, vector<int>> state_pairs;
	//for(int i=0; i<num_states; ++i)
	//	state_pairs[i] = vector<int>();


	// Per ogni coppia di stati devo creare una stringa che adopera il simbolo presente nella tabella del Table-filling
	// questo è il simbolo che contraddistingue i due stati in considerazione
	for(int i=0; i<num_states; ++i)
	{
		for(int j=i+1; j<num_states; ++j)
		{
			//V cout << "i: "<<i<<"; j: "<<j<<endl;

			// If current pair was checked yet, it goes over
			if( std::find(state_pairs[i].begin(), state_pairs[i].end(), j) != state_pairs[i].end() )
				continue;

			// Generated witness for the pairs of states
			vector<SYMBOL> wit;

			// Reed symbol during the execution
			SYMBOL reed_symbols;

			int icoppia = i;
			int jcoppia = j;


			while(1)
			{
				//Vcout << "Coppia: "<<icoppia << ","<<jcoppia<<";" <<endl;

				// Lo stato di arrivo letto nel dfa potrebbe essere una coppia del tipo (2,2)
//				if(icoppia == jcoppia){
//					cout << "NB: Smaller string" << endl;
//					break;
//				}

				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
				// quindi le coppie devono avere sempre i<j
				if(jcoppia < icoppia){
					int tmp = icoppia;
					icoppia = jcoppia;
					jcoppia = tmp;
				}


				// Add checked state
				state_pairs[icoppia].push_back(jcoppia);

				int n = num_states;
				int k = (n*(n-1)/2) - (n-icoppia)*((n-icoppia)-1)/2 + jcoppia - icoppia - 1;
				//Vcout << "Coppia: "<<icoppia << ","<<jcoppia<<"; vale: "<< distincts_table[k] <<endl;


				if(distincts_table[k] == DFA_TF_STATE_N){
					cout << "ERROR! Required conuterexample with equivalent states" << endl;
					exit(EXIT_FAILURE);
				}

				// Reed the simbol which specific the difference for the pair of states
				reed_symbols = distincts_table[k];

				//Vcout << "Reed: "<< reed_symbols << "; SIMBOL: "<<DFA_TF_STATE_X << endl;

				// Caso in cui lo start state dei 2, uno è accettante e l'altro no
				if(reed_symbols == DFA_TF_STATE_X){
					//Vcout << "Entrato!"<<endl;
					break;
				}

				//V cout << "stinga parziale: "<<wit<<" + " << input << endl;
				wit.push_back(reed_symbols);
				//V cout << "Vale: "<< input << ", ascii: "<< (int)input << ", diff_i:" << (int)input - 48 << endl;

				icoppia = ttable[icoppia][reed_symbols];
				jcoppia = ttable[jcoppia][reed_symbols];
			}



			// Add generated strings to characterization set
			if(wit.size() != 0)
			{
				// Prefix of characterizing strings is the access strings for analyzed states
				vector<SYMBOL> first_characterizing_strings = access_strings[i];
				vector<SYMBOL> second_characterizing_strings = access_strings[j];

				// Buil characterzing strings concatening access strings with witness for analyzed states
				first_characterizing_strings.insert( first_characterizing_strings.end(), wit.begin(), wit.end() );
				second_characterizing_strings.insert( second_characterizing_strings.end(), wit.begin(), wit.end() );

				// Check if current sample is in the set yet
				if( std::find(characterization_set.begin(), characterization_set.end(), first_characterizing_strings) == characterization_set.end() )
					characterization_set.push_back(first_characterizing_strings);
				if( std::find(characterization_set.begin(), characterization_set.end(), second_characterizing_strings) == characterization_set.end() )
					characterization_set.push_back(second_characterizing_strings);

				//Vcout << "Controesempio A aggiunto è: "<< first_characterizing_strings << endl;
				//Vcout << "Controesempio B aggiunto è: "<< second_characterizing_strings << endl;
			}

		}
	}

	return characterization_set;

}


vector<vector<SYMBOL>> gi::dfa::get_cover_set()
{

	//// Support structers
	// Record the already visited nodes
	vector<int> 							   	visited_nodes;

	// Queue of nodes to be checked
	list<int>										queue;

	// Structure to record access strings
	map<int, vector<SYMBOL>> 	access_strings;

	// Cover set
	vector<vector<SYMBOL>> 		cover_set(1, vector<SYMBOL>());  // "1" is for epsilon transition


	//// Init
	// Insert as first state the start state
	queue.push_back(start_state);
	visited_nodes.push_back(start_state);
	int current_node = start_state;


	//// BFS
	while(!queue.empty())
	{
		// Reference to the front of the queue
		current_node = queue.front();
		queue.pop_front();

		//cout << "Father node: "<< intTostring(current_node) << endl;

		// Cycle on linked node to current node
		for(int i=0; i<dim_alphabet; ++i)
		{
			int child_node = ttable[current_node][i];
			//cout << "Child node: "<< intTostring(child_node) << endl;


			// Add string for transition towards child nodes
			vector<SYMBOL> child_access_string = access_strings[current_node];
			child_access_string.push_back(i);
			cover_set.push_back(child_access_string);


			// If it is a not visited nodes
			if( std::find(visited_nodes.begin(), visited_nodes.end(), child_node) == visited_nodes.end() )
			{
				// Current access string
				vector<SYMBOL> current_access_string;

				if( access_strings.find(current_node) != access_strings.end()  )
					current_access_string = access_strings[current_node];
				current_access_string.push_back(i);
						//cout << "Stringa per il padre: "<<access_strings[current_node] << endl;
						//cout << "Stringa di accesso corrente: "<< current_access_string << endl;


				// Se non esiste la entry nella tabella la inserisco, se no confronto la lungzza, se è più corta la inserisco
				if (access_strings[child_node].empty())
					access_strings[child_node] = current_access_string;
				else if(access_strings[child_node].size() > current_access_string.size() )
					access_strings[child_node] = current_access_string;


				// Insert analyzed child node into queue and visited node set
				visited_nodes.push_back(child_node);
				queue.push_back(child_node);
					//cout << "Nodo "<< intTostring(child_node) << " inserito nella coda"<<endl;
			}
		}

	}


	return cover_set;
}



// It returns the cover set with strings concatenated to necessary prefix
// with regard to the difference of number of states among hypothesis and
// target automata
vector<vector<SYMBOL>> gi::dfa::get_augmented_characterization_set(int num_of_states_target_dfa, vector<vector<SYMBOL>>& aug_characterization_set)
{

	if(num_of_states_target_dfa == 0){
		cout << "ERR: Target DFA with no states" << endl;
		exit(EXIT_FAILURE);
	}

	int diff_of_states = num_states - num_of_states_target_dfa;

	//Vcout << "Differenza: "<<diff_of_states << endl;


	// Get simple characterization_set
	cout << "... START Simple char set creation..." << flush;
	vector<vector<SYMBOL>> characterization_set = get_characterization_set();
	cout << "END Simple char set. Size: "<< characterization_set.size() << flush;

	// New augmented_characterization_set
	//vector<vector<SYMBOL>> aug_characterization_set;

	// Prefixes to be concatened
	vector<vector<SYMBOL>> prefixes;



	///INIT
	for(int j=0; j<dim_alphabet; ++j){
		vector<SYMBOL> tmp;
		tmp.push_back(j);
		prefixes.push_back(tmp);
	}



	// Hypothesis have more states than target automaton
	if(diff_of_states != 0)
	{
		// If diff_of_states >0
		int prefix_length = diff_of_states;

		// If diff_of_states < 0
		if(diff_of_states < 0)
			prefix_length = 0;


		cout << "START creating prefixes to be concatenated..." << flush;
		// Create prefixes
		for(int i=0; i<prefix_length; ++i)
		{
			// Compute a new subset of strings
			vector<vector<SYMBOL>> new_subset;
			for(int j=0; j<dim_alphabet; ++j)
			{
				for(auto &it : prefixes)
				{
					vector<SYMBOL> tmp_prefix = it;
					tmp_prefix.push_back(j);
					new_subset.push_back(tmp_prefix);
					if(new_subset.size() > LIMIT_OF_TESTSET_W_METHOD)
						throw "TEST_SET_TOO_BIG";
				}
			}

			// It adds the new calculated sub set of strings at "prefixes"
			for(auto &it : new_subset)
				prefixes.push_back(it);
		}
		cout << "END prefixes to be concatenated. Size: " << prefixes.size() << endl;




		cout << "START creating final aug char set..." << flush;
		// It concatenates prefixes to the simple characterization_set
		for(auto &it1 : characterization_set)
		{
			aug_characterization_set.push_back(it1);
			//Vcout << endl << "Letto1: "<< it1 << endl;

			for(auto &it2 : prefixes){
				vector<SYMBOL> new_string = it2;

				new_string.insert( new_string.end(), it1.begin(), it1.end());

				aug_characterization_set.push_back(new_string);

				//Vcout << "Letto2: "<< it2 << endl;
				//Vcout << "Aggiunto: "<< new_string << endl;
			}
			//cout << "Current size of aug char set: "<<aug_characterization_set.size() << flush;
			if( aug_characterization_set.size() > LIMIT_OF_TESTSET_W_METHOD)
				throw "TEST_SET_TOO_BIG";

		}
		cout << "END final aug char set. Size: " << aug_characterization_set.size() << flush;


	} else {
		// When diff_of_states is 0 return the simple characterization set
		return characterization_set;
	}


	return aug_characterization_set;
}



// It returns the cover set with strings concatenated to necessary prefix
// with regard to the difference of number of states among hypothesis and
// target automata
// from
vector<vector<SYMBOL>> gi::dfa::get_w_method_test_set(int num_of_states_target_dfa)
{

	// Final test set
	vector<vector<SYMBOL>> w_method_test_set;


	cout << "START Cover set...";


	// (Transition) cover set - P
	vector<vector<SYMBOL>> cover_set = get_cover_set();


	cout << "END. Size: "<< cover_set.size() << endl;
	cout << "START aug char set ... " << flush;



	// Augmented characterization set - Z
	vector<vector<SYMBOL>> aug_characterization_set;

	try
	{
		get_augmented_characterization_set(num_of_states_target_dfa, aug_characterization_set);
	}catch(const std::bad_alloc&){
		cerr << "ERR: Too memory allocated" << endl;
		throw "TEST_SET_TOO_BIG";
	}

	cout << "END. Size: " << aug_characterization_set.size() << endl;



	// Limit the number of strings
	if(cover_set.size() * aug_characterization_set.size() > 100000000)
		throw "TEST_SET_TOO_BIG";



	cout << "START Test set...";
	////// Compute test set: concatenating cover set and aug. charac. set
	for(auto &it1 : cover_set)
	{
		// Add non-cancatened string of cover set to final set
		if(it1.size() != 0)
			w_method_test_set.push_back(it1);


		// Add concatened strings
		for(auto &it2 : aug_characterization_set)
		{
			vector<SYMBOL> tmp_string = it1;

			// Concatenating strings
			tmp_string.insert( tmp_string.end(), it2.begin(), it2.end() );

			// Add string to final set
			w_method_test_set.push_back(tmp_string);
		}
	}

	cout << "END. Size:" << w_method_test_set.size() << endl;

	return w_method_test_set;
}





map<int, vector<SYMBOL>>  gi::dfa::get_access_strings()
{

	//// Support structers
	// Record the already visited nodes
	vector<int> 							   	visited_nodes;

	// Queue of nodes to be checked
	list<int>										queue;

	// Structure to record access strings
	map<int, vector<SYMBOL>> 	access_strings;
	//map<int, string> 						mapped_access_strings;


	//// Init
	// Insert as first state the start state
	queue.push_back(start_state);
	visited_nodes.push_back(start_state);
	int current_node = start_state;


	// BFS
	while(!queue.empty())
	{
		// Reference to the front of the queue
		current_node = queue.front();
		queue.pop_front();

		//cout << "Father node: "<< intTostring(current_node) << endl;


		// Cycle on linked node to current node
		for(int i=0; i<dim_alphabet; ++i)
		{
			int child_node = ttable[current_node][i];
			//cout << "Child node: "<< intTostring(child_node) << endl;


			// If it is a not visited nodes
			if( std::find(visited_nodes.begin(), visited_nodes.end(), child_node) == visited_nodes.end() )
			{
				// Current access string
				vector<SYMBOL> current_access_string;

				if( access_strings.find(current_node) != access_strings.end()  )
					current_access_string = access_strings[current_node];
				current_access_string.push_back(i);
						//cout << "Stringa per il padre: "<<access_strings[current_node] << endl;
						//cout << "Stringa di accesso corrente: "<< current_access_string << endl;


				// Se non esiste la entry nella tabella la inserisco, se no confronto la lungzza, se è più corta la inserisco
				if (access_strings[child_node].empty())
					access_strings[child_node] = current_access_string;
				else if(access_strings[child_node].size() > current_access_string.size() ){
						//cout << "vecchia stringa di accesso: "<<access_strings[child_node] << " aggiornata!"<<endl;
					access_strings[child_node] = current_access_string;
				}


				// Insert analyzed child node into queue and visited node set
				visited_nodes.push_back(child_node);
				queue.push_back(child_node);
					//cout << "Nodo "<< intTostring(child_node) << " inserito nella coda"<<endl;
			}
		}

	}


	return access_strings;
}


vector<string>  gi::dfa::get_access_strings_with_alphabet_symbols()
{

	// Access strings using the internal mapped symbol (usually integer)
	map<int, vector<SYMBOL>> access_strings = get_access_strings();

	// Final access strings usign the DFA alphabet
	vector<string> mapped_access_strings;

	////// Map string through the used alphabet
	string mapped_string = "";

	for(auto it=access_strings.begin(); it != access_strings.end(); ++it)
	{
		for(auto &it2 :  it->second )
			mapped_string = mapped_string + alphabet[it2] + " ";
		cout << "stringa mappata: "<<mapped_string << endl;

		mapped_access_strings.push_back(mapped_string);
		mapped_string = "";
	}

	return mapped_access_strings;
}






//vector<SYMBOL>  gi::dfa::equivalence_query_approximate(dfa* dfa_hp, string samplestestpath)
//{
//	#ifdef DEBUG_2
//	cout << endl << "--------------------------" << endl;
//	cout << "EQUIVALENCE QUERY" << endl;
//	cout << "--------------------------" << endl;
//	#endif
//
//
//	// Build union DFA of target dfa (thisone) and dfa_hp
//	dfa* dfa_union = this->unionDFA(dfa_hp);
//	#ifdef DEBUG_2
//	dfa_union->print_dfa_ttable("DFA UNIONE");
//	#endif
//
//
//	// Table-filling algorithm on union dfa
//	SYMBOL* distincts_table = dfa_union->table_filling();
//
//	// Extract list of equivalent states from table of distinct states
//	vector<int>* equivalent_states_list = dfa_union->equivalent_states_list_from_table(distincts_table);
//
//	// Verify if start states of dfas are equivalent
//	bool eq_dfa = false;
//	for(int i=0; i<equivalent_states_list[0].size(); ++i)
//		if(equivalent_states_list[0][i] == this->get_num_states())
//			eq_dfa = true;
//
//
//	// If start states aren't equivalent then dfas aren't equivalent. Generate a witness (counterexample)
//	//  (index of first state of dfa_hp in states list is dfa_hp->num_state)
//	vector<SYMBOL> witness;
//	if( !eq_dfa )
//		witness = dfa_hp->witness_approximate(dfa_hp,samplestestpath);
//
//	// Free allocated memory
//	if(distincts_table != NULL)
//		delete[] distincts_table;
//
//	return witness;
//}


//vector<SYMBOL> gi::dfa::witness_approximate(dfa* dfa_hp, string path_samples)
//{
//	// Test samples
//	int dim_positive;
//	int dim_negative;
//	vector<SYMBOL> *positive=NULL;
//	vector<SYMBOL> *negative=NULL;
//
//	vector<SYMBOL> witness;
//
//	bool negative_bool=true;
//
//	read_example_file(path_samples, positive, &dim_positive, negative,  &dim_negative);
//
//
//	for(int i=0; i<dim_positive; ++i)
//		if(dfa_hp->membership_query(positive[i]) == 0){
//			witness = positive[i];
//			negative_bool=false;
//			break;
//		}
//
//	if(negative_bool)
//		for(int i=0; i<dim_negative; ++i)
//			if(dfa_hp->membership_query(negative[i]) == 1){
//				witness = negative[i];
//				break;
//			}
//
//	if(positive) delete[] positive;
//	if(negative) delete[] negative;
//
//	return witness;
//}


bool gi::dfa::membership_query(vector<SYMBOL> str){

	// Check if arrive_state is ND (for DFA without sink state)
	int arrive_state = get_arrive_state(str);
	if(arrive_state == ND)
		return false;

	if(ttable[arrive_state][dim_alphabet] == 1)
		return true;
	else
		return false;
}


// Call it in a DFA union
vector<SYMBOL> gi::dfa::witness_from_table(SYMBOL* distinct, int start_state_dfa_hp)
{
	// Se automi NON equivalenti, creo la witness
	vector<SYMBOL> wit;

	int icoppia = 0;
	int jcoppia = start_state_dfa_hp;

	SYMBOL input;

	#ifdef DEBUG_2
	cout << "--- Creo il CONTROESEMPIO --- " << endl;
	#endif

	while(1)
	{
		//V cout << "Coppia: "<<icoppia << ","<<jcoppia<<"; vale: "<<distinct[icoppia][jcoppia]<<endl;
		int n = num_states;
		int k = (n*(n-1)/2) - (n-icoppia)*((n-icoppia)-1)/2 + jcoppia - icoppia - 1;

		if(distinct[k] == DFA_TF_STATE_N)
			cout << "PROBLEMA! Richiesta di controesempio con automi equivalenti";

		input = distinct[k];
		// Caso in cui lo start state dei 2, uno è accettante e l'altro no
		if(input == DFA_TF_STATE_X)
			break;

		//V cout << "stinga parziale: "<<wit<<" + " << input << endl;
		// TODO: il push piuttosto che concatenazione è corretto? old: wit = wit + input;
		wit.push_back(input);
		//V cout << "Vale: "<< input << ", ascii: "<< (int)input << ", diff_i:" << (int)input - 48 << endl;

		icoppia = ttable[icoppia][input];
		jcoppia = ttable[jcoppia][input];

		if(distinct[k] == DFA_TF_STATE_X){
			//V cout << "fine stringa" << endl;
			break;
		}
	}

	#ifdef DEBUG_2
	cout << "Controesempio è: "<< wit << endl;
	#endif

	return wit;
}


vector<int>* gi::dfa::equivalent_states_list_from_table(SYMBOL* distincts)
{
	#ifdef DEBUG_2
	cout << endl << "--------------------------" << endl;
	cout << "Lista stati equivalenti:" << endl;
	cout << "--------------------------" << endl;
	#endif

	vector<int>* stati_equivalenti = new vector<int>[num_states];
	int n= num_states;

	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j){
			int k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			if(distincts[k] == DFA_TF_STATE_N)
				stati_equivalenti[i].push_back(j);
		}

	#ifdef DEBUG_2
	cout << "N di stati finali: " << num_states << endl;
	for(int i=0; i<num_states; ++i)
		cout << "S["<<i<<"] --> "<< stati_equivalenti[i] << endl;

	cout << "--------------------------" << endl;
	#endif

	return stati_equivalenti;
}


double gi::dfa::get_complexity()
{
	return num_states;
}




// Generazione random di samples
//bool gi::dfa::generate_pos_neg_samples(int n_pos_samples,int n_neg_samples, const char * file_path)
//{
//	srand (time(NULL));
//	int rand_trans =0;
//	map<vector<SYMBOL>, int> samples;										// Stringa e indice dello stato
//	typedef	map<string,int>::const_iterator It;
//
//	int freq_to_stop = 1+2*(dim_alphabet);							// The inverse is probability to stop
//
//	string input_simbol = "";
//	bool finded_one_minimun_accepted_str = false;
//	string tmp_sample = "";
//	string sample_pos = "";
//
//	bool go_next = true;
//	int current_state = 0;
//
//	ofstream myfile;
//	myfile.open(file_path);
//
//	int num_attempts = 0;
//
//	// Write positve samples
//	while(samples.size() < n_pos_samples)
//	{
//		rand_trans = rand() % dim_alphabet;									// Input for first transition
//
//		input_simbol = intTostring(rand_trans);								// Int to String
//
//		finded_one_minimun_accepted_str = false;
//		tmp_sample = "";
//		sample_pos = "";
//
//		// Make positive sample
//		go_next = true;
//		current_state = 0;
//		int num_iteration = 0;
//		while(go_next)
//		{
//			//this->print_dfa_ttable("l");
//			current_state = ttable[current_state][rand_trans];
//			input_simbol = intTostring(rand_trans);
////			cout << "state:"<<current_state<< ", input:"<<input_simbol<<endl;
//
//			// Increment the final random string
//			tmp_sample = tmp_sample + " " +input_simbol;
//
//
//			// Visit DFA with input symbol
//			if(ttable[current_state][dim_alphabet] == 1){
//				finded_one_minimun_accepted_str = true;
//				sample_pos = tmp_sample;
//			}
//
//			// TODO Gestire stati pozzo
//
//			// Keep attention if stop because the upper bound about iteration... maybe some problem!
//			//if(num_iteration > 2000)
//				//cout << "! Limitated generation !"<<endl;
//
//			// Stop condition
//			if( ((rand() % freq_to_stop) == 1 && finded_one_minimun_accepted_str) || num_iteration > 5000)	// If I stop i should have some good string
//				go_next = false;
//
//			// Generate new input symbol
//			rand_trans = rand() % dim_alphabet;
//
//			++num_iteration;
//		} // END while(go_next)
//
//		samples[sample_pos] = 1;
//
//		++num_attempts;
//
//		// Num of attempts to do before stop.. could be a problem!
//		if(num_attempts > 60000)
//			break;
//
//
//		// If you want check that every string is inside language
//		/*if( membership_query(sample_pos) == 1)
//			cout << "Correct"<<endl;
//		else{
//			cout << "Wrong"<< endl;
//			exit(1);
//		}*/
//	} // END while(sample.size...)
//
//	if(num_attempts > 60000)
//		cout << "! Poor set of samples !" << endl;
//
//	// Generation of NEGATIVE SAMPLES
//	while(samples.size() - n_pos_samples < n_neg_samples)
//	{
//		//cout << "DimTot: "<<samples.size() <<", Positivi: "<<n_pos_samples <<", Negativi: "<<n_neg_samples<<endl;
//		int mod_type = rand() % 3;											// 0: substituting, 1:inserting, 2:deleting
//		int n_of_editing = getPoisson(3); //Poisson distribution centred on 3
//		#ifdef DEBUG_2
//		cout << "Poissoin ha ritornato "<< n_of_editing << " modifiche"<<endl;
//		#endif
//
//		int pos_to_edit = 0;
//		int n_changes = 0;
//
//		tmp_sample = "";
//
//		int* clean_sample = NULL;
//		for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		{
//			if((*p1).second == 1)
//				tmp_sample = (*p1).first;
//
//			// Extract sample without white space
//			int dim_clean_sample = tmp_sample.length()/2;
//			if(clean_sample != NULL)
//				delete[] clean_sample;
//			clean_sample = new int[dim_clean_sample];
//
//
//			int c = 0;
//			for(int i=1; i<tmp_sample.length(); i=i+2){
//				clean_sample[c] = tmp_sample[i] - '0';
//				c++;
//			}
//
//			if(dim_clean_sample == 0)
//				continue;
//
//			// Do changes
//			n_changes = 0;
//			while(n_changes < n_of_editing)									// Make "n_changes" for the single string up to "n_of_editing"
//			{
//				pos_to_edit = rand() % dim_clean_sample;
//				//cout << endl << endl <<"Posizione: "<<pos_to_edit<<endl;
//
//				// Round shift of "pos_to_edit" to right
//				int shift = dim_clean_sample-pos_to_edit;
//
//				/*cout << "Originale"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				rotate(clean_sample, clean_sample + pos_to_edit, clean_sample + dim_clean_sample);
//
//				/*cout << "Rotata"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				mod_type = rand() % 3;
//
//				++n_changes;
//
//				if(mod_type == 0){											// Substitution
//					//cout << endl<<"SOSTITUZIONE"<<endl;
//
//					int new_char = rand() % dim_alphabet;
//					while(clean_sample[0] == new_char)
//						new_char = rand() % dim_alphabet;
//
//					clean_sample[0] = new_char;
//
//					/*STAMPAcout << "Dopo sostituzione:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//				}
//				else if(mod_type == 1){
//					//cout <<endl<< "INSERIMENTO"<<endl;
//					int new_char = rand() % dim_alphabet;
//
//					int* ScleansampleTMP = new int[dim_clean_sample+1];
//					ScleansampleTMP[0] = new_char;
//					for(int i=1; i<dim_clean_sample+1; ++i)
//						ScleansampleTMP[i] = clean_sample[i-1];
//
//					delete[] clean_sample;
//					clean_sample = ScleansampleTMP;
//
//					dim_clean_sample += 1;
//
//					/*STAMPAcout << "Dopo inserimento:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//
//				}
//				else if(mod_type == 2){
//					//cout << endl << "CANCELLO"<<endl;
//					if(dim_clean_sample <= 1)
//						continue;
//
//					int* ScleansampleTMP = new int[dim_clean_sample-1];
//					for(int i=0; i<dim_clean_sample-1; ++i){
//						ScleansampleTMP[i] = clean_sample[i+1];
//					}
//
//					delete[] clean_sample;
//					clean_sample = ScleansampleTMP;
//
//					dim_clean_sample -= 1;
//
//					/*STAMPAcout << "Dopo eliminazione:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//				}
//
//				/*cout << "Rotata dopo cambiamento"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				//cout << "Ritorno originale dopo cambiamento:";
//				rotate(clean_sample, clean_sample + (dim_clean_sample-pos_to_edit), clean_sample + dim_clean_sample);
//				/*for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//			}
//
//			tmp_sample = "";
//			for(int i=0; i<dim_clean_sample; ++i)
//				tmp_sample = tmp_sample+" "+intTostring(clean_sample[i]);
//
//			if(tmp_sample.length() == 0)
//				continue;
//
//			if(this->membership_query(tmp_sample) == 0){
//				//cout << "stringa Non accettata correttamente:"<<tmp_sample<<";"<<endl;
//				samples[tmp_sample] = 0;
//			}
//
//			if(samples.size() - n_pos_samples >=  n_neg_samples)
//				break;
//		} //EDN FOR
//		if(clean_sample != NULL)
//			delete[] clean_sample;
//	}
//
//	myfile << dim_alphabet << "\n";
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		if((*p1).second  == 1)
//			myfile << "+ "+(*p1).first+"\n";
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		if((*p1).second  == 0)
//			myfile << "- "+(*p1).first+"\n";
//
//	myfile.close();
//
//
//	// Verify if the set of samples is structuraly complete
//	bool** scomplete = new bool*[num_state];
//	for(int i=0; i<num_state; ++i)
//		scomplete[i] = new bool[dim_alphabet+2]; 				// Column "dim_alphabet" true if state visited,
//																//   -    "dim_alphabet+1" true if accepting state used
//	for(int i=0; i<num_state; ++i)
//		for(int j=0; j< dim_alphabet+2; ++j)
//			scomplete[i][j] = false;
//	scomplete[0][dim_alphabet] = true;			// In current set often there aren't some transition to firs state
//
//	 // Init every state that is not ACC to true in the "dim_alphabet+1" column
//	for(int i=0; i<num_state; ++i)
//		if(ttable[i][dim_alphabet] != ACC){
//			scomplete[i][dim_alphabet+1] = true;
//		}
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1){
//		if((*p1).second  == 1){
//			current_state=0;
//			for(int i=0; i<(*p1).first.length(); ++i)
//			{
//				if((*p1).first[i] == ' ')
//					continue;
//
//				int input_simbol = (*p1).first[i]-'0';
//
//				// Explorated transition
//				scomplete[current_state][input_simbol] = true;
//
//				current_state = ttable[current_state][input_simbol];
//
//				// Accepting state used?
//				if(ttable[current_state][dim_alphabet] == ACC && (i+1 == (*p1).first.length()))
//					scomplete[current_state][dim_alphabet+1] = true;
//
//				// Explorated state
//				scomplete[current_state][dim_alphabet] = true;
//
//
//			}
//		}
//	}
//
//	for(int i=0; i<num_state; ++i)
//		for(int j=0; j< dim_alphabet+2; ++j)
//			if(!scomplete[i][j]){
//				//cout << "Non  completo per " <<i<<","<<j<<endl;
//				return false;
//			}
//
//
//	// Free memory
//	if(scomplete != NULL){
//			for(int i=0; i<num_state; ++i)
//				if(scomplete[i] != NULL)
//					delete[] scomplete[i];
//
//			delete [] scomplete;
//	}
//
//
//	return true;
//}

