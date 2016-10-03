/*
 * observationpack.cpp
 *
 *  Created on: 12 feb 2016
 *      Author: nicola
 */

#include <observationpack.h>
#include <limits>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <assert.h>
#include <queue>
#include <cmath>
#include <math.h>

#include "messages.h"
#include "utilities.h"

gi::observationpack::observationpack(const dfa &targetdfa)
{
	target =  new dfa(targetdfa);
	ptrCounterexampleFunction = &bindclass::rivest_shapire;
	
	vector<SYMBOL> empty_short_prefix;
	// set alphabet size
	dim_alfabeto = target->get_dim_alphabet();

	// allocate vector from the alphabet
	alfabeto = new vector<SYMBOL> [dim_alfabeto]; //MIO Array della dimenssione dell'alfabeto
	for(SYMBOL i=0; i<dim_alfabeto; ++i)
		alfabeto[i].push_back(i);
			
	//Inizialitation
	//gi::observationpack::component* C_e= NULL;
    component *C_e = new component; //call default constructor
    components[empty_short_prefix] =C_e;
    
    suffixes.push_back(empty_short_prefix);
    //initialitation of table anf of prefixes of Component with the empty string as short prefix
    C_e->mq[empty_short_prefix] = make_pair(target->membership_query(empty_short_prefix),1); 
    C_e->pref.insert(empty_short_prefix);
    for(SYMBOL i=0; i<dim_alfabeto; ++i)
    {
		C_e->mq[alfabeto[i]] = make_pair(target->membership_query(alfabeto[i]),1);
		C_e->pref.insert(alfabeto[i]);
	}
	
	nodes.push_back(make_pair(empty_short_prefix,nodeType::NODE_TYPE_LEAF)); //the root has index 0; leaf node that represent component with access sequence= empty string 	
	edges.push_back({NIL,NIL}); //The root hasn't children
	
	#ifdef DEBUG_2
	    cout<<endl<<"After the initialization you have"<<endl<<endl;
	    print_components();
	    print_discrimination_tree();
	#endif
	
	/*#ifdef DEBUG_1
		n_memb_query = dim_alfabeto+1; //+1 because  there is empty_short_prefix
		n_eq_query = 0;
	#endif*/
	#ifdef DEBUG_1
	    int nn=(targetdfa.get_num_states());
        experiment = new esperimenti(nn,dim_alfabeto,dim_alfabeto+1,0); //+1 because  there is empty_short_prefix
    #endif 
}
 
gi::observationpack::observationpack(const dfa &targetdfa , string counterexample_function) :observationpack(targetdfa)
{
	
	boost::algorithm::to_lower(counterexample_function);
	#ifdef DEBUG_1
        experiment->setType(counterexample_function);
    #endif 
	map<string,array<vector<SYMBOL>,3> (bindclass::*) (const bool , const vector<SYMBOL> &)> choicefunction = {{"rivest-shapire",&bindclass::rivest_shapire} , {"find-exponential",&bindclass::find_exponential} , {"partition-search",&bindclass::partition_search} , {"rivest-shapire-eager",&bindclass::rivest_shapire_eager}};
	auto iteratore = choicefunction.find(counterexample_function);

	if ( iteratore != choicefunction.end() )
		ptrCounterexampleFunction = iteratore->second;//quando alla fine del costruttore la map viene deallocata si elimina il puntatore a funzione e non quello a cui punta. Perciò non ci dovrebbe essere bisogno di scomodare la memoria dinamica
	else
	{
		cerr<<MSG_WRONG_FUNCTION;
		exit(EXIT_FAILURE);
	}
	
}

gi::observationpack::~observationpack(){
	if(alfabeto != NULL)
		delete [] alfabeto;

	if(target != NULL)
		delete target;
	
	if(experiment != NULL)
	    delete experiment;
	
	for(auto i=components.begin() ; i!=components.end() ; i++)
	    if(i->second != NULL)
	        delete i->second;	
}
	
inline void gi::observationpack::set_mq(map<vector<SYMBOL>,pair<bool,unsigned short int> > *mqCompPtr , const vector<SYMBOL>& key , bool value, unsigned short int n)
{
	
	pair<bool,unsigned short int> *p_key = &((*mqCompPtr)[key]);
	p_key->first = value;
	p_key->second += n;
}	

int gi::observationpack::rivest_shapire_classic(int low , int hight , const bool lambda_H_witness ,const vector<SYMBOL> &witness)
{
	//if(low<0 || hight>witness.size() || hight<=low)
    //{
		//cerr<<MSG_WRONG_INDEX<<endl; //or throw exception
		//exit(EXIT_FAILURE);
	//}
	
	bool alfa_mid;
	int mid;
	int size_wit=witness.size();
	
	while( (hight-low)>1)
	{
		mid = (low+hight)/2; //integer division
		vector<SYMBOL> u_current(witness.begin(),witness.begin()+mid);  //if reallocation isn't elegant:  u_current.resize(mid);  and copy(witness.begin(),witness.begin()+mid,u_current.begin())    This last solution is slightly(very little) more efficient 
		vector<SYMBOL> u_current_H_temp = get_ShortPrefixFromComponents(u_current);
		
		vector<SYMBOL> rest_witness(witness.begin()+mid,witness.begin()+size_wit);//You take remaind string of witness
		
		alfa_mid= target->membership_query(append_vectors(&u_current_H_temp, &rest_witness)) == lambda_H_witness;
		#ifdef DEBUG_1
		   experiment->set_n_memb_query();
		#endif
		   
		if(!alfa_mid)
		   low=mid;
		else
		   hight=mid;
		   		   
	}
	//Here low is the right index for the decomposition
	#ifdef DEBUG_3
	   assert(("Indice errato trovato nella funzione di decomposizione del controesempio",(low>=0 && low<size_wit)));//low must be between 0 (included) and size_wit (excluding)
	#endif
	return low;
	
}

/*
 * Questa funzione rispetto alla  versione di Rivest-Shapire presentata a pag. 87 di An Abstract Framework for Counterexample
 * Analysis in Active Automata Learning non torna l'indice per il quale alfa(i) != alfa(i+1) ma direttamente i prefissi discriminati .
 * Si è fatta questa scelta perchè si è trovato un modo di trovare la decomposizione e i prefissi discriminati 
 * (u_H_a e u_0_current_H) in fieri cioè mentre si esegue il while che deve trovare l'indice low per il quale
 * alfa(low) != alfa(low+1).  Invece la procedura descritta nel paper suddetto  trova l'indice low e poi usa quest ultimo per effetuare la decomposizione,
 * e poi ancora bisogna fare ulteriori chiamate di funzioni e operazioni per trovare u_H_a e u_0_current_H.  Quindi questa procedura
 * è un pò più efficiente di quella classica. Però questa funzione è anche usata da find-exponential e partition-search. Queste funzioni chiamano
 * rivest-shapire su indici low e hight diversi possibilmente da 0 e lunghezza della witness. Ciò non garantisce sempre il corretto funzionamento.
 * Quindi viene implementata anche la funzione rivest-shapire classica presentata nel paper che torna l'indice low per il quale
 * alfa(low) != alfa(low+1) che sarà usata come funzione ausiliaria da partition-search e find-exponential. Invece quando sarà usata 
 * la ricerca binaria di rivest-shapire sarà usata questa funzione che è leggermente più efficiente.  Il costo da pagare è piccolo. E' necessario
 * avere due versioni di Rivest-Shapire il che comporterà una crescita in memoria del programma, tuttavia le funzioni sono piccole quindi
 * non è un problema.
 */	
array<vector<SYMBOL>,3> gi::observationpack::rivest_shapire(const bool lambda_H_witness ,const vector<SYMBOL> &witness)
{
	//if(low<0 || hight>witness.size() || hight<=low)
    //{
		//cerr<<MSG_WRONG_INDEX<<endl; //or throw exception
		//exit(EXIT_FAILURE);
	//}
	
	int m=witness.size();
	int low=0,hight=m;
	bool alfa_mid;
	int mid;
	vector<SYMBOL> u_current_H;
	vector<SYMBOL> u_0_current_H;
	vector<SYMBOL> u_H_a; 
	#ifdef DEBUG_2
	   vector<SYMBOL> u;
	#endif
	
	while( (hight-low)>1)
	{
		mid = (low+hight)/2; //integer division
		vector<SYMBOL> u_current(witness.begin(),witness.begin()+mid);  //if reallocation isn't elegant:  u_current.resize(mid);  and copy(witness.begin(),witness.begin()+mid,u_current.begin())    This last solution is slightly(very little) more efficient 
		vector<SYMBOL> u_current_H_temp = get_ShortPrefixFromComponents(u_current);
		
		vector<SYMBOL> rest_witness(witness.begin()+mid,witness.begin()+m);//You take remaind string of witness
		
		alfa_mid= target->membership_query(append_vectors(&u_current_H_temp, &rest_witness)) == lambda_H_witness;
		#ifdef DEBUG_1
		   experiment->set_n_memb_query();
		#endif
		   
		if(!alfa_mid)
		{
		   low=mid;
		   u_current_H = u_current_H_temp; //operator= for vector make copy. Save the short prefix for low. if the decomposition is witness=uaw u_current_H will contain the short prefix of u
	       #ifdef DEBUG_2
	          u = u_current;
	       #endif
	    }
		else
		{
		   hight=mid;
		   u_0_current_H = u_current_H_temp; //if the decomposition is w=uaw   u_0_current_H at the end will contain the short prefix of ua 
	    }		   
	}
	//Here low is the right index for the decomposition
	vector<SYMBOL> a(witness.begin()+low,witness.begin()+low+1); //if the decomposition is w=uaw this is a
	u_H_a = append_vectors(&u_current_H,&a);
	
	#ifdef DEBUG_3
	   assert(("Indice errato trovato nella funzione di decomposizione del controesempio",(low>=0 && low<m)));//low must be between 0 (included) and size_wit (excluding)
	   assert(("Lo short prefix di u concatenato ad a non appartiene allo short prefix di ua. Ciò contraddice la teoria. Errore" , (get_shortPrefix(u_H_a) == u_0_current_H)));
	   assert(("La dimensione di a non e' 1",a.size()==1));
	#endif
	
	vector<SYMBOL> v(witness.begin()+low+1 , witness.begin()+m); //This instruction if the counterexample has got length==1 (isn't possible that is 0) can seems ambiguous.Rather should be works (v will be the empty vector, rapresenting the empty string.
	#ifdef DEBUG_2
	    cout<<"La decomposizione del controesempio w = "<<witness<<" in"<<endl
	        <<"w = uav e\'"<<endl<<"u = "<<u<<endl<<"a = "<<a<<endl<<"v = "<<v<<endl
	        <<"Il prefisso "<<u_H_a<<"appartiene al componente con short prefix "<<u_0_current_H<<endl<<endl;
	         
	#endif 
	
	return {u_H_a , u_0_current_H , v};
}

array<vector<SYMBOL>,3> gi::observationpack::find_exponential(const bool lambda_H_witness,const vector<SYMBOL> &witness)
{
	int m=witness.size();
	int index,low=0,hight=m;
	int ofs=1;
	bool found = false;
	bool alfa_index;
	
	while((index=(hight-ofs))>0 && !found)
	{
		vector<SYMBOL> pref_witness(witness.begin(),witness.begin()+index); //prendo i primi "index" symbols of witness
		vector<SYMBOL> rest_witness(witness.begin()+index, witness.begin()+m);//You take remaind string of witness 
		vector<SYMBOL> pref_witness_H = get_ShortPrefixFromComponents(pref_witness); //ottengo short prefix di pref_witness
		
		alfa_index= target->membership_query(append_vectors(&pref_witness_H, &rest_witness)) == lambda_H_witness;
		#ifdef DEBUG_1
		   experiment->set_n_memb_query();
		#endif
		if(!alfa_index)
		{
			low = hight - ofs;
			found = true;
		}
		else
		{
			hight -= ofs;
			ofs *=2;
		}
	}
	
	low = rivest_shapire_classic(low , hight , lambda_H_witness , witness);
	vector<SYMBOL> u(witness.begin(),witness.begin()+low);//if the decomposition is w=uaw this is u
	vector<SYMBOL> a(witness.begin()+low,witness.begin()+low+1); //if the decomposition is w=uaw this is a
	vector<SYMBOL> v(witness.begin()+low+1 , witness.begin()+m); //This instruction if the counterexample has got length==1 (isn't possible that is 0) can seems ambiguous.Rather should be works (v will be t
	
	vector<SYMBOL> u_H=get_ShortPrefixFromComponents(u);
    vector<SYMBOL> u_H_a=append_vectors(&u_H,&a);
    vector<SYMBOL> u_a_H=get_ShortPrefixFromComponents( append_vectors(&u , &a) );
    
	#ifdef DEBUG_3
	  assert(("La dimensione di a non e' 1",a.size()==1));
	  assert(("Lo short prefix di u concatenato ad a non appartiene allo short prefix di ua. Ciò contraddice la teoria. Errore" , (get_shortPrefix(u_H_a) == u_a_H)));
	#endif
	
	#ifdef DEBUG_2
	    cout<<"La decomposizione del controesempio w = "<<witness<<" in"<<endl
	        <<"w = uav e\'"<<endl<<"u = "<<u<<endl<<"a = "<<a<<endl<<"v = "<<v<<endl
	        <<"Il prefisso "<<u_H_a<<"appartiene al componente con short prefix "<<u_a_H<<endl<<endl;
	         
	#endif 
	
    return {u_H_a , u_a_H , v};
}

array<vector<SYMBOL>,3> gi::observationpack::partition_search(const bool lambda_H_witness,const vector<SYMBOL> &witness)
{
	int m=witness.size();
	int index,low=0,hight=m;
	bool found = false;
	bool alfa_index;
	double interval = m / log2(m) ;
	int step=static_cast<int>(floor(interval));
	
	while((index=(hight-step))>low && !found)
	{
		vector<SYMBOL> pref_witness(witness.begin(),witness.begin()+index); //prendo i primi "index" symbols of witness
		vector<SYMBOL> rest_witness(witness.begin()+index, witness.begin()+m);//You take remaind string of witness 
		vector<SYMBOL> pref_witness_H = get_ShortPrefixFromComponents(pref_witness); //ottengo short prefix di pref_witness
		
		alfa_index= target->membership_query(append_vectors(&pref_witness_H, &rest_witness)) == lambda_H_witness;
		#ifdef DEBUG_1
		   experiment->set_n_memb_query();
		#endif
		if(!alfa_index)
		{
			low = hight - step;
			found = true;
			break;
		}
		else
		    hight -= step;
	}
	
	low = rivest_shapire_classic(low , hight , lambda_H_witness , witness);
	vector<SYMBOL> u(witness.begin(),witness.begin()+low);//if the decomposition is w=uaw this is u
	vector<SYMBOL> a(witness.begin()+low,witness.begin()+low+1); //if the decomposition is w=uaw this is a
	vector<SYMBOL> v(witness.begin()+low+1 , witness.begin()+m); //This instruction if the counterexample has got length==1 (isn't possible that is 0) can seems ambiguous.Rather should be works (v will be t
	
	vector<SYMBOL> u_H=get_ShortPrefixFromComponents(u);
    vector<SYMBOL> u_H_a=append_vectors(&u_H,&a);
    vector<SYMBOL> u_a_H=get_ShortPrefixFromComponents( append_vectors(&u , &a) );
    
	#ifdef DEBUG_3
	  assert(("La dimensione di a non e' 1",a.size()==1));
	  assert(("Lo short prefix di u concatenato ad a non appartiene allo short prefix di ua. Ciò contraddice la teoria. Errore" , (get_shortPrefix(u_H_a) == u_a_H)));
	#endif
	
	#ifdef DEBUG_2
	    cout<<"La decomposizione del controesempio w = "<<witness<<" in"<<endl
	        <<"w = uav e\'"<<endl<<"u = "<<u<<endl<<"a = "<<a<<endl<<"v = "<<v<<endl
	        <<"Il prefisso "<<u_H_a<<"appartiene al componente con short prefix "<<u_a_H<<endl<<endl;
	         
	#endif 
	
    return {u_H_a , u_a_H , v};
}

/*
 * Qui non conviene gestire le cose come in rivest-shapire cioè calcolando i prefissi discriminati dal suffisso
 * mentre si calcola l'indice che individua il suffisso. Perchè per farlo nell'else dove modifico hight
 * si imposta l'indice a mid-1 (io ho calcolato già alfa(mid) e alfa(mid+1) ma non alfa(mid-1) )   quindi dovrei
 * fare dei calcoli aggiuntivi che potrebbe essere controproducente.
 * Quindi prima si calcola l'indice  e poi si fa la decomposizione e il calcolo dei prefissi.
 */
array<vector<SYMBOL>,3> gi::observationpack::rivest_shapire_eager(const bool lambda_H_witness, const vector<SYMBOL> &witness)
{
	int m=witness.size();
	int low=0 , hight=m-1;
	int mid;
	int beta_mid;
	int index=low;
	bool alfa_mid,alfa_mid_plus1;
	
	while( hight>low )
	{
		mid = (low+hight)/2;
		vector<SYMBOL> u_current(witness.begin(),witness.begin()+mid); 
		vector<SYMBOL> u_current_H_temp = get_ShortPrefixFromComponents(u_current);
		vector<SYMBOL> rest_witness(witness.begin()+mid,witness.begin()+m);//You take remaind string of witness
		alfa_mid = target->membership_query(append_vectors(&u_current_H_temp, &rest_witness)) == lambda_H_witness;
		
		vector<SYMBOL> u_current_plus1(witness.begin(),witness.begin()+mid+1);
		vector<SYMBOL> u_current_H_temp_plus1 = get_ShortPrefixFromComponents(u_current_plus1);
		vector<SYMBOL> rest_witness_plus1(witness.begin()+mid+1,witness.begin()+m);
		alfa_mid_plus1 = target->membership_query(append_vectors(&u_current_H_temp_plus1, &rest_witness_plus1)) == lambda_H_witness;
		
		#ifdef DEBUG_1
		   experiment->set_n_memb_query(2);
		#endif
		beta_mid = static_cast<int>(alfa_mid) + alfa_mid_plus1;
		if(beta_mid == 1)
		{
			index = mid;
		    break;
		}
		else if(beta_mid == 0)
		{
		    low = mid +1;
		    index = low;
		}
		else //beta_mid == 2 (entrambi gli alfa a 1)
		    hight=mid-1;//perchè alfa(mid)==1  quindi beta(mid-1)>=1 (perchè comprende alfa(mid) )
	}
	
	vector<SYMBOL> u(witness.begin(),witness.begin()+index);//if the decomposition is w=uaw this is u
	vector<SYMBOL> a(witness.begin()+index,witness.begin()+index+1); //if the decomposition is w=uaw this is a
	vector<SYMBOL> v(witness.begin()+index+1 , witness.begin()+m); //This instruction if the counterexample has got length==1 (isn't possible that is 0) can seems ambiguous.Rather should be works (v will be t
	
	vector<SYMBOL> u_H=get_ShortPrefixFromComponents(u);
    vector<SYMBOL> u_H_a=append_vectors(&u_H,&a);
    vector<SYMBOL> u_a_H=get_ShortPrefixFromComponents( append_vectors(&u , &a) );
    
	#ifdef DEBUG_3
	  assert(("La dimensione di a non e' 1",a.size()==1));
	  assert(("Lo short prefix di u concatenato ad a non appartiene allo short prefix di ua. Ciò contraddice la teoria. Errore" , (get_shortPrefix(u_H_a) == u_a_H)));
	#endif
	
	#ifdef DEBUG_2
	    cout<<"La decomposizione del controesempio w = "<<witness<<" in"<<endl
	        <<"w = uav e\'"<<endl<<"u = "<<u<<endl<<"a = "<<a<<endl<<"v = "<<v<<endl
	        <<"Il prefisso "<<u_H_a<<"appartiene al componente con short prefix "<<u_a_H<<endl<<endl;
	         
	#endif 
	
    return {u_H_a , u_a_H , v};
}	


/*#ifdef DEBUG_1
    int gi::observationpack::get_n_memb_query()
    {
	    return n_memb_query;
    }

    int gi::observationpack::get_n_eq_query()
    {
	    return n_eq_query;
    }
#endif */


void gi::observationpack::print_components()
{
	cout <<       "---------------------------------------------";
	cout << endl<<"                 COMPONENTS " << endl;
	cout <<       "---------------------------------------------" << endl;
	
	int numComponent = 1;
	set<vector<SYMBOL> >  *tmpPref;
	map<vector<SYMBOL>, pair<bool,unsigned short int> >   *tmpMq;
	component *tmpComponent;
	
	vector<SYMBOL> tmpVect;
	
	for(auto i=components.begin() ; i!=components.end() ; i++)
	{
		cout<<"COMPONENT "<<numComponent++<<endl;
		cout<<"Short Prefix: "<<i->first<<endl;
		
		cout << "Suf: ";
		for(auto s=suffixes.begin() ; s!=suffixes.end() ; s++) //for ONEGlobally this execution give the same result ad any iteration. For clarity is this
		    cout<<*s <<", " ;
		cout<<endl;
		
		tmpComponent = i->second;
		tmpPref = &(tmpComponent->pref);
		tmpMq = &(tmpComponent->mq);
		
		for(auto p = tmpPref->begin() ; p!=tmpPref->end() ; p++)
		{
			tmpVect = *p;
			cout << tmpVect<<" | ";
			for(auto s2 = suffixes.begin() ; s2 != suffixes.end() ; s2++)
			{
			    cout<< ((*tmpMq)[append_vectors(&tmpVect, &(*s2))]).first ; 
			    #ifdef DEBUG_4
			        cout<<"\\"<<((*tmpMq)[append_vectors(&tmpVect, &(*s2))]).second<<" "; //stampa del contatore
			    #endif
			}
			cout<<endl;    
		}
		cout<<endl<<endl;    
	}
}

void gi::observationpack::print_discrimination_tree()
{
	cout <<left<< "---------------------------------------------";
	cout << endl<<"            DISCRIMINATION TREE " << endl;
	cout <<       "---------------------------------------------" << endl;
	
	//print edges
	cout<<setw(10)<<"NODE"<<setw(12)<<"Left_Child"<<"Right_Child"<<endl;
	for(int i = 0 ; i!=edges.size() ; i++)
	{
		cout<<setw(10)<<i;
		//così non si può fare
		//cout<<setw(12)<<(((numLeafNode=edges[i][0])==NIL) ? "NIL" : numLeafNode);
	    //cout<<edges[i][1]<<(((numLeafNode=edges[i][1])==NIL) ? "NIL" : numLeafNode);
	    
	    if(edges[i][0] == NIL)
	        cout<<setw(12)<<"NIL";
	    else
	        cout<<setw(12)<<edges[i][0];
	        
	    if(edges[i][1] == NIL)
	        cout<<setw(12)<<"NIL"<<endl;
	    else
	        cout<<setw(12)<<edges[i][1]<<endl;    
	    
	}
	
	//print for every node the label
	int numNode=0;
	cout<<endl<<endl<<setw(10)<<"NODE"<<setw(7)<<"LEAF"<<"LABEL"<<endl;
	for(auto i=nodes.begin() ; i!=nodes.end() ; i++)
		cout<<setw(10)<<numNode++<<setw(7)<<(((*i).second==(nodeType::NODE_TYPE_LEAF))?"YES":"NO")<<(*i).first<<endl;

}

//Questa funzione modifica il membro target dell'oggetto chiamante. Per il resto non è distruttiva
gi::dfa*  gi::observationpack::run_observationpack(bool approximate, string samplestestpath)
{
	dfa * dfahp = NULL;
	// Minimize tgdfa
	dfa * t1 = target->minimize_TF();

	if(target != NULL)
		delete target;

	target = t1; //MIO In pratica per l'oggetto corrente(il chiamante) modifico per il membro target ciò a cui punta(non viene chiamato nessun costruttore). Cioè 
                 //MIO il membro target non punterà più al DFA creato dal file passato(lstar.txt ad esempio) ma punterà al DFA minimizzato. Qui
//Quello sopra è l'unico membro modificato dell'oggetto chiamante
//Da adesso in poi si decide di fare una copia dell'oggetto in modo da non modificare niente nell'oggetto chiamante
//In questo modo si occulta al chiamante di vedere le componenti e il discrimination tree finale. A quest ultimo si da solo il DFA inferito.
//In questo modo il chiamante non dedurrà nulla di come funziona il processo di costruzione dell'ipotesi.
    
    observationpack* obpack = new observationpack(*target);
    #ifdef DEBUG_1
       obpack->experiment->setNumState(this->experiment->getNumState());
       obpack->experiment->setType(this->experiment->getType());
    #endif
    obpack->ptrCounterexampleFunction =  ptrCounterexampleFunction;  
    
    int count_generation = 0;
    bool lambda_H_witness;
    bool lambda_target_witness;
    bool eq;
    vector<SYMBOL> witness;
//Quando viene tornato il controesempio non so se è accettante nel target e rigettante nell'ipotesi oppure viceversa
//Siccome quest'informazione è necessaria( mi serve lambda_H_witness) devo fare la membership_query del controesempio nell'ipotesi

array<vector<SYMBOL>,3> decomposition={{{NIL},{NIL},{NIL}}};

    while(1)
    {
        #ifdef DEBUG_2
		    cout << endl << endl <<"********************************" << endl;
		    cout <<  "         GENERATION "<<count_generation << endl;
		    cout << "********************************" << endl<<endl;
		#endif
		
		obpack->closePack(decomposition);
        #ifdef DEBUG_2
            cout<<"Later closdness of components you have:"<<endl<<endl;
		    obpack->print_components();
            obpack->print_discrimination_tree();
        #endif
    
        if(dfahp != NULL)
            delete dfahp;
    
        dfahp = obpack->obk_to_dfa();
        #ifdef DEBUG_2
            cout<<endl<<endl<<"The current hypothesis DFA is:"<<endl;
		    dfahp->print_dfa_ttable("DFA_HP");
	    #endif 
	
	     
    
        if(!approximate)
        {
			if(count_generation != 0)
			{
				lambda_H_witness = dfahp->membership_query(witness);
				#ifdef DEBUG_1
		            obpack->experiment->set_n_memb_query(); // Though the membership query above is asked in the automaton of the hypothesis(no in the target) it should still count
		        #endif
		    }
	        if(count_generation == 0 || lambda_H_witness == lambda_target_witness)
	        {
				witness = target->equivalence_query(dfahp);
			    #ifdef DEBUG_1
                    obpack->experiment->set_n_eq_query();
                #endif 
			    eq = true;
			} 
			else 
			    eq = false;
		}
		else //is approximate
		 cout<<"!!!!!!!!!!!!!HAI FATTO UNA CHIAMATA APPROSSIMATA!!!!!!!!!!!!!!!"<<endl;
			//witness = target->equivalence_query_approximate(dfatmp,samplestestpath);
	    
	        
	        
      #ifdef DEBUG_4
      if(!approximate){
      string nome = "W"+intTostring(count_generation);
      string percorso="IpotesiIntermedie/";
      dfahp->print_dfa_dot_mapped_alphabet("OBSERVATIONPACK", (percorso  +nome+".dot" ).c_str());
	  }
	  #endif
    
        if(0 != witness.size())//In Observation Pack we can assume the length counterexample greater than 0. See 3.1.2, page 85 of An abstract framework for counterexample analysis in active automata learning (Isberner Steffen)
        {
			if(!approximate)
			{
				if( eq )
				{
					lambda_H_witness = dfahp->membership_query(witness); //MQ of new witness in the hypothesis
					lambda_target_witness = ( !lambda_H_witness );
					#ifdef DEBUG_1
					    obpack->experiment->setSizeCounterexamples(witness.size());
			            obpack->experiment->set_n_memb_query(); // Though the membership query above is asked in the automaton of the hypothesis(no in the target) it should still count
			        #endif
				}
				//else witness isn't changed. lambda_H_witness is setted above. lambda_target_witness is the same of previous iteration because
			    //the witness is the same and the target also
				#ifdef DEBUG_1
				else
				{
					    obpack->experiment->setSizeCounterexamples(witness.size()); //Viene riportata la dimensione, che sarà la stessa del controesempio precedente
				}
				#endif
				
		    }
			
			decomposition=((obpack->*ptrCounterexampleFunction)(lambda_H_witness,witness));
	        obpack->update_from_counterexample(decomposition , lambda_H_witness);//a[2] is the suffix found. membership query di a[0]+a[2] e' diverso da lambda_H_witness.Invece MQ(a[1]+a[2]) e' uguale a lambda_H_witness
	        //Passando anche lambda_H_witness posso risparmiare 2 MQ. Al massimo questa funzione sarà chiamata n(numero stati del target) volte. QUindi posso risparmiare fino a 2n MQ.
	        //Le 2 MQ già fatte saranno quelle presenti nel componente che ha come short prefix a[1]. Quindi gli if vanno inseriti solo in questo componente
	        #ifdef DEBUG_2
	            cout<<endl<<"Print the components after adding the suffix "<<decomposition[2]<<endl;
	            obpack->print_components();
	        #endif
	        
	        #ifdef DEBUG_1
	            obpack->experiment->setSizeSuffixes(decomposition[2].size());
	        #endif
	    }
	    else
	    {
	        cout << "OBSERVATION PACK: Automaton found!"<<endl;
		    //dfahp = dfahp->minimize_TF(); IL DFA inferito con obpack dovrebbe essere il dfa minimo.quindi non è necessario minimizzarlo  
            break;
	    }
	    ++count_generation; 
	    #ifdef DEBUG_3
	        assert(("The number of iterations overcome the numbers of states of the target",count_generation<=target->get_num_state())); //Ad ogni passo produco almeno un nuovo componente(stato). Il DFA è il minimo quindi in al massimo n (numeri stati del target)n passi devo terminare l'inferenza (upper bound)
	    #endif 
    }
    
    #ifdef DEBUG_1
    obpack->experiment->myPrint("../esperimenti/experiments.txt");
    #endif
    
    delete obpack;

    //#ifdef DEBUG_1
      //  dfa dfaInferred(*dfahp , *experiment);
    //#else    
        
    //#endif
    
    dfa dfaInferred(*dfahp); //MIO Chiamo il costruttore di copia per creare l'oggetto dfaInferred a partire dall'oggetto dfahp
    
	//if(dfahp != NULL)
		//delete dfahp;
  
	//return dfaInferred; //Ritorno dfaInferred e non dfahp perchè quest ultimo è memoria dinamica. Il client deve avere la sua copia privata
    return dfahp;
}

vector<SYMBOL> gi::observationpack::get_shortPrefix(const vector<SYMBOL> &prefix)
{
	vector<SYMBOL> notPresent={NIL};
    	
	for(auto i=components.begin() ; i!=components.end() ; i++)
	{
		if(i->second->pref.find(prefix) != i->second->pref.end() )
		   return i->first;  //copy constructor (of vector) called when there is  a return
	}
	
	return notPresent;
}

void gi::observationpack::update_from_counterexample(const array<vector<SYMBOL>,3> &decomposition , bool lambda_H_witness)
{
	set<vector<SYMBOL> >  *prefPtr;
	map<vector<SYMBOL>, pair<bool , unsigned short int> >   *mqPtr;
	component *componentPtr;
	vector<SYMBOL> tmpVect;
	//vector<SYMBOL> tmpVect2=decomposition[2];
	//vector<SYMBOL> tmpVect3;
	
	suffixes.push_back(decomposition[2]); //add suffix to suffix set
	for(auto i = components.begin() ; i != components.end() ; i++)
	{
		componentPtr = i->second;
	    prefPtr = &(componentPtr->pref);
		mqPtr = &(componentPtr->mq);
		if( get_shortPrefix(decomposition[1]) == i->first )
		{
			vector<SYMBOL> vect1 = append_vectors(&(decomposition[0]) , &(decomposition[2]));
			vector<SYMBOL> vect2 = append_vectors(&(decomposition[1]) , &(decomposition[2]));
			set_mq(mqPtr , vect1 , (!lambda_H_witness) , 1);
			set_mq(mqPtr , vect2 , lambda_H_witness , 1);
			
			for(auto pr = prefPtr->begin() ; pr!=prefPtr->end() ; pr++)
			{
				//tmpVect3 = (*pr);//this is necessary because the deferencing of a set iterator return a const vector and after will be mismatching with the function append_vectors
				//Rather I changed the prototype of append_vectors in utilities.h for avoid expensive copies here
				tmpVect = append_vectors(&(*pr) , &(decomposition[2])); //concateno il prefisso corrente con IL suffisso
				if( (*mqPtr).count(tmpVect) == 0 ) //Questa è necessaria perchè  per decomposition[0] e decomposition[1] non devo fare la MQ 
				{    //Se si usasse la soluzione 2) prospettata nell'e-mail il codice nell'else sotto sarebbe lo stesso e il codice potrebbe essere unificato
					set_mq(mqPtr , tmpVect , target->membership_query(tmpVect) , 1); // Here you are sure that is the first entry for tmpVect(because before there is the count)
					#ifdef DEBUG_1
					    experiment->set_n_memb_query();
					#endif
			    }
			    else
			    {
					if((tmpVect != vect1) && (tmpVect!= vect2))//perchè per vect1 e vect2 ho già settato sopra. Devo evitare di settare due volte
			            set_mq(mqPtr , tmpVect , ((*mqPtr)[tmpVect]).first , 1);//Se è presente devo comunque incrementare di uno
			    }
			}
			
		}
		else
		{
			for(auto pr = prefPtr->begin() ; pr!=prefPtr->end() ; pr++)
			{
				//tmpVect3 = (*pr);//this is necessary because the deferencing of a set iterator return a const vector and after will be mismatching with the function append_vectors
				//Rather I changed the prototype of append_vectors in utilities.h for avoid expensive copies here
				tmpVect = append_vectors(&(*pr) , &(decomposition[2]));
				//if( (*mqPtr).count(tmpVect) == 0 )   Questa è la soluzione 2) prospettata nell'e-mail. In pratica prima si controlla se l'elemento è già presente prima di fare la MQ 
				set_mq(mqPtr , tmpVect , target->membership_query(tmpVect) , 1);
				#ifdef DEBUG_1
				    experiment->set_n_memb_query();
				#endif
			}
			
	    }	
	}
}

gi::dfa* gi::observationpack::obk_to_dfa()
{
	map<vector<SYMBOL>, array<int,2> > states; //create the corrispondence between short prefix of a component and number of state (store also if the state that represents the component is accepting or not)
	int count_state = components.size();
	int index=1,indexEmptySP;
	vector<SYMBOL> empty_short_prefix;
	vector<SYMBOL> current_shortPr,stato_arrivo;
	
	dfa* dfaOBK = new dfa(count_state, target->get_dim_alphabet(), target->get_alphabet(), 0);
	int** dfaOBKtable = dfaOBK->get_ttable(); //MIO In questo modo modifico direttamente la tabella di transizione dell'oggetto dfaOBT modificando dfaOBKtable
	
	//assign index to components (unique constraint is that component with short prefix= empty string have to be index 0)
	for(auto i=components.begin() ; i!= components.end() ; i++,index++)
	{
		current_shortPr=i->first;
		if((current_shortPr) == empty_short_prefix)
		{
			indexEmptySP=index;
			states[current_shortPr][0] = 0;
			--index;
		}
		else
			states[current_shortPr][0] = index;
			
		states[current_shortPr][1] = ((i->second->mq)[current_shortPr]).first; //etablish if the state that rapresent the component is accepting or not
	} 
	
	index=1;
	int counter=1;
	//La struttura dati components non ha subito modifiche quindi l'ordine di restituzione degli elementi è lo stesso del ciclo sopra(cioè sono sicuro che il primo elemento è quello con indice 1 oppure 0 se è l'empty short prefix) 
	for(auto C=components.begin() ; C!= components.end() ; C++,index++,counter++)
	{
		current_shortPr=C->first;
		if(counter == indexEmptySP)
		{
			indexEmptySP=-1; //Di modo che non ci sia eguaglianza anche al ciclo successivo quello in cui risultano uguali counter e indexEmptySP
			index=0;
			--counter;
		}
		
		
		for(SYMBOL i=0; i<dim_alfabeto+1; ++i)
		{
			if(i==dim_alfabeto)
		        dfaOBKtable[index][dim_alfabeto] = states[current_shortPr][1];
		    else
		    {    
				stato_arrivo = append_vectors(&current_shortPr, alfabeto + i);
		        dfaOBKtable[index][i]  = states[get_shortPrefix(stato_arrivo)][0]; //Per costruzione stato_arrivo deve sempre appartenere a qualche componente. Potrei introdurre una assert. 
		    } 
		    
		}
		
		index=counter;
		
	}
	return dfaOBK;
}

void gi::observationpack::closePack(const array<vector<SYMBOL>,3> &decomposition)
{
	const vector<SYMBOL> nullo ={NIL};
	bool firstTime = true; //flag to indicate that is the first iteration
	bool firstTimeAbsolutely= (decomposition[2] == nullo) ; //flag to indicate that is the first iteration of closePack in  absolute
	bool isClosed = false;
	bool exitFor;
	queue<vector<SYMBOL> > W;
	vector<SYMBOL> empty_short_prefix,shortPr,rowShortPr , rowPref,pr,suf,currentComp,arriveComp;
	set<vector<SYMBOL> >  *prefC;
	component *cPtr;
	pair<vector<SYMBOL>::iterator,vector<SYMBOL>::iterator> result;
	int indexSuf;
	
	if(firstTimeAbsolutely)
	    W.push(empty_short_prefix); //la prima e unica componente
	//else  W is empty
	    
	
	while(1)
	{
		if((!firstTime) || firstTimeAbsolutely) //if isn't the first time or is the first time in absolute(the first time in absolute alone the parameters will be NULL
		{
			//Ora trovo  shortPr pr e suf  (u_0 u e v)
			exitFor=false;
			for( auto  c = components.begin() ; c!=components.end() ; c++)
			{
				shortPr = c->first;
				cPtr    = c->second; 
				prefC= &(cPtr->pref);
				rowShortPr = get_row(shortPr , &(cPtr->mq));
				for(auto p = prefC->begin() ; p!= prefC->end() ; p++)
				{
					//Si potrebbe fare(perchè tra i prefissi c'è anche lo short prefix) if(*p == shortPr) continue;   Ma forse non conviene. (Funziona in entrambi i casi è solo per efficienza)
					rowPref = get_row(*p , &(cPtr->mq));
					result = mismatch(rowShortPr.begin(),rowShortPr.end(),rowPref.begin()) ; //mismatch halt it when find the iterator to first element that differ (if there isn't this element return the end() iterator)
					if(result.first != rowShortPr.end() ) // i due vettori non hanno tutti gli elementi uguali
					{
						pr = *p;
					    indexSuf = result.first -  rowShortPr.begin();
					    suf = suffixes[indexSuf];
					    isClosed = false;
					    exitFor=true;
					    break; //IN REALTA' QUI LA COSA MIGLIORE ERA USARE UNA GOTO. UNP DEGLI ESEMPI CONSIGLIATI PER LA GOTO E' PROPRIO PER USCIRE DAI NESTED LOOP.  Senò si può inserire questo codice in una funzione e usare return isClosed.  Altrimenti come quà usare un flag (che comunque mi serve anche dopo)
					}
				}
				if(exitFor)
				    break;
			}
			if(!exitFor) //cioè non ho mai trovato due righe diverse e quindi sono arrivato qui senza usare dei break
			    isClosed=true;
		}
		//       Se sono qui è nell'else la prima volta ma non in assoluto che viene chiamato ClosePack. Qui sono anche sicuro che le component
		//       non sono chiuse e vanno splittate (short prefix, prefisso e suffisso ce li ho nel parametro decomposition)
		else
		{
			pr      = decomposition[0];
			shortPr = decomposition[1];
			suf     = decomposition[2];
		}
		
		if(isClosed)
		    if(!firstTimeAbsolutely) //Perchè la prima volta in assoluto anche se le componenti fossero chiuse devo comunque eseguire il codice a seguire. L'autore della tesi non ha bisogno di questo controllo perchè la prima volta in assoluto le componenti non sono complete e quindi è sicuro dato che manca la completezza di eseguire comunque il while e l'inner while. Io invece completo le componenti altrove. Quando chiamo ClosePack le componenti sono già complete e quindi se le componenti fossero anche chiuse  non eseguirei il ciclo while e l'inner while commettendo un errore.
		        break;
		    else
		        isClosed=false;
		else
		{
			split(shortPr,pr,suf);
			#ifdef DEBUG_4 //per una stampa ancora più approfondita
			cout<<endl<<"Componente con access sequence = "<<shortPr<<" splittato."<<endl;
			cout<<"Il nuovo componente ha access sequence = "<<pr<<endl<<"Il suffisso e\' = "<<suf<<endl;
			
			    print_components();
			    print_discrimination_tree();
			#endif
			W.push(pr);
		}
		
		firstTime = firstTimeAbsolutely=false;
		while(!W.empty())
		{
			currentComp = W.front();
			W.pop();
			for (SYMBOL i=0 ; i<dim_alfabeto ; i++)
			{
				arriveComp = append_vectors(&currentComp, alfabeto + i); 
				if( sift(arriveComp) )
				    W.push(arriveComp);
				    
				#ifdef DEBUG_4 //per una stampa ancora più approfondita
				    cout<<endl<<"Il sift di "<<arriveComp<<" ha prodotto le nuove seguenti strutture dati"<<endl;
			        print_components();
			        print_discrimination_tree();
			    #endif    
			}
		}
	}//close while(1)
	
}

vector<SYMBOL> gi::observationpack::get_row(vector<SYMBOL> prefix , map<vector<SYMBOL>, pair<bool,unsigned short int> > 	*mqComp){
	vector<SYMBOL> row;

	for(auto last_index=suffixes.begin(); last_index!=suffixes.end(); ++last_index){
		vector<SYMBOL> tmpVect = prefix;
		tmpVect.insert(tmpVect.end(), last_index->begin(), last_index->end()); //si poteva chiamare append. Così evito il costo di un'ulteriore chiamata di funzione

		row.push_back((((*mqComp)[tmpVect]).first)?DFA_STATE_ACCEPTING:DFA_STATE_REJECTING);
	}

	return row;
}

void gi::observationpack::split(const vector<SYMBOL> &shortPr ,const vector<SYMBOL> &pr ,const vector<SYMBOL> &suf)
{
	component *C = components[shortPr];
    component *newComponent = new component; //call default constructor
    components[pr] =newComponent; //Pr become the access sequence of new component (Pr is in the component with access sequence ShortPr and from there must be eliminated	
    
    set<vector<SYMBOL> >  *prefnewC  = &(newComponent->pref);
    map<vector<SYMBOL>, pair<bool,unsigned short int> >   *mqnewC = &(newComponent->mq);
    set<vector<SYMBOL> >  *prefC  = &(C->pref);
    map<vector<SYMBOL>, pair<bool,unsigned short int> >   *mqC = &(C->mq);
    pair <bool,unsigned short int> *pairC;
    vector<SYMBOL> tempVect;
   
    //Operations for split the component with short prefix=shortPr.  This component isn't closed
    bool mq_shortPr_suf = ((*mqC)[append_vectors(&shortPr,&suf)]).first;  //Mq stored in mq of concatenation of shortPr and Suf. 
    //Si dovrebbe controllare che mq_shortPr_suf sia diverso da ((*mqC)[append_vectors(pr,&suf)]).first  perchè la funzione richiede parametri siffatti in input
    for(auto i=prefC->begin() ; i!=prefC->end() ; )
    {
		if( mq_shortPr_suf !=  (((*mqC)[append_vectors(&(*i),&suf)]).first) )
		{
			(*prefnewC).insert(*i);
			for(auto s = suffixes.begin() ; s!=suffixes.end() ; s++)
			{
				tempVect = append_vectors(&(*i),&(*s));
				pairC =  &( (*mqC)[tempVect] );
				set_mq(mqnewC , tempVect , pairC->first ,1);
				--(pairC->second) ;
				if(pairC->second == 0)  //era l'unica entry, allora puoi eliminare
				    (*mqC).erase(tempVect); 	   
			}
			i = (*prefC).erase(i); //erase return iterator following the last removed element
		}
		else
		++i;
	}
	
	//Now the code for modify the discrimination tree
	
	//the node with short prefix = shortPr must be a leaf
	//this search is in a vector and has a linear cost but is performed just once for each call of the function.
	//The cost here not to use a map is amply repaid in the rest of this function and in sift function
	auto result = find(nodes.begin(),nodes.end(), pair<vector<SYMBOL>,nodeType> (shortPr,nodeType::NODE_TYPE_LEAF));
	#ifdef DEBUG_3
	    assert(("The leaf in the discrimination tree with this short prefix was not found",result!=nodes.end()));
	#endif
	int index = result - nodes.begin(); //index of nodes where there is the pair with access sequence=shortPr 
	int nodeNumber = nodes.size(); //Read the number of nodes in the discrimination tree
	bool shortPrLeftChild;
	nodes[index].first = suf; //set the label of node that before was with label=shortPr at suf (the discriminator)
	nodes[index].second= nodeType::NODE_TYPE_INNER; //the node described above was a leaf. Now become a inner node
	
	nodes.push_back(make_pair(shortPr,nodeType::NODE_TYPE_LEAF)); //insert in nodes the node leaf with access sequence=shortPr
	nodes.push_back(make_pair(pr,nodeType::NODE_TYPE_LEAF));  //insert in nodes the node leaf with access sequence=pr
	
	edges.push_back({NIL,NIL});//set the edges for node shortPt at NIL (It has not children). And you increment of 1 the dimension of row of edges
	edges.push_back({NIL,NIL});//set the edges for node pr at NIL (It has not children)
    //if(mq_shortPr_suf) 
		//shortPrLeftChild=false;
	//else
		//shortPrLeftChild=true;
	shortPrLeftChild= ( mq_shortPr_suf ? false : true ) ; //the node with short prefix=shortPr must be right child of inner node with shortPr=suf (else figlio sinistro)	
	edges[index][0] = nodeNumber + (!shortPrLeftChild); //a bool is automatically converted to 0 or 1 when you sum with a int
	edges[index][1] = nodeNumber + shortPrLeftChild;//a seconda del risultato della MQ il figlio destro e sinistro saranno invertiti (nodeNumber e nodeNumber+1 sono gli indici dei due nuovi nodi foglia. Bisogna solo stabilire chi va a destra e chi va a sinistra)
    //index of node with short prefix =shortPr è nodeNumber
}

bool gi::observationpack::sift(const vector<SYMBOL> &pr)
{
	unsigned short int indexNode=0; //0 is the index of the root
	unsigned short int child;
	map<vector<SYMBOL>, pair<bool,unsigned short int> >::iterator search;
	map<vector<SYMBOL>, pair<bool,unsigned short int> >   mqTemp; //Qui salvo le MQ effettuate durante il sift (la navigazione) del discrimination tree (in modo da non doverle rifare due volte)
	vector<SYMBOL> labelNode;
	vector<SYMBOL> tempVect;
	bool resultMq;
	
	while((nodes[indexNode].second) == nodeType::NODE_TYPE_INNER) //namely isn't a leaf. Like this ss the same thing
	{
		labelNode = nodes[indexNode].first;
		tempVect=append_vectors(&pr, &labelNode);
		resultMq = target->membership_query(tempVect);
		set_mq(&mqTemp , tempVect , resultMq , 1);  //Store the MQ
		#ifdef DEBUG_1
		   experiment->set_n_memb_query();
		#endif
		child = edges[indexNode][resultMq]; //resulMq is bool. Automatically converted to 0 or 1
	    if(child != NIL)
			indexNode = child;
		else
		{
			//Modify the discrimination tree
			nodes.push_back(make_pair(pr,nodeType::NODE_TYPE_LEAF)); //add new node leaf with access sequence pr
			nodes[indexNode].second= nodeType::NODE_TYPE_INNER; //indexNode potrebbe anche già essere un nodo interno come potrebbe essere una foglia. In entrambi i casi deve  diventare(o restare) un nodo interno
			edges[indexNode][resultMq] = nodes.size()-1; //l'arco del nodo indexNode che prima puntava a NIL adesso deve puntare al nuovo nodo aggiunto
		    edges.push_back({NIL,NIL});//add space and set children for new node added
		
		    //create new component with access sequence pr
		    component *newC = new component; //call default constructor
		    components[pr] =newC; //Pr become the access sequence of new component
		    (newC->pref).insert(pr); //add access sequence pr as prefix
		    //Now I must complete the table(mq) of newC (fill the holes)
		    //The following code  is repeated two times (see below) but non credo sia conveniente creare una funzione (la sift viene chiamata spesso. Al massimo se ne potrebbe fare una inline)(Inoltre il codice sotto può differe da quello immediatamente seguente se si sceglie di cercare anche nella tabella locale prima di fare la MQ (nel codice immediatamente seguente invece la tabella locale è vuota)
		    for(auto s = suffixes.begin() ; s!= suffixes.end() ; s++)
		    {
				tempVect = append_vectors(&pr, &(*s));
				search = mqTemp.find(tempVect);
				if(search != mqTemp.end() ) //the element is present
				    set_mq(&(newC->mq) , tempVect , (search->second).first , (search->second).second);
				else
				{
					set_mq(&(newC->mq) , tempVect , target->membership_query(tempVect) , 1);
				    #ifdef DEBUG_1
		                experiment->set_n_memb_query();
		            #endif
				}        
			}
		    return true;
		}
	}
	
	//Here indexNode is a leaf
	labelNode = nodes[indexNode].first;//labelNode is a short prefix
	component *C = components[labelNode];
	set<vector<SYMBOL> > *prPtr = &(C->pref);
	
	//now you must add pr to C if pr doesn't belong already to C's prefixes. Moreover If exist is already complete OK
	if(prPtr->find(pr) == prPtr->end() )
	{
		 prPtr->insert(pr);
		 //Here I complete the table (fill the holes).  For this prefix and for each suffix I must do a MQ (tranne quelle già salvate in mqTemp) 
	     for(auto s = suffixes.begin() ; s!= suffixes.end() ; s++)
	     {
			 tempVect = append_vectors(&pr, &(*s));
			 search = mqTemp.find(tempVect);
			 if(search != mqTemp.end() ) //the element is present
				 set_mq(&(C->mq) , tempVect , (search->second).first , (search->second).second);
			 else //the element isn't present. IN THIS CASE FIRST I COULD DO A SEARCH IN THE TABEL (mq) OF THE COMPONENT. IF If there is not even there I will do MQ. HERE I NO USE THIS STRATEGY (IS DIFFICULT FOUNT A ELEMENT IN THE LOCAL TABLE)
			 {	 
				 set_mq(&(C->mq) , tempVect , target->membership_query(tempVect) , 1);
				 #ifdef DEBUG_1
		             experiment->set_n_memb_query();
		         #endif
			 }
		 }
	}	 
		 
	return false; //no new component added 
}


/* Questa funzione simula il parsing del prefisso passato come parametro nel DFA. In realtà 
 * il parsing avviene usando le componenti.Si prende il primo simbolo del prefisso che è un 
 * simbolo dell'alfabeto che per costruzione deve stare in un componente e qi quest ultimo
 * si prende lo short prefix. Poi si concatena lo short prefix trovato al passo precedente
 * con il successivo simbolo del prefisso e si ottiene una nuova stringa della quale bisogna
 * trovare lo short prefix che sicuramente sarà presente per costruzione (perchè il metodo
 * di ricerca del prefisso ricalca il metodo di costruzione del DFA ipotesi dalle componenti
 */
vector<SYMBOL> gi::observationpack::get_ShortPrefixFromComponents(const vector<SYMBOL> &prefix)
{
	vector<SYMBOL> shortPrefix;
	#ifdef DEBUG_3
	    vector<SYMBOL> nilVector = {NIL};
	#endif
	
	for(auto i = prefix.begin(); i!=prefix.end() ; i++)
	{
		shortPrefix.push_back(*i);
		shortPrefix = get_shortPrefix(shortPrefix);
		#ifdef DEBUG_3
		   assert(("Il prefisso non e\' stato trovato in nessuna componente. Cosa impossibile per costruzione",(!(shortPrefix==nilVector))));
		#endif
	}
	
	return shortPrefix;
}

// Si può usare questa funzione in luogo di get_ShortPrefixFromComponents().
//Sulla carta dovrebbe essere più efficiente dato che parsa il prefisso di cui si vuole ottenere lo short prefix
//nel DFA ipotesi per ottenerne l'indice dello stato in cui va a finire. Da questo indice si risale allo short prefix
//se nel momento di creazione dell'ipotesi ci si è memorizzati la corrispondenza tra indice di stato e short prefix
//nella funzione obk_to_dfa. Tuttavia sperimentalmente si è visto che il tempo totale impiegato usando la funzione get_ShortPrefixFromComponents
//è dell'ordine delle decine di nanosecondi invece con get_shortPrefixFromHypothesis è dell'ordine delle migliaia di nanosecondi. Forse
//perchè in entrambi i casi il costo dell'esecuzione dell'algoritmo è talmente infimo che diventa preponderante la chiamata di funzione
//di un oggetto diverso da quello corrente.
//Per di piu' usare get_shortPrefixFromHypothesis comporta di sporcare il codice.  Fare diventare la map states in obk_to_dfa un membro
//della classe e passare il DFA ipotesi alla funzione di gestione del controesempio.
//Pertanto si usa get_ShortPrefixFromComponents()
/*
vector<SYMBOL> gi::observationpack::get_shortPrefixFromHypothesis(dfa *dfahp , const vector<SYMBOL> &prefix)
{
	int index_arriveState = dfahp->get_arrive_state(prefix); //index of the state where prefix arrives with the parsing into DFA hypothesis
	#ifdef DEBUG_3
	    assert(("Sono in get_shortPrefixFromHypothesis.Il prefisso arriva in uno stato non valido",index_arriveState!=ND)); //lancia la assert se get_arrive_state() torna ND
	#endif
	for(auto s=states.begin() ; s!=states.end() ; s++)
		if(s->second[0] == index_arriveState)
		   return s->first; //return the short prefix
	#ifdef DEBUG_3
	    assert(("Sono in get_shortPrefixFromHypothesis. L'indice dello stato d'arrivo nell'ipotesi non corrisponde a nessun short prefix",false)); //lancia la assert se get_arrive_state() torna ND
	#endif	   
		   
}
*/
