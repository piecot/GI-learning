#include <iostream>
#include <fstream>
#include "esperimenti.h"

gi::esperimenti::esperimenti() {}

/*gi::esperimenti::esperimenti(const string tipo) : esperimenti(tipo,0,0)
{
	setType(tipo);
	n_memb_query = 0;
	n_eq_query=0;
}*/

gi::esperimenti::esperimenti(int n_states,int dim_alphabet , unsigned long int nmq , int neq)
{
	setNumState(n_states);
	setDimAlphabet(dim_alphabet);
	n_memb_query = nmq;
	n_eq_query=neq;
}

//copy constructor
/*gi::esperimenti::esperimenti(const esperimenti &e1) :esperimenti(e1.type,e1.n_memb_query,e1.n_eq_query)
{
}*/

/*const gi::esperimenti &gi::esperimenti::operator=(const esperimenti &right_exp)
{
	if(&right_exp != this)
	{
		type = right_exp.type;
		n_memb_query = right_exp.n_memb_query;
		n_eq_query = right_exp.n_eq_query;
		size_counterexamples = right_exp.size_counterexamples; 
		
		if(! this->getType().compare("observationpack"))
		    size_suffixes = right_exp.size_suffixes;    
	}
	
	return *this;
}*/

void gi::esperimenti::setNumState(int n_states)
{
	numState = n_states;
}

int gi::esperimenti::getNumState()
{
	return numState;
}

void gi::esperimenti::setDimAlphabet(int da)
{
	dimAlphabet = da;
}

int gi::esperimenti::getDimAlphabet()
{
	return dimAlphabet;
}

void gi::esperimenti::set_n_memb_query ()
{
	++n_memb_query;
}

void gi::esperimenti::set_n_eq_query ()
{
	++n_eq_query;
}

void gi::esperimenti::set_n_memb_query (int num)
{
	n_memb_query += num;
}

unsigned long int gi::esperimenti::get_n_memb_query ()
{
	return n_memb_query;
}

int gi::esperimenti::get_n_eq_query ()
{
	return n_eq_query;
}

void gi::esperimenti::setType (const string tipo)
{
	if(!tipo.compare("rivest-shapire") || !tipo.compare("rivest-shapire-eager") || !tipo.compare("partition-search") || !tipo.compare("find-exponential") || !tipo.compare("lstar") )
		type = tipo;
	else
	{
	    cerr<<"Wrong class of experiments"<<endl;
		exit(EXIT_FAILURE);
	}
}

string gi::esperimenti::getType()
{
	return type;
}

void gi::esperimenti::setSizeCounterexamples(int length)
{
	size_counterexamples.push_back(length);
}

void gi::esperimenti::setSizeSuffixes(int length)
{
	size_suffixes.push_back(length);
}

void gi::esperimenti::myPrint(const string file)
{
        ofstream results;
        results.open(file, std::ios_base::app);
        results << getType() <<" "<<getNumState()<<" "<<getDimAlphabet()<<endl<<get_n_memb_query()<<" "<<get_n_eq_query()<<endl;
        for(auto i=size_counterexamples.begin(); i != size_counterexamples.end() ; i++)
        {
			results<<*i<<" ";
		}
		results<<endl;
		
		if(type.compare("lstar")) //se diverso la lstar
		{
		    for(auto i=size_suffixes.begin(); i != size_suffixes.end() ; i++)
            {
			    results<<*i<<" ";
		    }
		    results<<endl;
	    }
           
}

void gi::esperimenti::writeTime(const string file,double t)
{
	    ofstream results;
        results.open(file, std::ios_base::app);
        results<<t<<endl;	
}
