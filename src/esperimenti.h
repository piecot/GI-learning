#ifndef ESPERIMENTI_H_
#define ESPERIMENTI_H_

#include <string>
#include <vector>

using namespace std;

namespace gi
{
class esperimenti
{
	private:
	    unsigned long int n_memb_query;                                       
        int n_eq_query;
        int numState;
        int dimAlphabet;
        string type;
        vector<int> size_counterexamples; //store the dimensions of the counterexamples 
        vector<int> size_suffixes;
         
	public:
	    esperimenti();
	    //esperimenti(const string tipo);
	    esperimenti(int n_states,int dim_alphabet,unsigned long int nmq , int neq);
	    //esperimenti(const esperimenti &e);
	    void setSizeCounterexamples(int length);
	    void setSizeSuffixes(int length);
	    //const esperimenti &operator=(const esperimenti &right_exp);
	    void set_n_memb_query();
	    void set_n_memb_query(int n);
	    void set_n_eq_query(); 
	    void setNumState(int);
	    void setDimAlphabet(int da);
	    void setType(const string s);
	    void writeTime(const string file,double t);
	    unsigned long int get_n_memb_query();
	    int get_n_eq_query();
	    int getNumState();
	    int getDimAlphabet();
	    string getType();
	    void myPrint(const string file);
};
}

#endif
