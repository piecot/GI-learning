/*
 * observationpack.h
 *
 *  Created on: 12 feb 2016
 *      Author: nicola
 */

#ifndef OBSERVATIONPACK_H_
#define OBSERVATIONPACK_H_

#include <string>
#include <map>
#include <vector>
#include <set>

#include "dfa.h"
#include "esperimenti.h"

#define bindclass gi::observationpack 
#define NIL  numeric_limits<SYMBOL>::max()

using namespace std;

enum class nodeType
{
	NODE_TYPE_LEAF   = 1,
	NODE_TYPE_INNER   = 0
};


/*! \class lstar
    \brief Class for OBSERVATION PACK inference algorithm.

    Class for OBSERVATION PACK inference algorithm, with all members and methods for start an inference process through Learner and Teacher.
 */
namespace gi
{
class observationpack
{
    private:
    
    /*! \class component
    \brief Class for store the individual components.
    */
    class component //classe privata accessibile solo da observationpack
    {
		public:
		    //default constructor
		    
		    //FORSE È UNA BUONA IDEA DICHIARARLA COME multimap<vector<SYMBOL>,bool>  DA VALUTARE
		    //The value of map (the second field) can't be a bool. The reason is thin. In split i must erase from the table mq(prefix Concatenated a suffix) because the prefix goes to a new component. But can happen that another prefix created  the same prefix concatenated with suffix. Then I must erase if the second field of the pair  of mq is 1 .
		    map<vector<SYMBOL>, pair<bool,unsigned short int> > 	mq; /*!< Function that realizes the mapping between strings(the first field of mq) and true,false */
		    set<vector<SYMBOL> > pref;              /*!< Prefixes of a component. Each component has own prefixes */
	};
	
	map<vector<SYMBOL> , component* >  components;               /*!< Set of components. There is the access sequence of component and the pointer to component */
	/*Un'altra possibilità che evita la dichiarazione della classe component è
	 * typedef vector<SYMBOL> prefissi;
     * typedef vector<SYMBOL> shortprefix;
     * typedef  map<vector<SYMBOL>,bool>  observationTable;
     * typedef pair<observationTable , prefissi > componente;
     * map<shortprefix ,  componente> componenti;  
     * Dal punto di vista dell'efficenza non credo ci siano grosse differenze
     * dato che comunque viene definito un oggetto pair. Di contro l'accesso
     * risulterebbe più farraginoso*/
     
    //sono globali e uguali per tutte le componenti con ONEGlobally 
    vector< vector<SYMBOL> > suffixes;                          /*!< Suffixes of all components */
    
    //now the data structures for discrimination tree
    vector< pair<vector<SYMBOL> , nodeType> > nodes;            /*!< Nodes of the discrimination tree. There is a label and the nodetype for every node */
    //isn't possible write vector< unsigned short int [2]> >edges because classical array has not copy constructor
    // è un array bidimensionale che puo' crescere come numero di righe e non come colonne (un nodo ha due archi al massimo)
    vector< array< unsigned short int, 2 > > edges;             /*!< Edges for each node of discrimination tree. It' a binary tree then there are two value the first for left child and the second for right child */
                                               
    dfa*    target;	                                            /*!< Target DFA of inference process. Used from Teacher. */
    unsigned short int    dim_alfabeto;							/*!< Alphabet size */	
	vector<SYMBOL>*   	  alfabeto;                             /*!< Alphabet */
	array<vector<SYMBOL>,3> (gi::observationpack::*ptrCounterexampleFunction) (const bool , const vector<SYMBOL> &);	/*!< Pointer to function used for find the suffix(discriminator) from counterexample */		
    
    #ifdef DEBUG_1
        esperimenti* experiment;     /*!< Experiments object */
    #endif
    //#ifdef DEBUG_1
      //  int n_memb_query;                                       /*!< Number of membership query  carried out during the inference process */
        //int n_eq_query;                                         /*!< Number of equivalences query  carried out during the inference process */
        /*//un altro dei parametri è la lunghezza del controesempio di lunghezza massima
    #endif*/ 
    //This function is for make more readable the code.
    //set the second field of mq of the component to value and n
    
    /**
	 * Utility function. This function is for make more readable the code.
	 * 
	 * Set the second field of the \seealso gi::observationpack::component::mq of the component @p mqCompPtr for the string @p key to @p value and @p n
	 * @param mqCompPtr Pointer to a component
	 * @param key String of \seealso gi::observationpack::component::mq
	 * @param value Is 0 or 1, is the result of MQ on the target for @p key string.
	 * @param n is the increment or decrement. It points how many times @p key appears in the \seealso gi::observationpack::component::mq
	 */
    inline void set_mq(map<vector<SYMBOL>,pair<bool,unsigned short int> > *mqCompPtr , const vector<SYMBOL>& key , bool value, unsigned short int n); //I can declare it in the class component but later the implementation must be inserted in this file
   
    public:	
     
        /**
		 * Instantiates an object with all the members and methods for OBSERVATION PACK inference process. It set the rivest-shapire like counterexample decomposition method
		 * @param targetdfa It's the target DFA for the inference process, used from Teacher.
		 */
        observationpack(const dfa &targetdfa); //Con rivest-shapire per default
        
        /**
		 * Instantiates an object with all the members and methods for OBSERVATION PACK inference process
		 * @param targetdfa It's the target DFA for the inference process, used from Teacher.
		 * @param counterexample_function is a string indicating the counterexample decomposition method. It isn't case-sensitive. Admissible values are : 
		 *     - rivest-shapire
		 *     - find-exponential
		 *     - partition-search
		 *     - rivest-shapire-eager 
		 */
        observationpack(const dfa &targetdfa , string counterexample_function);
        
        /**
		 * Destroy an OBSERVATION PACK object, freeing the memory.
		 */
        ~observationpack();
        
        /**
		 * Print the set of components
		 */
        void print_components();
        
        /**
		 * Print the discrimination tree
		 */
        void print_discrimination_tree();
        
        /**
		 * Start an OBSERVATION PACK inference process.
		 * 
		 * Minimize ,in a destructive way, the DFA member \seealso gi::observationpack::target of the calling object
		 * @param approximate If it's true, set an OBSERVATION PACK inference process using an approximation of target DFA.
		 * @param samplestestpath Samples approximating target DFA. Its value must be the empty string "" if @p approximate  is false
		 * @return Return the minimal dfa inferred from OBSERVATION PACK algorithm
		 */
        dfa*  run_observationpack(bool approximate, string samplestestpath);
        //#ifdef DEBUG_1
          //  /**
		   //  * Return the number of membership query in the OBSERVATION PACK inference process
		   //  * @return Return the number of membership query in the OBSERVATION PACK inference process
		   //  */
            //int get_n_memb_query();
            
            // /**
		     //* Return the number of equivalence query in the OBSERVATION PACK inference process
		     //* @return Return the number of equivalence query in the OBSERVATION PACK inference process
		     //*/
            //int get_n_eq_query();
        //#endif
        
        
        /**
		 * Find a decomposition of the @p witness that make unclosed a component
		 * 
		 * For more details <http://jmlr.org/proceedings/papers/v34/isberner14a.html>
		 * @param lambda_H_witness Is true if the @p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param witness Counterexample returned from the teacher
		 * @return Return the decomposition of the counterexample in a array struct of 3 elements: In order you have the first prefix, the second prefix that is an access sequence in some component (both prefixes are part of the same component), the suffix which discriminates the 2 prefixes
		 */
        array<vector<SYMBOL>,3> rivest_shapire(const bool lambda_H_witness,const vector<SYMBOL> &witness);
        
        /**
		 * Find a decomposition of the @p witness that make unclosed a component
		 * 
		 * For more details <http://jmlr.org/proceedings/papers/v34/isberner14a.html>
		 * @param lambda_H_witness Is true if the @p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param witness Counterexample returned from the teacher
		 * @return Return the decomposition of the counterexample in a array struct of 3 elements: In order you have the first prefix, the second prefix that is an access sequence in some component (both prefixes are part of the same component), the suffix which discriminates the 2 prefixes
		 */
        array<vector<SYMBOL>,3> find_exponential(const bool lambda_H_witness,const vector<SYMBOL> &witness);
        
        /**
		 * Find a decomposition of the @p witness that make unclosed a component
		 * @param lambda_H_witness Is true if the @p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param witness Counterexample returned from the teacher
		 * @return Return the decomposition of the counterexample in a array struct of 3 elements: In order you have the first prefix, the second prefix that is an access sequence in some component (both prefixes are part of the same component), the suffix which discriminates the 2 prefixes
		 */
        array<vector<SYMBOL>,3> partition_search(const bool lambda_H_witness,const vector<SYMBOL> &witness);
        
        /**
		 * Find a decomposition of the \p witness that make unclosed a component
		 * 
		 * For more details <http://jmlr.org/proceedings/papers/v34/isberner14a.html>
		 * @param lambda_H_witness Is true if the \p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param witness Counterexample returned from the teacher
		 * @return Return the decomposition of the counterexample in a array struct of 3 elements: In order you have the first prefix, the second prefix that is an access sequence in some component (both prefixes are part of the same component), the suffix which discriminates the 2 prefixes
		 */
        array<vector<SYMBOL>,3> rivest_shapire_eager(const bool lambda_H_witness,const vector<SYMBOL> &witness);
        
        /**
		 * Find a index i such that the symbols from i+1 onwards of witness are a suffix that make unclosed a component.
		 * 
		 * For more details <http://jmlr.org/proceedings/papers/v34/isberner14a.html>
		 * @param  low Index of the @p witness where start the binary search.
		 * @param  hight Index of the @p witness where end the binary search.
		 * @param  lambda_H_witness Is true if the @p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param  witness Counterexample returned from the teacher
		 * @return Return a index i such that the symbols from i+1 onwards of witness are a suffix that make unclosed a component.
		 */
        int rivest_shapire_classic(int low , int hight , const bool lambda_H_witness ,const vector<SYMBOL> &witness);
        
        
        /**
		 * Utility function allowing you to find a row of the observation table for a given prefix and component.
		 * 
		 * E.g. if the prefix is pr and the suffixes are j and jk the row builted is @f$ MQ[prj] \cdot MQ[prjk] @f$ where MQ are
		 * the membership queries realized in the target DFA 
		 * @param  low Index of the @p witness where start the binary search.
		 * @param  hight Index of the @p witness where end the binary search.
		 * @param  lambda_H_witness Is true if the @p witness ends in an accepting state of the  hypothesis DFA H.
		 * @param  witness Counterexample returned from the teacher
		 * @return Return a index i such that the symbols from i+1 onwards of witness are a suffix that make unclosed a component.
		 */
        vector<SYMBOL> get_row(vector<SYMBOL> prefix , map<vector<SYMBOL>, pair<bool,unsigned short int> > 	*mqComp);
        
        /**
		 * Find the component where is contained the \p prefix.
		 * @param  prefix The prefix that you want to find the component.
		 * @return Return the access sequence of the component where is contained the @p prefix. Return {NIL} if is not present.
		 */
        vector<SYMBOL> get_shortPrefix(const vector<SYMBOL> &prefix);
		 
		/**
		 * Find the short prefix corresponding at the state in the hypothesis where arrives the \p prefix.
		 * 
		 * Use this method when you do not have the assurance that the prefix is contained in the set of components (else you can use \seealso gi::observationpack::get_shortPrefix() ) 
		 * This method simulates the building of the hypothesis on the fly 
		 * @param  prefix The prefix that you want to find the access sequence.
		 * @return Return the access sequence of the state in the hypothesis where arrive the @p prefix.
		 */ 
        vector<SYMBOL> get_ShortPrefixFromComponents(const vector<SYMBOL> &prefix); 
        //vector<SYMBOL> get_shortPrefixFromHypothesis(dfa *dfahp , const vector<SYMBOL> &prefix);
        
        
        /**
		 * Add the suffix returned from decomposition counterexample function to  \seealso gi::observationpack::suffixes and complete components
		 * 
		 * Moreover the MQ for the 2 prefixes returned from decomposition counterexample function aren't effectuated because the outcome is already known.
		 * @param  decomposition The decomposition of the counterexample(witness) obtained from a method of decomposition counterexample.
		 * @param  lambda_H_witness Is true if the  counterexample(witness) ends in an accepting state of the  hypothesis DFA H. Here is used for no makes the MQ for the 2 prefixes of the decomposition too.
		 */   
        void update_from_counterexample(const array<vector<SYMBOL>,3> &decomposition , bool lambda_H_witness);
        
        /**
		 * Build a DFA from \seealso gi::observationpack::components.
		 * @return Pointer to builded DFA. 
		 */
        dfa* obk_to_dfa();
        
        /**
		 * Make closed the observation pack.
		 * @param decomposition From counterexample decomposition function you know which prefixes make not closed the components' set: the two prefixes of the decomposition by means the siffix (the third value of @p decomposition).The first time in absolute that \seealso gi::observationpack::closePack() is calles \p decomposition has to be {{{NIL},{NIL},{NIL}}}
		 */
        void closePack(const array<vector<SYMBOL>,3> &decomposition);
        
        /**
		 * Modify the observation pack splitting the component and the leaf of the discrimination tree where there is unclosdeness.
		 * @param shortPr Identify the component where there is unclosdness. Is a prefix with which \seealso gi::observationpack::get_row(@p shortPr) \f$ \neq \f$ \seealso gi::observationpack::get_row(@p pr)
		 * @param pr Prefix that holds unclosdness with shortPr
		 * @param suf the suffix that discriminates the two prefixes  
		 */
        void split(const vector<SYMBOL> &shortPr ,const vector<SYMBOL> &pr ,const vector<SYMBOL> &suf); //this funciont can be void  rather than return a component because in closePack you save the access sequences in W(and you know already the access sequence in closePack)
        
        /**
		 * Sink a prefix in the observation pack
		 * @param pr Prefix inserted in the observation pack
		 * @return Return true if is created a new component (and leaf) else false 
		 */
        bool sift(const vector<SYMBOL> &pr);
        
        
};
}
#endif
