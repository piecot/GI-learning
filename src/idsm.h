/*
 * idsm.h
 */


#ifndef IDSM_H_
#define IDSM_H_

#include "bluefringe.h"


using namespace std;



namespace gi
{
/*! \class IDSM
    \brief Class for IDSM inference algorithm.

    Class for IDSM inference algorithm, with all members and methods for start an inference process.
 */
class idsm : public bluefringe
{
private:

	double	alpha = 0.01;
	double	delta = 0.005;
	//double	delta = 1000;


	// STATISTICS
	double error_rate_final_dfa;						/*!< Error rate for the dfa inferred respect to samples*/



	void set_acceptor_and_rejector_states(dfaEDSM* dfa1, vector<SYMBOL>* positive, const int dim_positive, vector<SYMBOL>* negative, const int dim_negative, int* &wp, int* &wn);


	/**
	 * Return a score for a dfa.
	 * It's calculated by statistical heuristic of IDSM, minimimazing the ratio between second type error and earn states.
	 *
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return A score for a DFA, greater is better. "MINF" value if the merge is not acceptable.
	 */
	double merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, double earn_states);

	double merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, int earn_states){};


	int  merge_heuristic_score(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* wp , int* wn){};


	/**
	 * Check if the merge is statisticaly acceptable.
	 *
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return Return if a merge is statistically acceptable.
	 */
	bool is_mergeable(double error_rate_before, double error_rate_after, int dim_strings, double alpha);


	/**
	 * Return the number of strings, of the Training Set, not recognized from final automaton.
	 * Note, that is not a generalization index but an error index, because it's calculated on the
	 * training set rather than test set.
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return
	 */
	static double error_rate(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* &wp, int* &wn, const int tot_wp_w, const int tot_wn_w);

	/**
	 * Return a score for a dfa. It's calculated by heuristic emerged during Abbadingo Competition (1998).
	 *
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return A score for a DFA, greater is better.
	 */
	//int  promotion_heuristic_score(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative);


public:

	/**
	 * Instance an object with all the members and methods for EDSM inference process.
	 * @param path It's the path where find positive and negative samples
	 */
	idsm(const char* path, double alpha_value, double delta_value);							// Take the path where find samples

	/**
	 * Destroy and EDSM object, freeing the memory.
	 */
	~idsm();


	/**
	 * Start an EDSM inference process.
	 * @param path It's the base path where create all the output files of EDSM.
	 * @return Inferred DFA
	 */
	dfa* run(string path);											// Argument is the base path where create files


	double get_error_rate_final_dfa();


	static void compute_ir_stats(dfaEDSM* dfa1, ir_statistical_measures &stats, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* &wp, int* &wn);

};

}


#endif /* IDSM_H_ */
