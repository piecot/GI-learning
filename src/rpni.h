

#ifndef RPNI_H_
#define RPNI_H_

#include "bluefringe.h"


using namespace std;


namespace gi
{
/*! \class rpni
    \brief Class for RPNI inference algorithm.

    Class for RPNI inference algorithm, with all members and methods for start an inference process.
 */
class rpni : public bluefringe
{
private:


	/**
	 * NOT USED IN RPNI! DELETE IT! !?!?!?!?! Return a score for a dfa. It's calculated by heuristic emerged during Abbadingo Competition (1998).
	 *
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return A score for a DFA, greater is better.
	 */
	int merge_heuristic_score(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* wp , int* wn);

	double merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, int earn_states){};

	double merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, double earn_states){};


public:

	/**
	 * Instance an object with all the members and methods for RPNI inference process.
	 * @param path It's the path where find positive and negative samples
	 */
	rpni(const char* path);							// Take the path where find samples

	/**
	 * Destroy and RPNI object, freeing the memory.
	 */
	~rpni();


	/**
	 * Start an RPNI inference process.
	 * @param path It's the base path where create all the output files of RPNI.
	 * @return Inferred DFA
	 */
	dfa* run(string path);						// Argument is the base path where create files


	double run_elapsed_time(string base_path, dfa** res);


};

}


#endif /* RPNI_H_ */
