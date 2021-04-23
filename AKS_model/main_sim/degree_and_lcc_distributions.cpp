#include "../../GraphOS/src/matrices/A_matrix.hpp"
#include "../../GraphOS/src/aux_math.hpp"
#include "../../GraphOS/src/matrices/col_vector.tpp"
#include "../../GraphOS/src/matrices/matrix.tpp"
#include "../../GraphOS/src/matrices/symm_matrix.tpp"

#define N_NODES 200
#define N_PAIRS (N_NODES*(N_NODES-1)/2.)
#define N_ITERS_PER_LINK 100
#define N_ITERS ((N_NODES*(N_NODES-1))/2*N_ITERS_PER_LINK)
#define N_AVRGING 10000

#define Ld .3
#define mu 3
#define lambda 1
//#define tau1 ( 1/2.*log(Ld/(1-Ld)) - mu*lambda*(1-pow(1-Ld/(N_NODES*lambda), N_NODES-2)) )
#define tau1 -1.1950138739587435

#define INITIALIZE_RANDOMLY false

using namespace std;

int main() {
  prng rnd(RAND_MAX);
  cout << string("AKS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__tau1_" + to_string(tau1) + "__mu_" + to_string(mu) + "__lambda_" + to_string(lambda) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';
    
  ofstream ofs_deg(string("../data/AKS_model/degree_distributions/AKS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__tau1_" + to_string(tau1) + "__mu_" + to_string(mu) + "__lambda_" + to_string(lambda) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  ofstream ofs_cc(string("../data/AKS_model/clustering_distributions/AKS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__tau1_" + to_string(tau1) + "__mu_" + to_string(mu) + "__lambda_" + to_string(lambda) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");

  A_matrix A(N_NODES);
  A.set_Erdos_Renyi(Ld, rnd);

  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_AKS_model_with_single_link_Metropolis(N_ITERS*10, rnd, tau1, mu, lambda);
  
  for(unsigned int i=0; i<N_AVRGING; ++i) {
    A.sample_AKS_model_with_single_link_Metropolis(N_ITERS, rnd, tau1, mu, lambda, INITIALIZE_RANDOMLY);
    ofs_deg << A.degree_sequence_col_vec();
    ofs_deg.flush();
    ofs_cc << A.loc_clust_coeff_col_vec();
    ofs_cc.flush();
    cerr << i << '\t'; 
  }
  
  ofs_deg.close();
  ofs_cc.close();
  return 0;
}
