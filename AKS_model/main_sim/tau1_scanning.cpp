#include "../../GraphOS/src/matrices/A_matrix.hpp"
#include "../../GraphOS/src/aux_math.hpp"
#include "../../GraphOS/src/matrices/col_vector.tpp"
#include "../../GraphOS/src/matrices/matrix.tpp"
#include "../../GraphOS/src/matrices/symm_matrix.tpp"

#define N_NODES 200
#define N_ITERS_INIT (N_NODES*(N_NODES-1)/2*10)
#define N_ITERS_PER_LINK 100
#define N_ITERS (N_NODES*(N_NODES-1)/2*N_ITERS_PER_LINK)
#define N_DYNAMIC_PAIRS 1 //(N_NODES*(N_NODES-1)/2)
#define N_AVRGING 10
#define INITIALIZE_RANDOMLY false

#define tau1_step .1
#define tau1_start 1.5
#define tau1_fin 11.5

#define mu -200
#define lambda .04


using namespace std;

double t[3] = {tau1_start, mu, lambda};

double N_pairs = N_NODES*(N_NODES-1)/2;


int main() {
  prng rnd(RAND_MAX);

  A_matrix A(N_NODES);
  
  cout << string("AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__tau1_start_" + to_string(tau1_start) + "__tau1_fin_" + to_string(tau1_fin) + "__tau1_step_" + to_string(tau1_step) + "__mu_" + to_string(mu) + "__lambda_" + to_string(lambda) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';

  
  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_AKS_model_with_single_link_Metropolis(N_ITERS_INIT, rnd, t[0], mu, lambda, INITIALIZE_RANDOMLY);
  
  ofstream ofs(string("../data/AKS_model/tau1_scanning/AKS___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__tau1_start_" + to_string(tau1_start) + "__tau1_fin_" + to_string(tau1_fin) + "__tau1_step_" + to_string(tau1_step) + "__mu_" + to_string(mu) + "__lambda_" + to_string(lambda) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
  for(t[0]=tau1_start; t[0]<tau1_fin; t[0]+=tau1_step) { // t[0] is tau1 from the draft
    cerr << t[0] << '\t';
    for(unsigned int i=0; i<N_AVRGING; ++i) {
      A.sample_AKS_model_with_single_link_Metropolis(N_ITERS, rnd, t[0], mu, lambda, INITIALIZE_RANDOMLY);
      ofs << t[0] << ',' << A.num_links()/N_pairs << ',' << A.average_p_stars(2) << ',' << A.num_triangles() << ',' << A.loc_clust_coeff_col_vec().avrg() << "\n";
      ofs.flush();
    }
  }
  ofs.close();
  return 0;
}
