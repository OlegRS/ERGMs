#include "../../GraphOS/src/graph.hpp"
#include <ctime>

#define N_NODES 500
#define N_PAIRS (N_NODES*(N_NODES-1)/2.)
#define N_ITERS_PER_LINK 10
#define N_ITERS ((N_NODES*(N_NODES-1))/2*N_ITERS_PER_LINK)

#define Ld .3
#define t2_MF -75
#define t3_MF 24
#define t1_MF ( 1/2.*log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld )

#define h     t1_MF
#define sigma (N_NODES/(double)(N_NODES-1) * t2_MF) // Conversion from the MF couplings to triad
#define tau   (N_NODES/(double)(N_NODES-2) * t3_MF) // Conversion from the MF couplings to triad

#define INITIALIZE_RANDOMLY false

// #define MEAN_FIELD

using namespace std;;

#ifdef MEAN_FIELD
double t_rescaled[3] = {2*t1_MF, 2*t2_MF/N_PAIRS, 2*t3_MF/(N_PAIRS*N_PAIRS)};

double Hamiltonian(const unsigned int &L) {
  return -t_rescaled[0]*L - t_rescaled[1]*L*L - t_rescaled[2]*L*L*L;
}

#endif

int main() {
  prng rnd(RAND_MAX);
#ifdef MEAN_FIELD
  col_vector<double> params(3, (double[3]){t1_MF, t2_MF, t3_MF});
  cout << string("MF__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';

  string str1 = "mean_field_model";
  string str2 = "MF";
  
  A_matrix A(N_NODES);

  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS*10, rnd);
  

  A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS, rnd, INITIALIZE_RANDOMLY);

#else
  col_vector<double> params(3, (double[3]){h, sigma, tau});
  cout << params << '\n';

  cout << string("MATRIX_BASED_TRIAD__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';

  string str1 = "triad_model";
  string str2 = "TRIAD";

  A_matrix A(N_NODES);

  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_triad_model_with_single_link_Metropolis(N_ITERS*10, rnd, h, sigma, tau);
  
  
  A.sample_triad_model_with_single_link_Metropolis(N_ITERS, rnd, h, sigma, tau, INITIALIZE_RANDOMLY);
#endif

  std::time_t time = std::time(nullptr);
  // Saving the graph in GraphOS format
  graph(A).save(string("../data/") + str1 + "/graphs/" + str2 + "__N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + "__TIME_" + to_string(time));
  // Saving the graph in Cytoscape recognised format
  graph(A).save_to_sif_and_attrs(string("../data/") + str1 + "/graphs/" + str2 + "__N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + "__TIME_" + to_string(time));

  return 0;
}
