#include "../../GraphOS/src/matrices/A_matrix.hpp"
#include "../../GraphOS/src/aux_math.hpp"
#include "../../GraphOS/src/matrices/col_vector.tpp"
#include "../../GraphOS/src/matrices/matrix.tpp"
#include "../../GraphOS/src/matrices/symm_matrix.tpp"

unsigned int N_NODES = 3;
#define N_ITERS_INIT (N_NODES*(N_NODES-1)/2*10)
#define N_ITERS_PER_LINK 10
#define N_ITERS (N_NODES*(N_NODES-1)/2*N_ITERS_PER_LINK)
#define INITIALIZE_RANDOMLY false
#define N_DYNAMIC_PAIRS (N_NODES*(N_NODES-1)/2)
#define N_AVRGING 10

#define p  3
#define N_NODES_fin 1000
#define t1_MF 20
#define t2_MF -75
#define t3_MF 25

double sigma2; 
double theta;

// #define MEAN_FIELD

using namespace std;

double N_pairs = N_NODES*(N_NODES-1)/2;

#ifdef MEAN_FIELD
double Hamiltonian(const unsigned int &L) {
  double t_MF_rescaled[p] = {2*t1_MF, 2*t2_MF/N_pairs, 2*t3_MF/(N_pairs*N_pairs)};
  col_vector<double> t_rescaled(p, t_MF_rescaled);
  double H=0;
  for(unsigned int i=0; i<p; ++i)
    H-=t_rescaled[i]*pow(L, i+1); //+1 because t[0]==t1
  return H;
}
#endif

int main() {
  prng rnd(RAND_MAX);
#ifndef MEAN_FIELD
  cout << string("MATRIX_BASED__TRIAD__CONNECTANCE_2_STARS_and_LC__N_NODES_start_") + to_string(N_NODES) + "__N_NODES_fin_" + to_string(N_NODES_fin) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INIT_RANDOM_" + to_string(INITIALIZE_RANDOMLY) + "__N_DYNAMIC_PAIRS_1.csv" << '\n';
  ofstream ofs(string("../data/triad_model/n_scanning/MATRIX_BASED__TRIAD__CONNECTANCE_2_STARS_and_LC__N_NODES_start_") + to_string(N_NODES) + "__N_NODES_fin_" + to_string(N_NODES_fin) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INIT_RANDOM_" + to_string(INITIALIZE_RANDOMLY) + "__N_DYNAMIC_PAIRS_1.csv");
#else
  cout << string("MATRIX_BASED__MF__CONNECTANCE_2_STARS_and_LC__N_NODES_start_") + to_string(N_NODES) + "__N_NODES_fin_" + to_string(N_NODES_fin) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INIT_RANDOM_" + to_string(INITIALIZE_RANDOMLY) + "__N_DYNAMIC_PAIRS_1.csv" << '\n';
  ofstream ofs(string("../data/triad_model/n_scanning/MATRIX_BASED__MF__CONNECTANCE_2_STARS_and_LC__N_NODES_start_") + to_string(N_NODES) + "__N_NODES_fin_" + to_string(N_NODES_fin) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INIT_RANDOM_" + to_string(INITIALIZE_RANDOMLY) + "__N_DYNAMIC_PAIRS_1.csv");
#endif

  for(;N_NODES < N_NODES_fin; ++N_NODES) {
    N_pairs = N_NODES*(N_NODES-1)/2;
    sigma2 = 2*N_NODES/(double)(N_NODES-2.) * t2_MF; // Conversion from the MF couplings to triad
    theta = N_NODES/(double)(N_NODES-2.) * t3_MF; // Conversion from the MF couplings to triad
    
    cerr << N_NODES << '\t';
    A_matrix A(N_NODES);
#ifndef MEAN_FIELD
    if(!INITIALIZE_RANDOMLY) // Initial thermalisation
      A.sample_triad_model_with_single_link_Metropolis(N_ITERS_INIT, rnd, t1_MF, sigma2, theta, INITIALIZE_RANDOMLY);
#else
    if(!INITIALIZE_RANDOMLY) // Initial thermalisation
      A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS_INIT, rnd, N_DYNAMIC_PAIRS);
#endif
      
      for(unsigned int i=0; i<N_AVRGING; ++i) {
#ifdef MEAN_FIELD
        A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS, rnd, INITIALIZE_RANDOMLY);
#else
        A.sample_triad_model_with_single_link_Metropolis(N_ITERS, rnd, t1_MF, sigma2, theta, INITIALIZE_RANDOMLY);
#endif
        ofs << N_NODES << ',' << A.num_links()/N_pairs << ',' << A.average_p_stars(2) << ',' << A.num_triangles() << ',' << A.loc_clust_coeff_col_vec().avrg() << "\n";
        ofs.flush();
      }
  }
    ofs.close();
    return 0;
}
