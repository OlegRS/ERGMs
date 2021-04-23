#include "../../GraphOS/src/matrices/A_matrix.hpp"
#include "../../GraphOS/src/aux_math.hpp"
#include "../../GraphOS/src/matrices/col_vector.tpp"
#include "../../GraphOS/src/matrices/matrix.tpp"
#include "../../GraphOS/src/matrices/symm_matrix.tpp"


#define N_NODES 1000
#define N_ITERS_INIT (N_NODES*(N_NODES-1)/2*10)
#define N_ITERS_PER_LINK 10
#define N_ITERS (N_NODES*(N_NODES-1)/2*N_ITERS_PER_LINK)
#define N_DYNAMIC_PAIRS 1 //(N_NODES*(N_NODES-1)/2)
#define N_AVRGING 10
#define INITIALIZE_RANDOMLY false

#define p  3
#define t1_step .01
#define t1_start -5.5
#define t1_fin -1
#define t2_MF 2
#define t3_MF 1

#define sigma2 (2*N_NODES/(double)(N_NODES-2.) * t2_MF) // Conversion from the MF couplings to triad
#define theta  (N_NODES/(double)(N_NODES-2.) * t3_MF) // Conversion from the MF couplings to triad

//#define MEAN_FIELD

using namespace std;

double t[p] = {t1_start, sigma2, theta};
double t_MF[p] = {t1_start, t2_MF, t3_MF};

double N_pairs = N_NODES*(N_NODES-1)/2;


#ifdef MEAN_FIELD
double t_MF_rescaled[p] = {2*t1_start, 2*t2_MF/N_pairs, 2*t3_MF/(N_pairs*N_pairs)};
col_vector<double> t_rescaled(p, t_MF_rescaled);

double Hamiltonian(const unsigned int &L) {
  t_rescaled[0] = 2*t[0];
  double H=0;
  for(unsigned int i=0; i<p; ++i)
    H-=t_rescaled[i]*pow(L, i+1); //+1 because t[0]==t1
  return H;
}
#endif

int main() {
  prng rnd(RAND_MAX);

  A_matrix A(N_NODES);
#ifndef MEAN_FIELD
  cout << string("TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';

  
  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_triad_model_with_single_link_Metropolis(N_ITERS_INIT, rnd, t[0], sigma2, theta, INITIALIZE_RANDOMLY);
  
  ofstream ofs(string("../data/triad_model/t1_scanning/TRIAD___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
#else
  cout << string("MF___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_DYNAMIC_PAIRS_" + to_string(N_DYNAMIC_PAIRS) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';


  if(!INITIALIZE_RANDOMLY)
    A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS_INIT, rnd, N_DYNAMIC_PAIRS);

  ofstream ofs(string("../data/mean_field_model/t1_scanning/MF___CONNECTANCE__2_STARS__N_TRIANGLES__LCC_AVRG___") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_DYNAMIC_PAIRS_" + to_string(N_DYNAMIC_PAIRS) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__ENSEMBLE_AVRG_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
#endif
  for(t[0]=t1_start; t[0]<t1_fin; t[0]+=t1_step) { // t[0] is t1 from the draft
    cerr << t[0] << '\t';
    for(unsigned int i=0; i<N_AVRGING; ++i) {
#ifdef MEAN_FIELD
      A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS, rnd, INITIALIZE_RANDOMLY);
#else
      A.sample_triad_model_with_single_link_Metropolis(N_ITERS, rnd, t[0], sigma2, theta, INITIALIZE_RANDOMLY);
#endif
      ofs << t[0] << ',' << A.num_links()/N_pairs << ',' << A.average_p_stars(2) << ',' << A.num_triangles() << ',' << A.loc_clust_coeff_col_vec().avrg() << "\n";
      ofs.flush();
    }
  }
  ofs.close();
  return 0;
}
