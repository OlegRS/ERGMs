///////////////////////////////////////////////////////////
// This code saves connectances and the numbers of stars //
// for the different values of t1 at fixed t2 and t3.    //
///////////////////////////////////////////////////////////

#include "../../GraphOS/src/matrices/A_matrix.hpp"
#include "../../GraphOS/src/aux_math.hpp"
#include "../../GraphOS/src/matrices/col_vector.tpp"
#include "../../GraphOS/src/matrices/matrix.tpp"
#include "../../GraphOS/src/matrices/symm_matrix.tpp"


#define N_NODES 1000
#define N_ITERS_INIT (N_NODES*(N_NODES-1)/2*100)
#define N_ITERS_PER_LINK 10
#define N_ITERS (N_NODES*(N_NODES-1)/2*N_ITERS_PER_LINK)
#define INITIALIZE_RANDOMLY false
#define N_DYNAMIC_PAIRS 1 //(N_NODES*(N_NODES-1)/2)
#define N_AVRGING 10

#define p  3
#define t1_step .01
#define t1_start 2
#define t1_fin 10
#define t2_MF -15.
#define t3_MF 9.

#define t2 (2*N_NODES/(double)(N_NODES-2) * t2_MF) // Conversion from MF to PS couplings
#define t3 (6*N_NODES*N_NODES/(double)(N_NODES*N_NODES-5*N_NODES+6) * t3_MF) // Conversion from MF to PS couplings

//#define MEAN_FIELD

using namespace std;

double t[p] = {t1_start, t2, t3};
double N_pairs = N_NODES*(N_NODES-1)/2.;


#ifdef MEAN_FIELD
double t_MF[p] = {t1_start, t2_MF, t3_MF};
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
  cout << string("PS__CONNECTANCE_and_STARS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALIZE_RANDOMLY) + ".csv"<< '\n';

  
  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_p_star_model_with_single_link_Metropolis(N_ITERS_INIT, rnd, col_vector<double>(p,t), INITIALIZE_RANDOMLY);
  
  ofstream ofs(string("../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
#else
  cout << string("MF__CONNECTANCE_and_STARS__") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_DYNAMIC_PAIRS_" + to_string(N_DYNAMIC_PAIRS) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';


  if(!INITIALIZE_RANDOMLY)
    A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS_INIT, rnd, N_DYNAMIC_PAIRS);

  ofstream ofs(string("../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__") + "N_NODES_" + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_start_" + to_string(t1_start) + "__t1_fin_" + to_string(t1_fin) + "__t1_step_" + to_string(t1_step) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_DYNAMIC_PAIRS_" + to_string(N_DYNAMIC_PAIRS) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
#endif
  for(t[0]=t1_start; t[0]<t1_fin; t[0]+=t1_step) { // t[0] is t1 from the draft
    cerr << t[0] << '\t';
    for(unsigned int i=0; i<N_AVRGING; ++i) {
#ifdef MEAN_FIELD
      A.MF_GB_Metropolis_generator(Hamiltonian, N_ITERS, rnd, N_DYNAMIC_PAIRS, INITIALIZE_RANDOMLY);
#else
      A.sample_p_star_model(N_ITERS, rnd, col_vector<double>(p,t), N_DYNAMIC_PAIRS, INITIALIZE_RANDOMLY);
#endif
      ofs << t[0] << ',' << A.num_links()/N_pairs;
      for(unsigned int s=2; s<=p; ++s)
        ofs << ',' << A.average_p_stars(s);
      ofs << '\n';
      ofs.flush();
    }
  }
  ofs.close();
  return 0;
}
