#!/bin/bash

rm ../bin/*

[ -d "../data" ] || mkdir "../data"

[ -d "../data/triad_model" ] || mkdir "../data/triad_model"

[ -d "../data/triad_model/t1_scanning" ] || mkdir "../data/triad_model/t1_scanning"
[ -d "../data/triad_model/t2_scanning" ] || mkdir "../data/triad_model/t2_scanning"
[ -d "../data/triad_model/t3_scanning" ] || mkdir "../data/triad_model/t3_scanning"
[ -d "../data/triad_model/n_scanning" ] || mkdir "../data/triad_model/n_scanning"
[ -d "../data/triad_model/degree_distributions" ] || mkdir "../data/triad_model/degree_distributions"
[ -d "../data/triad_model/clustering_distributions" ] || mkdir "../data/triad_model/clustering_distributions"

[ -d "../data/mean_field_model" ] || mkdir "../data/mean_field_model"

[ -d "../data/mean_field_model/t1_scanning" ] || mkdir "../data/mean_field_model/t1_scanning"
[ -d "../data/mean_field_model/t2_scanning" ] || mkdir "../data/mean_field_model/t2_scanning"
[ -d "../data/mean_field_model/t3_scanning" ] || mkdir "../data/mean_field_model/t3_scanning"
[ -d "../data/mean_field_model/n_scanning" ] || mkdir "../data/mean_field_model/n_scanning"
[ -d "../data/mean_field_model/degree_distributions" ] || mkdir "../data/mean_field_model/degree_distributions"
[ -d "../data/mean_field_model/clustering_distributions" ] || mkdir "../data/mean_field_model/clustering_distributions"

[ -d "../bin" ] || mkdir "../bin"


g++ -std=c++11 t1_scanning.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp  -Ofast -o ../bin/t1_scanning
g++ -std=c++11 ensemble_avrging_t1_scan.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp  -Ofast -o ../bin/ensemble_avrging_t1_scan
g++ -std=c++11 degree_and_cc_distr.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp  -Ofast -o ../bin/degree_and_cc_distr
g++ -std=c++11 n_scanning.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp  -Ofast -o ../bin/n_scanning
