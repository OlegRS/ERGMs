#!/bin/bash

rm ../bin/*

[ -d "../data" ] || mkdir "../data"

[ -d "../data/p_star_model" ] || mkdir "../data/p_star_model"

[ -d "../data/p_star_model/t1_scanning" ] || mkdir "../data/p_star_model/t1_scanning"
[ -d "../data/p_star_model/t2_scanning" ] || mkdir "../data/p_star_model/t2_scanning"
[ -d "../data/p_star_model/t3_scanning" ] || mkdir "../data/p_star_model/t3_scanning"
[ -d "../data/p_star_model/degree_distributions" ] || mkdir "../data/p_star_model/degree_distributions"
[ -d "../data/p_star_model/clustering_distributions" ] || mkdir "../data/p_star_model/clustering_distributions"

[ -d "../data/mean_field_model" ] || mkdir "../data/mean_field_model"

[ -d "../data/mean_field_model/t1_scanning" ] || mkdir "../data/mean_field_model/t1_scanning"
[ -d "../data/mean_field_model/t2_scanning" ] || mkdir "../data/mean_field_model/t2_scanning"
[ -d "../data/mean_field_model/t3_scanning" ] || mkdir "../data/mean_field_model/t3_scanning"
[ -d "../data/mean_field_model/degree_distributions" ] || mkdir "../data/mean_field_model/degree_distributions"
[ -d "../data/mean_field_model/clustering_distributions" ] || mkdir "../data/mean_field_model/clustering_distributions"

[ -d "../bin" ] || mkdir "../bin"


g++ -std=c++11 t1_scanning.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp -Ofast -o ../bin/t1_scanning

g++ -std=c++11 deg_distrib.cpp  ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp -Ofast -o ../bin/deg_distrib

g++ -std=c++11 degree_and_lcc_distributions.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp  -Ofast -o ../bin/degree_and_lcc_distributions
