#!/bin/bash

rm ../bin/*

[ -d "../data" ] || mkdir "../data"

[ -d "../data/AKS_model" ] || mkdir "../data/AKS_model"

[ -d "../data/AKS_model/tau1_scanning" ] || mkdir "../data/AKS_model/tau1_scanning"
[ -d "../data/AKS_model/mu_scanning" ] || mkdir "../data/AKS_model/mu_scanning"
[ -d "../data/AKS_model/lambda_scanning" ] || mkdir "../data/AKS_model/lambda_scanning"
[ -d "../data/AKS_model/degree_distributions" ] || mkdir "../data/AKS_model/degree_distributions"
[ -d "../data/AKS_model/clustering_distributions" ] || mkdir "../data/AKS_model/clustering_distributions"


[ -d "../bin" ] || mkdir "../bin"


g++ -std=c++11 tau1_scanning.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp -Ofast -o ../bin/tau1_scanning
g++ -std=c++11 degree_and_lcc_distributions.cpp ../../GraphOS/src/aux_math.cpp ../../GraphOS/src/matrices/A_matrix.cpp -Ofast -o ../bin/degree_and_lcc_distributions

