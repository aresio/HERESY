#!/bin/bash
nvcc kernel.cu RSParser.cpp -std=c++11 -gencode=arch=compute_20,code=compute_20 -gencode=arch=compute_30,code=compute_30 -gencode=arch=compute_50,code=compute_50  -gencode=arch=compute_35,code=compute_35 -O3 -o RSsim --use_fast_math  --cudart static


