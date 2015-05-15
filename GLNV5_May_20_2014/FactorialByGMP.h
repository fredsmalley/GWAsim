//
//  FactorialByGMP.h
//  GLNV5_Apr_14_2014
//
//  Created by Hua Zhong on 5/19/14.
//  Copyright (c) 2014 ZH. All rights reserved.
//
//Hua added, May 19 2014
#ifndef GLNV5_Apr_14_2014_FactorialByGMP_h
#define GLNV5_Apr_14_2014_FactorialByGMP_h

#include <gmp.h>
#include <gmpxx.h>

void initialize_mpf_t(mpf_t p, int precision){
    mpf_init_set_ui(p, 1.0);
    mpf_class (p, precision);
};

void factorial_GMP(mpf_t p, int n){
    for (size_t i=1; i <= n ; i++){
        mpf_mul_ui(p,p,i); /* p = p * i */
    }
};

double factorialFormatToDouble(mpf_t p){
    return mpf_get_d(p);
};

#endif
