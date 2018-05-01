//
//  rsa.cpp
//  KRY-RSA
//
//  Created by Jordán Jarolím on 30.04.18.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#include "rsa.hpp"

using namespace std;


void Rsa::getN(mpz_t N, mpz_t p, mpz_t q){
    mpz_mul(N, p, q);
}

void Rsa::getPhi(mpz_t phi, mpz_t p, mpz_t q){
    mpz_t pMinus; mpz_init(pMinus);
    mpz_t qMinus; mpz_init(qMinus);
    mpz_sub_ui(pMinus, p, 1);
    mpz_sub_ui(qMinus, q, 1);
    mpz_mul(phi, pMinus, qMinus);
}
