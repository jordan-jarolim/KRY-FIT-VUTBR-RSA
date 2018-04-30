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
