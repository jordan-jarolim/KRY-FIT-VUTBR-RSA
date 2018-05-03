//
//  solovay.hpp
//  KRY-RSA
//
//  Created by Jordán Jarolím on 28.04.18.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#ifndef solovay_hpp
#define solovay_hpp

#include <stdio.h>
#include <tuple>
#include <gmp.h>

typedef struct {
    mpz_t p;
    mpz_t q;
} tRandoms;


class Solovay{
public:
    int length;
    Solovay(int length);
    tRandoms getPrimes(gmp_randstate_t rstate);
    
private:
    void getRandom(mpz_t p, mpz_t q, gmp_randstate_t rstate);
    bool testPrime(mpz_t prime);
    int jacob(mpz_t a, mpz_t n);
};
#endif /* solovay_hpp */
