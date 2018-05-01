//
//  rsa.hpp
//  KRY-RSA
//
//  Created by Jordán Jarolím on 30.04.18.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#ifndef rsa_hpp
#define rsa_hpp

#include <stdio.h>
#include <gmp.h>


class Rsa{
public:
    void getN(mpz_t N, mpz_t p, mpz_t q);
    void getPhi(mpz_t phi, mpz_t p, mpz_t q);
};
#endif /* rsa_hpp */
