//
//  rsa.cpp
//  KRY-RSA
//
//  Created by Jordán Jarolím on 30.04.18.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#include "rsa.hpp"
#include <iostream>
#include <string.h>



using namespace std;

/**
 * Get N = p*q
 */
void Rsa::getN(mpz_t N, mpz_t p, mpz_t q){
    mpz_mul(N, p, q);
}

/**
 * Get Phi=(p-1)*(q-1)
 */
void Rsa::getPhi(mpz_t phi, mpz_t p, mpz_t q){
    mpz_t pMinus; mpz_init(pMinus);
    mpz_t qMinus; mpz_init(qMinus);
    mpz_sub_ui(pMinus, p, 1);
    mpz_sub_ui(qMinus, q, 1);
    mpz_mul(phi, pMinus, qMinus);
}

/**
 * e : gcd(e, phi) = 1 && e > 1 && e < phi
 */
void Rsa::getPublicExponent(mpz_t publicExponent, mpz_t phi, gmp_randstate_t rstate){
    bool isOk = false;
    char* str;
    str = mpz_get_str(NULL, 2, phi);
    mpz_t testExp; mpz_init(testExp);
    mpz_t testPhi; mpz_init(testPhi);
    mpz_t helper; mpz_init(helper);


    // is GCD 1?
    while (!isOk){
        mpz_urandomb(publicExponent, rstate, strlen(str));
        mpz_set(testExp, publicExponent);
        mpz_set(testPhi, phi);
        
        // geenrate random which meets requirements
        while(mpz_cmp(phi, publicExponent) <= 0 || mpz_cmp_ui(publicExponent, 1) <= 0){
            mpz_urandomb(publicExponent, rstate, strlen(str));
            mpz_set(testExp, publicExponent);
        }
    
        // sub until get gcd
        /*
        while(mpz_cmp(testExp, testPhi) != 0)
        {
            if(mpz_cmp(testExp, testPhi) > 0){
                mpz_sub(testExp, testExp, testPhi);
            }
            else{
                mpz_sub(testPhi, testPhi, testExp);
            }
        }*/
        this->euclid(testExp, testPhi);
        
        // gcd is 1
        if (mpz_cmp_ui(testExp, 1) == 0){
            isOk = true;
        }
        
        /* mpz_gcd(helper, phi, publicExponent);
        if ((isOk && mpz_cmp_ui(helper, 1) != 0) || (!isOk && mpz_cmp_ui(helper, 1) == 0)){
            cout << "ERROR KAMO";
            exit(1);
        }*/
    }
}

/**
 * euclid to find gcd
 */
void Rsa::euclid(mpz_t u, mpz_t w){
    mpz_t r; mpz_init(r);
    while(mpz_cmp_ui(w, 0) != 0){
        mpz_mod(r, u, w);
        mpz_set(u, w);
        mpz_set(w, r);
    }
}

/**
 * function extended_gcd(a, b)
 * s := 0;    old_s := 1
 * t := 1;    old_t := 0
 * r := b;    old_r := a
 * while r ≠ 0
 *  quotient := old_r div r
 *  (old_r, r) := (r, old_r - quotient * r)
 *  (old_s, s) := (s, old_s - quotient * s)
 *  (old_t, t) := (t, old_t - quotient * t)
 * output "Bézout coefficients:", (old_s, old_t)
 * output "greatest common divisor:", old_r
 * output "quotients by the gcd:", (t, s)
 */
void Rsa::extendedEuclid(mpz_t result, mpz_t a, mpz_t b){
    mpz_t s; mpz_init_set_ui(s, 0);
    mpz_t t; mpz_init_set_ui(t, 1);
    mpz_t r; mpz_init_set(r, b);
    mpz_t old_S; mpz_init_set_ui(old_S, 1);
    mpz_t old_T; mpz_init_set_ui(old_T, 0);
    mpz_t old_R; mpz_init_set(old_R, a);
    mpz_t quotient; mpz_init(quotient);
    mpz_t helper1; mpz_init(helper1);
    mpz_t helper2; mpz_init(helper2);

    while(mpz_cmp_ui(r, 0) != 0){
        // quotient := old_r div r
        mpz_div(quotient, old_R, r);
        
        // (old_r, r) := (r, old_r - quotient * r)
        mpz_set(helper1, r);
        mpz_mul(helper2, quotient, helper1);
        mpz_sub(r, old_R, helper2);
        mpz_set(old_R, helper1);

        // (old_s, s) := (s, old_s - quotient * s)
        mpz_set(helper1, s);
        mpz_mul(helper2, quotient, helper1);
        mpz_sub(s, old_S, helper2);
        mpz_set(old_S, helper1);
        
        // (old_t, t) := (t, old_t - quotient * t)
        mpz_set(helper1, t);
        mpz_mul(helper2, quotient, helper1);
        mpz_sub(t, old_T, helper2);
        mpz_set(old_T, helper1);
    }
    
    mpz_set(result, old_S);
    
    // invert negative value to its positive representation
    // https://crypto.stackexchange.com/questions/5889/calculating-rsa-private-exponent-when-given-public-exponent-and-the-modulus-fact/5894
    if (mpz_cmp_ui(result, 0) < 0){
        mpz_add(result, result, b);
    }
}


/**
 * RSA cypher
 */
void Rsa::cypher(mpz_t cyphered, mpz_t message, mpz_t publicExponent, mpz_t N){
    mpz_powm(cyphered, message, publicExponent, N);
}

/**
 * RSA decypher
 */
void Rsa::decypher(mpz_t message, mpz_t cyphered, mpz_t privateD, mpz_t N){
    mpz_powm(message, cyphered, privateD, N);
}






