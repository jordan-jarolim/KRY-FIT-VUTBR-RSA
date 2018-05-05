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

#define BREAK_ITERATIONS 1000000

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
    mpz_clear(pMinus);
    mpz_clear(qMinus);
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
    }
    mpz_clear(testExp);
    mpz_clear(testPhi);
    mpz_clear(helper);
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
    mpz_clear(r);
}

/**
 * modular pow
 */
void Rsa::modPow(mpz_t result, mpz_t base, mpz_t exponent, mpz_t mod){
    mpz_t x; mpz_init_set_ui(x, 1);
    mpz_t y; mpz_init_set(y, base);
    mpz_t helper; mpz_init(helper);
    mpz_t expo; mpz_init_set(expo, exponent);
    
    while(mpz_cmp_ui(expo, 0) > 0)
    {
        mpz_mod_ui(helper, expo, 2);
        if (mpz_cmp_ui(helper, 1) == 0){
            mpz_mul(helper, x, y);
            mpz_mod(x, helper, mod);
        }
        mpz_mul(helper, y, y);
        mpz_mod(y, helper, mod);
        mpz_div_ui(expo, expo, 2);
    }
    mpz_mod(result, x, mod);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(helper);
    mpz_clear(expo);
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
    mpz_clear(s);
    mpz_clear(r);
    mpz_clear(t);
    mpz_clear(old_S);
    mpz_clear(old_R);
    mpz_clear(old_T);
    mpz_clear(quotient);
    mpz_clear(helper1);
    mpz_clear(helper2);
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

/**
 * RSA breaker
 */
void Rsa::breakIt(mpz_t p, mpz_t N, gmp_randstate_t rstate){
    bool found = false;
    int counter = 0;
    mpz_t result; mpz_init(result);
    mpz_set_ui(p, 2);
    
    // is even?
    mpz_mod(result, N, p);
    if (mpz_cmp_ui(result, 0) == 0){
        return;
    }
    
    // set odd
    mpz_set_ui(p, 1);
    
    /* Try brute force */
    while (!found && counter < BREAK_ITERATIONS){
        counter++;

        mpz_add_ui(p, p, 2);
        mpz_mod(result, N, p);
        if (mpz_cmp_ui(result, 0) == 0){
            found = true;
        }
    }
    
    mpz_t x; mpz_init(x);
    mpz_t y; mpz_init(y);
    mpz_t c; mpz_init(c);
    mpz_t d; mpz_init(d);
    mpz_t helper; mpz_init(helper);
    mpz_t helper2; mpz_init(helper2);

    /* Pollard's Rho */
    while (!found){
        /* get random x 2 - N */
        mpz_sub_ui(helper, N, 2);
        mpz_urandomm(x, rstate, helper);
        mpz_add_ui(x, x, 2);
        mpz_set(y, x);

        /* get random c 1 - N */
        mpz_sub_ui(helper, N, 1);
        mpz_urandomm(c, rstate, helper);
        mpz_add_ui(c, c, 1);

        /* init candidate */
        mpz_set_ui(d, 1);
        
        /* GCD is 1 */
        while (mpz_cmp_ui(d, 1) == 0)
        {
            /* Tortoise */
            mpz_set_ui(helper, 2);
            modPow(helper2, x, helper, N);
            mpz_add(x, helper2, c);
            mpz_add(x, x, N);
            mpz_mod(x, x, N);

            /* Hare Move */
            modPow(helper2, y, helper, N);
            mpz_add(y, helper2, c);
            mpz_add(y, y, N);
            mpz_mod(y, y, N);
            
            /* Hare Move */
            modPow(helper2, y, helper, N);
            mpz_add(y, helper2, c);
            mpz_add(y, y, N);
            mpz_mod(y, y, N);

            /* abs x-y */
            mpz_sub(helper, x, y);
            if (mpz_cmp_ui(helper, 0) < 0){
                mpz_neg(helper, helper);
            }
            
            /* get gcd */
            mpz_set(helper2, N);
            this->euclid(helper, helper2);

            /* set value */
            mpz_set(d, helper);
            
            if (mpz_cmp(d, N) != 0){
                found = true;
            }
        }
        mpz_set(p, d);
    }
    
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(c);
    mpz_clear(d);
    mpz_clear(helper);
    mpz_clear(helper2);
}






