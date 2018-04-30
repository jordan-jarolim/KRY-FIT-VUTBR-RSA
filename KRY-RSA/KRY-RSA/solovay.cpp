//
//  solovay.cpp
//  KRY-RSA
//
//  Created by Jordán Jarolím on 28.04.18.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#include "solovay.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>
#include <cstdlib>

// number of iterations for Solovay-Strassen
#define ITERS 50

using namespace std;

/**
 * Constructor
 */
Solovay::Solovay(int length):length(length){}

/**
 * Random number generator
 */
void Solovay::getRandom(mpz_t p, mpz_t q, gmp_randstate_t rstate){
    // https://crypto.stackexchange.com/questions/20474/maximum-length-in-bits-of-the-product-n-pq
    int opLength = this->length / 2;
    
    // generate p, q
    mpz_urandomb(p, rstate, opLength);
    mpz_urandomb(q, rstate, opLength);
    
    // Odds
    mpz_setbit(p, 0);
    mpz_setbit(q, 0);

    // Msb
    if (this->length >= 10){
        mpz_setbit(p, opLength-1);
        mpz_setbit(q, opLength-1);
    }
   
    // mpz_setbit(p, opLength-2);
    // mpz_setbit(q, opLength-2);
}


/**
 * Test if number is prime by Solovay-Strassen method
 */
bool Solovay::testPrime(mpz_t prime){
    int i = 0;
    int opLength = this->length / 2;
    
    // init of random generator
    unsigned long seed;
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, seed);
    
    // init random
    mpz_t veryRandomNumber; mpz_init(veryRandomNumber);
    // random a
    mpz_t a; mpz_init(a);
    // left side
    mpz_t leftCongruent; mpz_init(leftCongruent);
    // right side
    mpz_t rightCongruent; mpz_init(rightCongruent);
    // helper var
    mpz_t helper; mpz_init(helper);
    
    // prime < 2
    if (mpz_cmp_ui(prime, 2) < 0){
        return false;
    }
    
    // prime % 2
    mpz_mod_ui(helper, prime, 2);
    
    // prime != 2 || prime % 2 == 0
    if (mpz_cmp_ui(prime, 2) != 0 && mpz_cmp_ui(helper, 0) == 0){
        return false;
    }

    // repeat multiple times to increase probability
    while(i < ITERS)
    {
        i++;
        /* generate random a smaller than prime */
        mpz_urandomb(veryRandomNumber, rstate, opLength);
        mpz_sub_ui(helper, prime, 1);
        mpz_mod(a, veryRandomNumber, helper);
        mpz_add_ui(a, a, 1);
        
        /* get left side of congruence */
        mpz_t x; mpz_init_set_ui(x, 1);
        mpz_t y; mpz_init_set(y, a);
        mpz_t exponent; mpz_init(exponent);
        mpz_sub_ui(exponent, prime, 1);
        mpz_div_ui(exponent, exponent, 2);
        
        // modular exponentiation
        while(mpz_cmp_ui(exponent, 0) > 0)
        {
            mpz_mod_ui(helper, exponent, 2);
            if (mpz_cmp_ui(helper, 1) == 0){
                mpz_mul(helper, x, y);
                mpz_mod(x, helper, prime);
            }
            mpz_mul(helper, y, y);
            mpz_mod(y, helper, prime);
            mpz_div_ui(exponent, exponent, 2);
        }
        mpz_mod(leftCongruent, x, prime);
        
        /* get right side of congruence */
        // get jacob numb
        int myJac = this->jacob(a, prime);
        mpz_t jacobVar; mpz_init(jacobVar);
        if (myJac < 0){
            mpz_set_ui(jacobVar, -myJac);
            mpz_sub(rightCongruent, prime, jacobVar);
        }else{
            mpz_add_ui(rightCongruent, prime, myJac);
        }
        mpz_mod(rightCongruent, rightCongruent, prime);
        
        /* are left and right side in congruence? */
        if (mpz_cmp_ui(rightCongruent, 0) == 0 || mpz_cmp(rightCongruent, leftCongruent) != 0)
        {
            return false;
        }
    }
    return true;
}

/**
 * get Jacob symbol result
 */
int Solovay::jacob(mpz_t a_in, mpz_t n_in)
{
    int jacobValue = 1;
    mpz_t helper; mpz_init(helper);
    mpz_t helper2; mpz_init(helper2);
    // do not operate with params as they are passed as reference
    mpz_t a; mpz_init_set(a, a_in);
    mpz_t n; mpz_init_set(n, n_in);
    
    
    if (mpz_cmp_ui(a, 0) == 0){
        return 0;
    }
    if (mpz_cmp_ui(a, 0) < 0)
    {
        mpz_neg(a, a);
        mpz_mod_ui(helper, n, 4);
        if (mpz_cmp_ui(helper, 3) == 0){
            jacobValue = -jacobValue;
        }
    }
    if (mpz_cmp_ui(a, 1) == 0){
        return jacobValue;
    }
    
    while (mpz_cmp_ui(a, 0) != 0)
    {
        if (mpz_cmp_ui(a, 0) < 0)
        {
            mpz_neg(a, a);
            mpz_mod_ui(helper, n, 4);
            if (mpz_cmp_ui(helper, 3) == 0 ){
                jacobValue = -jacobValue;
            }
        }
        mpz_mod_ui(helper, a, 2);
        while (mpz_cmp_ui(helper, 0) == 0)
        {
            mpz_div_ui(a, a, 2);
            mpz_mod_ui(helper, n, 8);
            if (mpz_cmp_ui(helper, 3) == 0 || mpz_cmp_ui(helper, 5) == 0){
                jacobValue = -jacobValue;
            }
            mpz_mod_ui(helper, a, 2);
        }
        mpz_swap(a, n);
        mpz_mod_ui(helper, a, 4);
        mpz_mod_ui(helper2, n, 4);
        if (mpz_cmp_ui(helper, 3) == 0 && mpz_cmp_ui(helper2, 3) == 0){
            jacobValue = -jacobValue;
        }
        mpz_mod(a, a, n);
        mpz_div_ui(helper, n, 2);
        if (mpz_cmp(a, helper) > 0){

            mpz_sub(a, a, n);
        }
        mpz_cmp_ui(a, 0);
    }
    if (mpz_cmp_ui(n, 1) == 0){
        return jacobValue;
    }
    return 0;
}




tRandoms Solovay::getPrimes(gmp_randstate_t rstate){
    // init long vars
    mpz_t p; mpz_init(p);
    mpz_t q; mpz_init(q);
    double helper;

    
    
    
    // prime result
    bool myPrime = false;
    bool oficPrime;
    this->getRandom(p, q, rstate);

    /* Iter until you get prime p */
    while (myPrime == false){
        // gmp_printf ("%s input %Zd\n", "here", p);
        myPrime = this->testPrime(p);
        /*
        helper = pow(2, (this->length/2) - 0.5);
        cout << "helper1: "<<helper<<"\n";
        cout << "myprime: "<<myPrime<<"\n";
        if (myPrime && mpz_cmp_ui(p, floor(helper)) >= 0){
            helper = pow(2, (this->length/2));
            cout << "helper2: "<<helper<<"\n";
            if (mpz_cmp_ui(p, ceil(helper)) < 0){
                myPrime = true;
            }else{
                myPrime = false;
            }
            
        } else {
            myPrime = false;
        }
        */
        oficPrime = mpz_probab_prime_p(p, ITERS);
        // cout  << ". my: " << myPrime << ", ofic: " << oficPrime << "\n";
        if ((myPrime > 0 && oficPrime == 0) || (myPrime == 0 && oficPrime > 0)){
            cout << "ERROR KAMO";
            exit(1);
        }
        
        if (myPrime){
            // cout << "has it!!";
            // exit(0);
        }else{
            // add 2 because random numb is odd and prime number cannot be even
            mpz_add_ui(p, p, 2);
        }
        // cout << "----------------------\n";
        
    }
    
    myPrime = false;
    /* Iter until you get prime q */
    while (myPrime == false){
        // gmp_printf ("%s input %Zd\n", "here", q);
        myPrime = this->testPrime(q);
        /*
        helper = pow(2, (this->length/2) - 0.5);
        cout << "helper1: "<<helper<<"\n";
        cout << "myprime: "<<myPrime<<"\n";
        if (myPrime && mpz_cmp_ui(q, floor(helper)) >= 0){
            helper = pow(2, (this->length/2));
            cout << "helper2: "<<helper<<"\n";

            if (mpz_cmp_ui(q, ceil(helper)) < 0){
                myPrime = true;
            }else{
                myPrime = false;
            }
            
        } else {
            myPrime = false;
        }
        */
        oficPrime = mpz_probab_prime_p(q, ITERS);
        // cout  << ". my: " << myPrime << ", ofic: " << oficPrime << "\n";
        if ((myPrime > 0 && oficPrime == 0) || (myPrime == 0 && oficPrime > 0)){
            cout << "ERROR KAMO";
            exit(1);
        }
        
        if (myPrime){
            cout << "has it!!";
            // exit(0);
        }else{
            // add 2 because random numb is odd and prime number cannot be even
            mpz_add_ui(q, q, 2);
            if (mpz_cmp(p,q) == 0){
                mpz_add_ui(q, q, 2);
            }
        }
        // cout << "----------------------\n";
    }
    
    // gmp_printf ("%s is an mpz %Zd\n", "here", p);
    // gmp_printf ("%s is an mpz %Zd\n", "here", q);
    
    tRandoms randoms;
    mpz_init(randoms.p);
    mpz_init(randoms.q);
    mpz_set(randoms.p, p);
    mpz_set(randoms.q, q);
    
    return randoms;
    
};



/*
 bool Solovay::solovoyCombined(ll p, int iteration)
 {
 mpz_t a;
 mpz_init(a);
 mpz_t jacobVar; mpz_init(jacobVar);
 
 
 if (p < 2)
 return false;
 if (p != 2 && p % 2 == 0)
 return false;
 for (int i = 0; i < iteration; i++)
 {
 ll kokos = rand();
 ll aInt = kokos % (p - 1) + 1;
 ll jacobian = (p + calculateJacobian(aInt, p)) % p;
 ll mod = modulo(aInt, (p - 1) / 2, p);
 
 mpz_t leftCongruent;
 mpz_init(leftCongruent);
 // right side
 mpz_t rightCongruent;
 mpz_init(rightCongruent);
 // helper var
 mpz_t helper;
 mpz_init(helper);
 
 mpz_t a;
 mpz_init_set_ui(a, kokos);
 
 mpz_t prime;
 mpz_init_set_ui(prime, p);
 
 mpz_sub_ui(helper, prime, 1);
 gmp_printf ("%s is my (p-1) %Zd\n", "here", helper);
 cout << "true (p-1) " << (p-1) << "\n";
 
 
 mpz_mod(a, a, helper);
 gmp_printf ("%s is my rand %(p-1) %Zd\n", "here", a);
 cout << "true (p-1) " << kokos % (p-1) << "\n";
 
 mpz_add_ui(a, a, 1);
 
 
 gmp_printf ("%s is my random %Zd\n", "here", a);
 cout << "true random " << aInt << "\n";
 
 
 
 mpz_t x; mpz_init_set_ui(x, 1);
 mpz_t y; mpz_init_set(y, a);
 mpz_t exponent; mpz_init(exponent);
 mpz_sub_ui(exponent, prime, 1);
 mpz_div_ui(exponent, exponent, 2);
 
 
 while(mpz_cmp_ui(exponent, 0) > 0)
 {
 mpz_mod_ui(helper, exponent, 2);
 if (mpz_cmp_ui(helper, 1) == 0){
 mpz_mul(helper, x, y);
 mpz_mod(x, helper, prime);
 }
 mpz_mul(helper, y, y);
 mpz_mod(y, helper, prime);
 mpz_div_ui(exponent, exponent, 2);
 }
 
 mpz_mod(leftCongruent, x, prime);
 gmp_printf ("%s is my mod %Zd\n", "here", leftCongruent);
 
 
 cout << "other mod " << mod << "\n";
 fflush(stdout);
 
 
 if (mpz_cmp_ui(leftCongruent, mod) != 0){
 cout << "ERROR";
 exit(1);
 }
 
 
 int pureJac = calculateJacobian(aInt, p);
 int myJac = this->jacob(a, prime);
 cout << "pure original jacob: " << pureJac << "\n";
 cout << "my jacob: "<< myJac << "\n";
 
 if (pureJac != myJac){
 cout << "jiny jacob!!!\n";
 exit(1);
 }
 
 int extendedJacobOrig = (p + calculateJacobian(aInt, p));
 
 cout << "soucet orig: " << extendedJacobOrig << "\n";
 cout << "p v souctu: "<<p<<"\n";
 
 
 if (myJac < 0){
 mpz_set_ui(jacobVar, -myJac);
 gmp_printf ("my jac inverse: %Zd\n", jacobVar);
 
 mpz_sub(rightCongruent, prime, jacobVar);
 }else{
 cout << " >= 0\n";
 mpz_add_ui(rightCongruent, prime, myJac);
 }
 
 gmp_printf ("%s soucet my: %Zd\n", "here", rightCongruent);
 gmp_printf ("prime v souctu: %Zd\n", prime);
 
 
 if (mpz_cmp_ui(rightCongruent, extendedJacobOrig) != 0){
 cout << "++++ ERROR";
 exit(1);
 }
 
 
 mpz_mod(rightCongruent, rightCongruent, prime);
 extendedJacobOrig = extendedJacobOrig % p;
 
 cout << "extended Orig: "<<extendedJacobOrig<<"\n";
 gmp_printf ("%s extended my %Zd\n", "here", rightCongruent);
 
 if (mpz_cmp_ui(rightCongruent, extendedJacobOrig) != 0){
 cout << "\n mod ERROR\n";
 exit(1);
 }
 
 cout << p << ". je to ok\n";
 
 bool origLast = !jacobian || mod != jacobian;
 
 bool myLast = mpz_cmp_ui(rightCongruent, 0) == 0 || mpz_cmp(rightCongruent, leftCongruent) != 0;
 
 if (origLast != myLast){
 cout << "WTFFF";
 exit(1);
 }
 
 if (!jacobian || mod != jacobian)
 {
 return false;
 }
 }
 
 return true;
 }
 */



