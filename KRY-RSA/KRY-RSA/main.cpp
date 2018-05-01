//
//  main.cpp
//  KRY-RSA
//
//  Created by Jordán Jarolím SDE on 20/04/2018.
//  Copyright © 2018 FIT VUTBR. All rights reserved.
//

#include <iostream>
#include <getopt.h>
#include <tuple>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include <gmp.h>

#include "solovay.hpp"
#include "rsa.hpp"

using namespace std;

// https://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

/**
 * Parse args
 */
int getOptions(int argc, char *argv[]) {
    int c;
    int length = 1024;
    opterr = 0;
    while ((c = getopt(argc, argv, "g:")) != -1) {
        switch (c) {
            case 'g':
                length = atoi(optarg);
                break;
            default:
                cerr << "unknown option\n";
                exit(1);
        }
    }
    return length;
}


int main(int argc, char ** argv) {
    int keyLength = getOptions(argc, argv);
    /* Init rand */
    unsigned long seed = mix(clock(), time(NULL), getpid());
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, seed);
    
    /* Init SolovayStrassen class to get primes */
    Solovay* solovay = new Solovay(keyLength);
    Rsa* rsa = new Rsa();
    tRandoms randoms;
    mpz_t N; mpz_init(N);
    mpz_t phi; mpz_init(phi);
    
    /* Generate primes, get N and check its length */
    char* str = mpz_get_str(NULL, 2, N);
    while((unsigned)strlen(str) != keyLength){
        randoms = solovay->getPrimes(rstate);
        rsa->getN(N, randoms.p, randoms.q);
        str = mpz_get_str(NULL, 2, N);
        cout << str << "\n";
        cout <<(unsigned)strlen(str)<<"\n";
        
    }
    cout << str << "\n";
    cout <<(unsigned)strlen(str)<<"\n";
    
    rsa->getPhi(phi, randoms.p, randoms.q);
    gmp_printf ("p: %Zd\n", randoms.p);
    gmp_printf ("q: %Zd\n", randoms.q);
    gmp_printf ("N: %Zd\n", N);
    gmp_printf ("Phi:  %Zd\n", phi);

    return 0;
}
