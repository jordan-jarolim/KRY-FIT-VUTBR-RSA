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

#include <gmp.h>

#include "solovay.hpp"
#include "rsa.hpp"





using namespace std;

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
    // init rand
    unsigned long seed;
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, seed);
    
    Solovay* solovay = new Solovay(keyLength);
    tRandoms randoms;
    mpz_t N;
    mpz_init(N);
    
    char* str = mpz_get_str(NULL, 2, N);

    while((unsigned)strlen(str) != keyLength){
        randoms = solovay->getPrimes(rstate);
        
        // gmp_printf ("%s is an mpz in main %Zd\n", "here", randoms.p);
        // gmp_printf ("%s is an mpz in main %Zd\n", "here", randoms.q);
        
        Rsa* rsa = new Rsa();
        rsa->getN(N, randoms.p, randoms.q);
        // gmp_printf ("%s is an mpz in main %Zd\n", "here", N);
        
        str = mpz_get_str(NULL, 2, N);
        
    }
    cout << str << "\n";
    cout <<(unsigned)strlen(str)<<"\n";
    

    return 0;
}
