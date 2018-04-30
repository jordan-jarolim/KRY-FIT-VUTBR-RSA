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
    Solovay* solovay = new Solovay(keyLength);
    tRandoms randoms;
    mpz_t N;
    mpz_init(N);
    randoms = solovay->getPrimes();
    
    gmp_printf ("%s is an mpz in main %Zd\n", "here", randoms.q);
    gmp_printf ("%s is an mpz in main %Zd\n", "here", randoms.p);
    
    Rsa* rsa = new Rsa();
    rsa->getN(N, randoms.p, randoms.q);
    gmp_printf ("%s is an mpz in main %Zd\n", "here", N);
    
    char* str = mpz_get_str(NULL, 2, N);
    cout << str;

    return 0;
}
