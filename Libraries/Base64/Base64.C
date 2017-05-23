/*
 Copyright (C) 2015-2016 M. Govoni 
 This file is distributed under the terms of the
 GNU General Public License. See the file `License'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt .

 This file is part of WEST.

 Contributors to this file: 
 Huihuo Zheng
 */

#include <iostream>
#include "Base64Transcoder.h"
#include <complex>
#include "Base64.h"

using namespace std; 

void Base64::encode(double *a, int n, char *b) {
  int nb = sizeof(double)*n; 
  bt.encode(nb, (unsigned char*) &a[0], b);
}

void Base64::encode(complex<double> *a, int n, char *b) {
  int nb = sizeof(complex<double>)*n; 
  bt.encode(nb, (unsigned char*) &a[0], b);
}


void Base64::decode(char *b, int n, double *a) {
  int nb = bt.nchars(sizeof(double)*n); 
  bt.decode(nb, b, (unsigned char*) &a[0]); 
}

void Base64::decode(char *b, int n, complex<double> *a) {
  int nb = bt.nchars(sizeof(complex<double>)*n); 
  bt.decode(nb, b, (unsigned char*) &a[0]); 
}
void Base64::byteswap_double(int n, double* const x) {
  bt.byteswap_double(n, x); 
}
void Base64::byteswap_complex(int n, complex<double>* const x) {
  bt.byteswap_double(n*2, (double*) x); 
}


extern "C" {
  #include "Base64.h"
  Base64 *Base64__new() {
    return new Base64(); 
  }
  
  void Base64__encode_double(Base64 *This, double *a, int n, char* to) {
    This->encode(a, n, to); 
  }

  void Base64__encode_complex(Base64 *This, complex<double> *a, int n, char* to) {
    This->encode(a, n, to); 
  }

  void Base64__decode_double(Base64 *This, char* from, int n, double *a) {
    This->decode(from, n, a); 
  }

  void Base64__decode_complex(Base64 *This, char* from, int n, complex<double> *a) {
    This->decode(from, n, a); 
  }
  void Base64__byteswap_double(Base64 *This, int n, double *x) {
    This->byteswap_double(n, x); 
  }
  void Base64__byteswap_complex(Base64 *This, int n, complex<double> *x) {
    This->byteswap_complex(n, x); 
  }
}
