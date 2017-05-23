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
#ifndef BASE64__H__
#define BASE64__H__
#include <iostream>
#include "Base64Transcoder.h"
#include <complex>
using namespace std; 
class Base64 {
  Base64Transcoder bt;
 public:
  void encode(double *a, int n, char *b); 
  void decode(char *b, int n, double *); 
  void encode(complex<double> *a, int n, char *b); 
  void decode(char *b, int n, complex<double> *); 
  void byteswap_double(int n, double* const); 
  void byteswap_complex(int n, complex<double>* const); 
}; 
#endif
