/*
! Copyright (C) 2015-2024 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST. It is created based on Base64Transcoder.h in Qbox
!
! Contributors to this file:
! Huihuo Zheng
*/
#include <complex.h>
typedef unsigned char byte;

char etable[64];  // encode table
byte dtable[256]; // decode table

void b64init();

int encode(int nbytes, const byte* const from, char* const to);

int encode_complex(double complex* from, int n, char* const to) {
  // n -- dimension of complex array
  int nbytes=sizeof(double complex)*n;
  return encode(nbytes, (byte*) from, to);
}

int encode_double(double* from, int n, char* const to) {
  // n -- dimension of array
  int nbytes=sizeof(double)*n;
  return encode(nbytes, (byte*) from, to);
}

int decode(int nc, const char* const from, byte* const to );

void byteswap_double(int nbytes, double* const x);

void byteswap_complex(int nbytes, double complex* const x) {
  byteswap_double(nbytes, (double*) x);
};

int nchars(int nbytes) { return 4*((nbytes + 2) / 3 ); }

int nbytes(int nc) { return 3*nc/4; }

int decode_double(const char* const from, int n, double* to ) {
  int nb = sizeof(double)*n;
  int nc=nchars(nb);
  return decode(nc, from, (byte*) to);
};

int decode_complex(const char* const from, int n, double complex* to ) {
  int nc=nchars(sizeof(double complex)*n);
  return decode(nc, from, (byte*) to);
};
