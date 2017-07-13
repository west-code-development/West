/*
! Copyright (C) 2015-2017 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST. It is created based on Base64Transcoder.C in Qbox
!
! Contributors to this file: 
! Huihuo Zheng
*/

#include "cbase64.h"

void b64init() {
  int i; 
  /// this is to initialized the encode and decode tables
  for (i = 0; i < 26; i++) {
      etable[i] = 'A' + i;
      etable[26 + i] = 'a' + i;
  }
  for (i = 0; i < 10; i++) {
      etable[52 + i] = '0' + i;
  }
  etable[62] = '+';
  etable[63] = '/';
  
  for (i = 0; i < 255; i++) {
      dtable[i] = 0x80;
  }
  for (i = 'A'; i <= 'Z'; i++) {
      dtable[i] = 0 + (i - 'A');
  }
  for (i = 'a'; i <= 'z'; i++) {
      dtable[i] = 26 + (i - 'a');
  }
  for (i = '0'; i <= '9'; i++) {
      dtable[i] = 52 + (i - '0');
  }
  dtable['+'] = 62;
  dtable['/'] = 63;
  dtable['='] = 0;
}

int encode(int nbytes, const byte* const from, char* const to) {
  const byte* fptr = from;
  char* tptr = to;

  int n3 = nbytes / 3; // number of groups of three bytes

  while ( n3-- > 0 )
  {
    byte ig0 = *fptr++;
    byte ig1 = *fptr++;
    byte ig2 = *fptr++;

    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = etable[((ig1 & 0xF) << 2) | (ig2 >> 6)];
    *tptr++ = etable[ig2 & 0x3F];
  }

  int nr = nbytes % 3; // remaining bytes

  if ( nr == 2 )
  {
    byte ig0 = *fptr++;
    byte ig1 = *fptr++;
    byte ig2 = 0;

    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = etable[((ig1 & 0xF) << 2) | (ig2 >> 6)];
    *tptr++ = '=';
  }
  else if ( nr == 1 )
  {
    byte ig0 = *fptr++;
    byte ig1 = 0;

    *tptr++ = etable[ig0 >> 2];
    *tptr++ = etable[((ig0 & 3) << 4) | (ig1 >> 4)];
    *tptr++ = '=';
    *tptr++ = '=';
  }
  return 0;
}

int decode(int nchars, const char* const from, byte* const to )
{
    // Decode Base64 chars in array "from" into bytes in array "to"
  // White space and new lines are skipped
  // extra characters at end that do not form a valid group of 4 chars are
  // ignored.
  // nchars: number of chars in array "from"
  // the number of bytes successfully translated is returned

  byte a0,a1,a2,a3,b0,b1,b2,b3;
  int c;
  const char* fptr = from;
  const char* const fptr_end = from+nchars+1;
  byte* tptr = to;

  while ( fptr < fptr_end-4 )
  {
    // get 4 valid characters from input string
    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c0 "
           << endl;
#endif
      break;
    }
    a0 = (byte) c;
    b0 = (byte) dtable[c];

    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c1 "
           << endl;
#endif
      break;
    }
    a1 = (byte) c;
    b1 = (byte) dtable[c];

    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c2 "
           << endl;
#endif
      break;
    }
    a2 = (byte) c;
    b2 = (byte) dtable[c];

    do
    {
      c = *fptr++;
    }
    while ( (c <= ' ') && (fptr < fptr_end) );
    if ( (c <= ' ') && fptr >= fptr_end )
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: end of string reached reading c3\n"
           << " (without reading a valid c3) " << endl;
#endif
      break;
    }
    a3 = (byte) c;
    b3 = (byte) dtable[c];

    if ((b0|b1|b2|b3) & 0x80)
    {
#ifdef DEBUG
      cerr << " Base64Transcoder::decode: Illegal character in input: "
           << endl;
#endif
      return tptr - to;
    }

    if ( a3 == '=' )
    {
      if ( a2 == '=' )
      {
        // write 1 byte only
        *tptr++ = (b0 << 2) | (b1 >> 4);
      }
      else
      {
        // write 2 bytes only
        *tptr++ = (b0 << 2) | (b1 >> 4);
        *tptr++ = (b1 << 4) | (b2 >> 2);
      }
    }
    else
    {
      // write 3 bytes
      *tptr++ = (b0 << 2) | (b1 >> 4);
      *tptr++ = (b1 << 4) | (b2 >> 2);
      *tptr++ = (b2 << 6) | b3;
    }

  }
#ifdef DEBUG
  if ( fptr >= fptr_end )
  {
    cerr << " Base64Transcoder::decode: end of string reached in input: "
         << endl;
  }
#endif

  return tptr - to;
}

void byteswap_double(int n, double* const x) {
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 ) {
      unsigned char tmp;
      tmp = c[7]; c[7] = c[0]; c[0] = tmp;
      tmp = c[6]; c[6] = c[1]; c[1] = tmp;
      tmp = c[5]; c[5] = c[2]; c[2] = tmp;
      tmp = c[4]; c[4] = c[3]; c[3] = tmp;
    c+=8;
  }
}
