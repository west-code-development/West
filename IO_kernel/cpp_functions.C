/*
 Copyright (C) 2015-2016 M. Govoni 
 This file is distributed under the terms of the
 GNU General Public License. See the file `License'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt .

 This file is part of WEST.

 Contributors to this file: 
 He Ma
 */

#include <string>
#include <sys/stat.h>
#include <unistd.h>

using namespace std; 

extern "C" {
  void c_sleep(int* seconds) {
    usleep(useconds_t(*seconds));
  }
  
  void c_wait_for_file(int* max_seconds, bool* file_exists, const char* lockfilename) {
    struct stat statbuf;
    int status;
    *file_exists = false;
    for ( int i = 0; i < *max_seconds; i++ )
    {
      status = stat(lockfilename, &statbuf);
      if (status == 0)
      {
        *file_exists = true;
        break;
      }
      else
        usleep(1000000);
    }
  }
}
