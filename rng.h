//Copyright (C) 2015 Thomas Brown, Xavier Didelot, Daniel J. Wilson, Nicola De Maio
//
// This file is part of SimBac.
//
// SimBac is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SimBac is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SimBac.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __RNG_H__
#define __RNG_H__

#include <gsl/gsl_rng.h>
#include <time.h>
#include <getopt.h>

gsl_rng * rng;

void makerng() {
    const gsl_rng_type *rng_type;
    long int rng_seed;
    gsl_rng_env_setup();
    rng_type = gsl_rng_default;
    rng_seed = gsl_rng_default_seed;
    rng = gsl_rng_alloc (rng_type);
    unsigned int seed;
    FILE *devrandom;

 if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
   seed = (unsigned long) time(NULL);
   printf("Got seed %u from time()\n",seed);
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   printf("Got seed %u from /dev/urandom\n",seed);
   fclose(devrandom);
 }
    //seed=0;//This is for debugging purposes only
    gsl_rng_set(rng,seed);
}


#endif // __RNG_H__
