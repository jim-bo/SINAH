/*
 * SINAH.h
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

#ifndef SINAH_H_
#define SINAH_H_

// system
#include <fstream>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// stl
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <utility>

// logging and debugging.
#define DEBUG 1
#define VERBOSE 1
char * DEBUGTXT = new char[5012];

static void WRITE_OUT(const char * txt){
	if( VERBOSE == 1 ) {
		fprintf(stdout,txt);
		fflush(stdout);
	}
}

static void WRITE_ERR(const char * txt){
	fprintf(stderr,txt);
	fflush(stdout);
}

static bool ACTIVITY(int i, int j, int k){
	if(i != 0 && (i % k) == 0){
		if( VERBOSE == 1 ) fprintf(stdout,"%d of %d\n", i, j);
		fflush(stdout);
		if( DEBUG == 1 ) {
			return true;
		}
	}
	return false;
}

#endif /* SINAH_H_ */
