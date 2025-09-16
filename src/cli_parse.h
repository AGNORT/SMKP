#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct {
    const char *prog_name;
	const char *input_file;		//the input file
	const char *output_file;	//the output file
	const char *method;         //the method to solve the problem, default is DP or BPC, specifed as Gurobi to use the gurobi solver
} Args;

int parse_args(Args* args, int argc, char** argv);

#ifdef __cplusplus
}
#endif

#endif
