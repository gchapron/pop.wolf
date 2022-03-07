/* main.c
 *
 * Copyright (C) 2011-2022 Guillaume Chapron.
 * guillaume.chapron@slu.se
 * with contributions from Camilla Wikenros, Olof Liberg, Øystein Flagstad,
 * Cyril Milleret, Johan Månsson, Linn Svensson, Barbara Zimmermann,
 * Mikael Åkesson, Petter Wabakken & Håkan Sand
 *
 * This file is part of 'pop.wolf', a R package to simulate wolf populations
 *
 * 'pop.wolf' is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * 'pop.wolf' is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 'pop.wolf'. If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "mc.h"
#include "tools.h"
#include "pop.h"

#define extern
# include "globals.h"
#undef extern

/*******************************************************************************
 Simple Monte Carlo
 *******************************************************************************/

SEXP C_montecarlo (
				   SEXP SEXP_years,
				   SEXP SEXP_runs,
				   SEXP SEXP_pp_surviving,
				   SEXP SEXP_sb_surviving,
				   SEXP SEXP_vg_surviving,
				   SEXP SEXP_ad_surviving,
				   SEXP SEXP_dispersing_weib_shape,
				   SEXP SEXP_dispersing_weib_scale,
				   SEXP SEXP_settling_weib_shape,
				   SEXP SEXP_settling_weib_scale,
				   SEXP SEXP_pair1breed,
				   SEXP SEXP_litter_size,
				   SEXP SEXP_quota,
				   SEXP SEXP_initial_packs,
				   SEXP SEXP_initial_vagrants) {

	int nprot = 0;

	R_number_of_years = INTEGER(SEXP_years)[0];
	R_number_mc_runs = INTEGER(SEXP_runs)[0];

	number_of_months = 12*R_number_of_years;
	number_of_months++;

    R_survival_av_PUP = REAL(SEXP_pp_surviving)[0];
	R_survival_av_SUBADULT = REAL(SEXP_sb_surviving)[0];
	R_survival_av_VAGRANT = REAL(SEXP_vg_surviving)[0];
	R_survival_av_ALPHA = REAL(SEXP_ad_surviving)[0];

	R_survival_sd_PUP = REAL(SEXP_pp_surviving)[1];
	R_survival_sd_SUBADULT = REAL(SEXP_sb_surviving)[1];
	R_survival_sd_VAGRANT = REAL(SEXP_vg_surviving)[1];
	R_survival_sd_ALPHA = REAL(SEXP_ad_surviving)[1];

	R_litter_size_av = REAL(SEXP_litter_size)[0];
	R_litter_size_sd = REAL(SEXP_litter_size)[1];

	R_dispersing_weib_shape_av = REAL(SEXP_dispersing_weib_shape)[0];
	R_dispersing_weib_scale_av = REAL(SEXP_dispersing_weib_scale)[0];
	R_settling_weib_shape_av = REAL(SEXP_settling_weib_shape)[0];
	R_settling_weib_scale_av = REAL(SEXP_settling_weib_scale)[0];

	R_dispersing_weib_shape_sd = REAL(SEXP_dispersing_weib_shape)[1];
	R_dispersing_weib_scale_sd = REAL(SEXP_dispersing_weib_scale)[1];
	R_settling_weib_shape_sd = REAL(SEXP_settling_weib_shape)[1];
	R_settling_weib_scale_sd = REAL(SEXP_settling_weib_scale)[1];

	R_pair1breed_av = REAL(SEXP_pair1breed)[0];
	R_pair1breed_sd = REAL(SEXP_pair1breed)[1];

	R_initial_vagrant_number = INTEGER(SEXP_initial_vagrants)[0];

	R_quota = malloc(LENGTH(SEXP_quota)/5 * sizeof(int*));
	for (int i = 0; i < LENGTH(SEXP_quota)/5; i++) {
		R_quota[i] = malloc(5 * sizeof(int));
		for (int j = 0; j < 5; j++) {
			R_quota[i][j] = INTEGER(SEXP_quota)[i + j * LENGTH(SEXP_quota)/5];
		}
	}

	R_initial_population = malloc(LENGTH(SEXP_initial_packs)/3 * sizeof(int*));
	for (int i = 0; i < LENGTH(SEXP_initial_packs)/3; i++) {
		R_initial_population[i] = malloc(3 * sizeof(int));
		for (int j = 0; j < 3; j++) {
			R_initial_population[i][j] = INTEGER(SEXP_initial_packs)[i + j * LENGTH(SEXP_initial_packs)/3];
		}
	}

	R_initial_pack_number = LENGTH(SEXP_initial_packs)/3;

	stats = malloc(sizeof(struct statistics));

	mc_allocate_statistics(stats);

	monte_carlo(stats);

	SEXP R_runs;

	PROTECT(R_runs = allocVector(REALSXP, number_of_months * R_number_mc_runs * NUMBER_OF_STATS)); nprot++;

	for (long i = 0; i < R_number_mc_runs; i++) {
		for (long j = 0; j < number_of_months; j++) {
			for (long k = 0; k < NUMBER_OF_STATS; k++) {
				REAL(R_runs)[ i*number_of_months*NUMBER_OF_STATS + j*NUMBER_OF_STATS + k] = stats->runs[i][j][k];
			}
		}
	}

	SEXP R_individuals;

	PROTECT(R_individuals = allocVector(REALSXP, NUMBER_OF_EVENTS * stats->number_mc_indiv_ever)); nprot++;
	for (long i = 0; i < stats->number_mc_indiv_ever; i++) {
		for (long j = 0; j < NUMBER_OF_EVENTS; j++) {
			REAL(R_individuals)[ i * NUMBER_OF_EVENTS + j] = stats->individuals[i][j];
		}
	}

	char *names[2] = {"runs", "individuals"};
	SEXP list, list_names;

	PROTECT(list_names = allocVector(STRSXP,2)); nprot++;
	for (int i = 0; i < 2; i++) {
		SET_STRING_ELT(list_names, i , mkChar(names[i]));
	}

	PROTECT(list = allocVector(VECSXP, 2)); nprot++;

	SET_VECTOR_ELT(list, 0, R_runs);

	SET_VECTOR_ELT(list, 1, R_individuals);

	setAttrib(list, R_NamesSymbol, list_names);

	UNPROTECT(nprot);

	mc_free_results(stats);

	return(list);

}


