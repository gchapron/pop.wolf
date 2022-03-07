/* mc.c
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
#include <stdlib.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "mc.h"
#include "globals.h"

/////////////////////////////////////////////////////////////////////////////////
// Allocate statistics
/////////////////////////////////////////////////////////////////////////////////

void mc_allocate_statistics(struct statistics *stats) {

	stats->runs = malloc(R_number_mc_runs * sizeof(double **));

	for (long i = 0; i < R_number_mc_runs; i++) {
		stats->runs[i] = malloc(number_of_months * sizeof(double *));
		for (long j = 0; j < number_of_months; j++) {
			stats->runs[i][j] = malloc(NUMBER_OF_STATS * sizeof(double));
			for (long k = 0; k < NUMBER_OF_STATS; k++) {
					stats->runs[i][j][k] = 0.0;
			}
		}
	}

	stats->individuals = malloc(R_number_mc_runs * MAX_INDIV * sizeof(double*));

	for (int i = 0; i < R_number_mc_runs * MAX_INDIV; i++) {
		stats->individuals[i] = malloc(5 * sizeof(double));
		for (int j = 0; j < 5; j++) {
			stats->individuals[i][j] = 0;
		}
	}

	stats->number_mc_indiv_ever = 0;
	stats->number_mc_indiv_eversize = R_number_mc_runs * MAX_INDIV;

}

/////////////////////////////////////////////////////////////////////////////////
// Free results
/////////////////////////////////////////////////////////////////////////////////

void mc_free_results(struct statistics *stats) {

	for (long i = 0; i < R_number_mc_runs; i++) {
		for (long j = 0; j < number_of_months; j++) {
			free(stats->runs[i][j]);
		}
		free(stats->runs[i]);
	}
	free(stats->runs);

	for (long i = 0; i < stats->number_mc_indiv_eversize; i++) {
		free(stats->individuals[i]);
	}
	free(stats->individuals);

	free(stats);

}

/////////////////////////////////////////////////////////////////////////////////
// MONTE CARLO
/////////////////////////////////////////////////////////////////////////////////

void monte_carlo(struct statistics *stats) {

	GetRNGstate();

	long steps = R_number_mc_runs/50;

	long incr_mc_alloc = 0;

    Rprintf("\n|");

	for (long i = 0; i < R_number_mc_runs; i++) {

		t_population *pop = malloc(sizeof(t_population));

		set_constant_parameters(pop);

		set_deterministic_parameters(pop);

		create_population(pop);

        do_statistics(pop, i, 0, stats);

		for (long j = 1; j <= R_number_of_years; j++) {

//            set_stochastic_parameters(pop);

			cycle_year(pop, i, j, stats);

			if (pop->number_indiv == 0) {
				break;
			}

		}

		if ( stats->number_mc_indiv_ever + pop->number_indiv_history > stats->number_mc_indiv_eversize ) {

			incr_mc_alloc = (int)(stats->number_mc_indiv_ever + pop->number_indiv_history) / stats->number_mc_indiv_eversize + 1;

			stats->individuals = realloc(stats->individuals, incr_mc_alloc * stats->number_mc_indiv_eversize * sizeof(double *));

			for (long k = stats->number_mc_indiv_eversize; k < incr_mc_alloc * stats->number_mc_indiv_eversize; k++) {
				stats->individuals[k] = malloc(NUMBER_OF_EVENTS * sizeof(double));
				for (int l = 0; l < NUMBER_OF_EVENTS; l++) {
					stats->individuals[k][l] = 0;
				}
			}

			stats->number_mc_indiv_eversize = incr_mc_alloc * stats->number_mc_indiv_eversize;

		}

		for (int k = 0; k < pop->number_indiv_history; k++) {
			for (int l = 0; l < NUMBER_OF_EVENTS - 1; l++) {
				stats->individuals[stats->number_mc_indiv_ever + k][l] = pop->history_indiv[k][l];
			}
			stats->individuals[stats->number_mc_indiv_ever + k][NUMBER_OF_EVENTS - 1] = i;
		}

		stats->number_mc_indiv_ever += pop->number_indiv_history;

		if (steps > 0) {
			if (i % steps == 0) {
                    Rprintf("*");
			}
		}

		free_population(pop);
		free(pop);

	}

	if (steps > 0) {
            Rprintf("|");
	}

	PutRNGstate();

}
