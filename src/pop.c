/* pop.c
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

#include "pop.h"
#include "globals.h"

#define PUP 0
#define SUBADULT 1
#define VAGRANT 2
#define ALPHA 3

#define FEMALE 0
#define MALE 1

#define DISPERSED 0
#define SETTLED 1
#define FIRSTBRED 2
#define DIED 3

#define MAX_PACK_SIZE 20

#define MIN_AGE_DISPERSED 10
#define MAX_AGE_ALIVE 132
#define MAX_AGE_VAGRANTS 48
#define SEX_RATIO 0.5

t_individual* create_individual(t_population *pop, int the_sex, int the_age, int the_stage);
t_pack* create_pack_empty(t_population *pop);
t_pack* create_pack_filled(t_population *pop, int age_M, int age_F, int psize);
void individual_joins_pack(t_individual *the_indiv, int the_stage, t_pack *the_pack);
void survival_of_individuals(t_population *pop);
void remove_individuals(t_population *pop);
void free_individual(t_individual *idv);
void reproduction_of_individuals(t_population *pop);
void dispersal_of_individuals(t_population *pop);
void transition_of_individuals(t_population *pop);
void ageing_of_individuals(t_population *pop);
void remove_packs(t_population *pop);
void settle_in_packs(t_population *pop);
void hunt_individuals(t_population *pop, int *quota_month);
void cycle_month(t_population *pop);

/*******************************************************************************
 Assign population parameters
 *******************************************************************************/

void set_constant_parameters(t_population *pop) {

	pop->dispersing_weib_shape = R_dispersing_weib_shape_av;

	pop->dispersing_weib_scale = R_dispersing_weib_scale_av;

	pop->settling_weib_shape = R_settling_weib_shape_av;

	pop->settling_weib_scale = R_settling_weib_scale_av;

}

/*******************************************************************************
 Assign deterministic population parameters
 *******************************************************************************/

void set_deterministic_parameters(t_population *pop) {

	pop->survival[PUP] = pow(R_survival_av_PUP, (double)1/12);

	pop->survival[SUBADULT] = pow(R_survival_av_SUBADULT, (double)1/12);

	pop->survival[VAGRANT] = pow(R_survival_av_VAGRANT, (double)1/12);

	pop->survival[ALPHA] = pow(R_survival_av_ALPHA, (double)1/12);

	pop->litter_size = R_litter_size_av;

	pop->pair1breed = R_pair1breed_av;

	pop->dispersing_weib_shape = R_dispersing_weib_shape_av;

	pop->dispersing_weib_scale = R_dispersing_weib_scale_av;

	pop->settling_weib_shape = R_settling_weib_shape_av;

	pop->settling_weib_scale = R_settling_weib_scale_av;

}

/*******************************************************************************
 Assign stochastic population parameters
 *******************************************************************************/

void set_stochastic_parameters(t_population *pop) {

	double mu = 0;
	double sigma = 0;
	double randomv = 0;

	mu = R_survival_av_PUP;
	sigma = R_survival_sd_PUP;
	if (sigma == 0) {
		randomv = mu;
	} else {
		randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
	}
    pop->survival[PUP] = pow(randomv, (double)1/12);

	mu = R_survival_av_SUBADULT;
	sigma = R_survival_sd_SUBADULT;
	if (sigma == 0) {
		randomv = mu;
	} else {
		randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
	}
    pop->survival[SUBADULT] = pow(randomv, (double)1/12);

	mu = R_survival_av_VAGRANT;
	sigma = R_survival_sd_VAGRANT;
	if (sigma == 0) {
		randomv = mu;
	} else {
		randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
	}
    pop->survival[VAGRANT] = pow(randomv, (double)1/12);

	mu = R_survival_av_ALPHA;
	sigma = R_survival_sd_ALPHA;
	if (sigma == 0) {
		randomv = mu;
	} else {
		randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
	}
    pop->survival[ALPHA] = pow(randomv, (double)1/12);

	mu = R_litter_size_av;
	sigma = R_litter_size_sd;
	pop->litter_size = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));

	mu = R_pair1breed_av;
	sigma = R_pair1breed_sd;
	if (sigma == 0) {
		randomv = mu;
	} else {
		randomv = rbeta(beta_shape(mu, sigma), beta_rate(mu, sigma));
	}
	pop->pair1breed = randomv;

	mu = R_dispersing_weib_shape_av;
	sigma = R_dispersing_weib_shape_sd;
	pop->dispersing_weib_shape = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));

	mu = R_dispersing_weib_scale_av;
	sigma = R_dispersing_weib_scale_sd;
	pop->dispersing_weib_scale = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));

	mu = R_settling_weib_shape_av;
	sigma = R_settling_weib_shape_sd;
	pop->settling_weib_shape = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));

	mu = R_settling_weib_scale_av;
	sigma = R_settling_weib_scale_sd;
	pop->settling_weib_scale = rgamma(gamma_shape(mu, sigma), 1/gamma_rate(mu, sigma));

}

/*******************************************************************************
 Create individual in population
 *******************************************************************************/

t_individual* create_individual(t_population *pop, int the_sex, int the_age, int the_stage) {

	t_individual *new_idv = malloc(sizeof(t_individual));

	pop->number_indiv++;
	pop->number_indiv_history++;

	new_idv->unique = pop->number_indiv_history;
	new_idv->alive = 1;
	new_idv->sex = the_sex;
	new_idv->stage = the_stage;

	new_idv->age = the_age;
	new_idv->age_disperse = 0;
	new_idv->age_settle = 0;

	if (new_idv->stage == PUP) {
		new_idv->age = 7;
		new_idv->age_disperse = MIN_AGE_DISPERSED + rweibull(pop->dispersing_weib_shape, pop->dispersing_weib_scale);
		new_idv->age_settle = new_idv->age_disperse + rweibull(pop->settling_weib_shape, pop->settling_weib_scale);

	}

	if (new_idv->stage == SUBADULT) {
		new_idv->age = 19;
		new_idv->age_disperse = MIN_AGE_DISPERSED + rweibull(pop->dispersing_weib_shape, pop->dispersing_weib_scale);
		new_idv->age_disperse = fmax2(new_idv->age, new_idv->age_disperse);
		new_idv->age_settle = new_idv->age_disperse + rweibull(pop->settling_weib_shape, pop->settling_weib_scale);

	}

	if (new_idv->stage == VAGRANT) {
		new_idv->age = 19;
		new_idv->age_settle = fmax2(new_idv->age, rweibull(pop->settling_weib_shape, pop->settling_weib_scale));

	}

	if (new_idv->stage == ALPHA) {

	}

	if (pop->number_indiv == 1) {
		new_idv->previous = NULL;
		new_idv->next = NULL;
		pop->all_indiv = new_idv;
	}

	else {

		new_idv->previous = NULL;
		new_idv->next = pop->all_indiv;
		new_idv->next->previous = new_idv;
		pop->all_indiv = new_idv;

	}

	return(new_idv);

}

/*******************************************************************************
 Create empty pack in population
 *******************************************************************************/

t_pack* create_pack_empty(t_population *pop) {

	t_pack *a_pack = malloc(sizeof(t_pack));

	a_pack->all_members = g_ptr_array_sized_new(MAX_PACK_SIZE);
	a_pack->alphaF = NULL;
	a_pack->alphaM = NULL;

	pop->number_packs++;

	if (pop->number_packs == 1) {
		a_pack->previous = NULL;
		a_pack->next = NULL;
		pop->all_packs = a_pack;
	}

	else {
		a_pack->previous = NULL;
		a_pack->next = pop->all_packs;
		a_pack->next->previous = a_pack;
		pop->all_packs = a_pack;
	}

	return(a_pack);

}

/*******************************************************************************
 Create a pack with given ages of alphas and size
 *******************************************************************************/

t_pack* create_pack_filled(t_population *pop, int age_M, int age_F, int psize) {

	t_pack *a_pack = malloc(sizeof(t_pack));
	a_pack->all_members = g_ptr_array_sized_new(MAX_PACK_SIZE);

	t_individual *alpha_f = create_individual(pop, FEMALE, age_F, ALPHA);
	t_individual *alpha_m = create_individual(pop, MALE, age_M, ALPHA);

	individual_joins_pack(alpha_f, ALPHA, a_pack);
	individual_joins_pack(alpha_m, ALPHA, a_pack);

	int f = fmax2(0, psize-2);

	for (int l = 0; l < f; l++) {

		t_individual *new_idv = create_individual(pop,
												  (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE,
												  0,
												  (rbinom(1, 0.75) == 1)?PUP:SUBADULT);

		individual_joins_pack(new_idv, new_idv->stage, a_pack);

	}

	a_pack->did_bred_ever = (psize > 2)?1:0;
	a_pack->did_bred_yearbefore = a_pack->did_bred_ever;
	a_pack->just_bred = 1;
	a_pack->together = 6;

	pop->number_packs++;

	if (pop->number_packs == 1) {
		a_pack->previous = NULL;
		a_pack->next = NULL;
		pop->all_packs = a_pack;
	}

	else {
		a_pack->previous = NULL;
		a_pack->next = pop->all_packs;
		a_pack->next->previous = a_pack;
		pop->all_packs = a_pack;
	}

	return(a_pack);

}

/*******************************************************************************
 Individual joins a pack with a given stage
 *******************************************************************************/

void individual_joins_pack(t_individual *the_indiv, int the_stage, t_pack *the_pack) {

	the_indiv->pack = the_pack;
	the_indiv->stage = the_stage;

	g_ptr_array_add(the_pack->all_members, the_indiv);

	if (the_stage == ALPHA) {

		if (the_indiv->sex == FEMALE) {
			the_pack->alphaF = the_indiv;
		}

		if (the_indiv->sex == MALE) {
			the_pack->alphaM = the_indiv;
		}

		the_pack->did_bred_ever = 0;
		the_pack->did_bred_yearbefore = 0;
		the_pack->just_bred = 0;
		the_pack->together = 0;

	}

}

/*******************************************************************************
 Create Scandinavian population
 *******************************************************************************/

void create_population(t_population *pop) {

	pop->number_packs = 0;
	pop->all_packs = NULL;

	pop->number_indiv = 0;
	pop->number_indiv_history = 0;
	pop->all_indiv = NULL;

	pop->history_indiv = malloc(MAX_INDIV * sizeof(double*));
	for (int i = 0; i < MAX_INDIV; i++) {
		pop->history_indiv[i] = malloc((NUMBER_OF_EVENTS - 1)*sizeof(double));
		for (int j = 0; j < (NUMBER_OF_EVENTS - 1); j++) {
			pop->history_indiv[i][j] = 0;
		}
	}

	for (int i = 0; i < R_initial_pack_number; i++) {

		create_pack_filled(pop,
						   12*R_initial_population[i][0]+7,
						   12*R_initial_population[i][1]+7,
						   R_initial_population[i][2]);

	}

	int trans = R_initial_vagrant_number;

	for (int l = 0; l < trans; l++) {

		create_individual(pop,
						  (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE,
						  0,
						  VAGRANT);

	}

	pop->number_initial_indiv = pop->number_indiv;

}

/*******************************************************************************
 Survival of individuals
 *******************************************************************************/

void survival_of_individuals(t_population *pop) {

	t_individual *current_idv = pop->all_indiv;

	double the_surv;

	while (current_idv != NULL) {

		the_surv = pop->survival[current_idv->stage];

		if ((current_idv->stage == PUP) & (current_idv->age >= 4)) {

			the_surv = pop->survival[SUBADULT];

		}

		current_idv->alive = rbinom(1, the_surv);

		if ((current_idv->stage == VAGRANT) & (current_idv->age >= MAX_AGE_VAGRANTS)) {
			current_idv->alive = 0;
		}

		if (current_idv->age >= MAX_AGE_ALIVE) {
			current_idv->alive = 0;
		}

		current_idv = current_idv->next;

	}

}

/*******************************************************************************
 Removal of individuals
 *******************************************************************************/

void remove_individuals(t_population *pop) {

	t_individual *current_idv = pop->all_indiv;
	t_individual *next_idv = NULL;

	while (current_idv != NULL) {

		if (current_idv->alive == 0) {

			pop->history_indiv[current_idv->unique-1][DIED] = current_idv->age;

			if (current_idv->stage == ALPHA) {

				if (current_idv->sex == FEMALE) {
					current_idv->pack->alphaF = NULL;
				}

				if (current_idv->sex == MALE) {
					current_idv->pack->alphaM = NULL;
				}

				g_ptr_array_remove_fast(current_idv->pack->all_members, current_idv);
				current_idv->pack = NULL;
			}

			if (current_idv->stage == PUP || current_idv->stage == SUBADULT) {

				g_ptr_array_remove_fast(current_idv->pack->all_members, current_idv);
				current_idv->pack = NULL;

			}

			next_idv = current_idv->next;

			if ( (current_idv->previous == NULL) & (current_idv->next == NULL) ) {

				pop->all_indiv = NULL;

			}

			else if ( (current_idv->previous == NULL) & (current_idv->next != NULL) ) {

				current_idv->next->previous = NULL;
				pop->all_indiv = current_idv->next;

			}

			else if ( (current_idv->previous != NULL) & (current_idv->next == NULL) ) {

				current_idv->previous->next = NULL;

			}

			else {

				current_idv->next->previous = current_idv->previous;
				current_idv->previous->next = current_idv->next;

			}

			free_individual(current_idv);

			pop->number_indiv--;

			current_idv = next_idv;

		} else {

			current_idv = current_idv->next;

		}

	}

}

/*******************************************************************************
 Free an individual
 *******************************************************************************/

void free_individual(t_individual *idv) {

	free(idv);

}

/*******************************************************************************
 Reproduction of individuals
 *******************************************************************************/

void reproduction_of_individuals(t_population *pop) {

	t_pack *current_pack = pop->all_packs;

	int f;
	int will_breed = 0;

	while (current_pack != NULL) {

		if (current_pack->did_bred_ever == 1) {
			will_breed = 1;
		}

		if (current_pack->just_bred == 1) {
			current_pack->did_bred_yearbefore = 1;
		}

		if (current_pack->did_bred_ever == 0) {
			will_breed = rbinom(1, pop->pair1breed);
		}

		if (current_pack->together <= 3) {
			will_breed = 0;
		}

		if ((current_pack->alphaF != NULL)
			& (current_pack->alphaM != NULL)
			& (will_breed == 1) ) {

			f = fmax2(0, pop->litter_size);
			f = rpois(f);

			if (f > 0) {

				if (current_pack->did_bred_ever == 0) {

					if (pop->history_indiv[current_pack->alphaF->unique-1][FIRSTBRED] == 0) {
						pop->history_indiv[current_pack->alphaF->unique-1][FIRSTBRED] = current_pack->alphaF->age;
					}

					if (pop->history_indiv[current_pack->alphaM->unique-1][FIRSTBRED] == 0) {
						pop->history_indiv[current_pack->alphaM->unique-1][FIRSTBRED] = current_pack->alphaM->age;
					}

				}

				current_pack->did_bred_ever = 1;
				current_pack->just_bred = 1;

				for (int l = 0; l < f; l++) {

					t_individual *new_idv = malloc(sizeof(t_individual));
					new_idv->alive = 1;
					new_idv->sex = (rbinom(1, SEX_RATIO) == 1)?MALE:FEMALE;
					new_idv->age = 0;
					new_idv->stage = PUP;
					new_idv->pack = current_pack;

					new_idv->age_disperse = MIN_AGE_DISPERSED + rweibull(pop->dispersing_weib_shape, pop->dispersing_weib_scale);
					new_idv->age_settle = new_idv->age_disperse + 1 + rweibull(pop->settling_weib_shape, pop->settling_weib_scale);

					new_idv->previous = NULL;
					new_idv->next = pop->all_indiv;
					new_idv->next->previous = new_idv;
					pop->all_indiv = new_idv;

					pop->number_indiv++;
					pop->number_indiv_history++;

					g_ptr_array_add(current_pack->all_members, new_idv);

					new_idv->unique = pop->number_indiv_history;

					if (pop->number_indiv_history/MAX_INDIV == (double)pop->number_indiv_history/MAX_INDIV) {
						pop->history_indiv = realloc(pop->history_indiv, (pop->number_indiv_history + MAX_INDIV) * sizeof(double *));

						for (int i = 0; i < MAX_INDIV; i++) {
							pop->history_indiv[pop->number_indiv_history+i] = malloc((NUMBER_OF_EVENTS - 1)*sizeof(double));
							for (int j = 0; j < (NUMBER_OF_EVENTS - 1); j++) {
								pop->history_indiv[pop->number_indiv_history+i][j] = 0;
							}
						}
					}

				}

			}

		}

		current_pack = current_pack->next;

	}

}

/*******************************************************************************
 Dispersal of individuals
 *******************************************************************************/

void dispersal_of_individuals(t_population *pop) {

	t_individual *current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		if (current_idv->stage == PUP || current_idv->stage == SUBADULT) {

			if ( current_idv->age == current_idv->age_disperse) {

				current_idv->stage = VAGRANT;
				pop->history_indiv[current_idv->unique-1][DISPERSED] = current_idv->age;

				g_ptr_array_remove_fast(current_idv->pack->all_members, current_idv);

				current_idv->pack = NULL;

			}

		}

		current_idv = current_idv->next;

	}

}

/*******************************************************************************
 Transition of individuals
 *******************************************************************************/

void transition_of_individuals(t_population *pop) {

	t_individual *current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		if ( (current_idv->age == 12) & (current_idv->stage == PUP) ) {

            current_idv->stage = SUBADULT;

        }

		current_idv = current_idv->next;

	}

}

/*******************************************************************************
 Ageing of individuals
 *******************************************************************************/

void ageing_of_individuals(t_population *pop) {

    t_individual *current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		current_idv->age++;

		current_idv = current_idv->next;

	}

	t_pack *current_pack = pop->all_packs;

	while (current_pack != NULL) {

		if ((current_pack->alphaF != NULL)
			& (current_pack->alphaM != NULL)) {

			current_pack->together++;

		}

		current_pack = current_pack->next;

	}

}

/*******************************************************************************
 Update number of packs
 *******************************************************************************/

void remove_packs(t_population *pop) {

	t_pack *current_pack = pop->all_packs;
	t_pack *next_pack = NULL;

	t_individual *current_idv = NULL;

	while (current_pack != NULL) {

		if (current_pack->alphaF == NULL && current_pack->alphaM == NULL) {

			for (int i = 0; i < current_pack->all_members->len; i++) {

				current_idv = (t_individual*)g_ptr_array_index(current_pack->all_members, i);
				current_idv->pack = NULL;

				if (current_idv->age < 7) {
					current_idv->alive = 0;
				}

				current_idv->stage = VAGRANT;

			}

			g_ptr_array_free(current_pack->all_members);

			next_pack = current_pack->next;

			if ( (current_pack->previous == NULL) & (current_pack->next == NULL) ) {
				pop->all_packs = NULL;
			}

			else if ( (current_pack->previous == NULL) & (current_pack->next != NULL) ) {
				current_pack->next->previous = NULL;
				pop->all_packs = current_pack->next;
			}

			else if ( (current_pack->previous != NULL) & (current_pack->next == NULL) ) {
				current_pack->previous->next = NULL;
			}

			else {
				current_pack->next->previous = current_pack->previous;
				current_pack->previous->next = current_pack->next;
			}

			free(current_pack);

			pop->number_packs--;

			current_pack = next_pack;

		} else {

			current_pack = current_pack->next;

		}

	}

	remove_individuals(pop);

}


/*******************************************************************************
 Settle in packs
 *******************************************************************************/

void settle_in_packs(t_population *pop) {

	t_individual *mate_M = NULL;
	t_individual *mate_F = NULL;
	t_individual *current_idv = pop->all_indiv;

	GPtrArray *array_vagrant_males = g_ptr_array_sized_new(pop->number_indiv);
	GPtrArray *array_vagrant_females = g_ptr_array_sized_new(pop->number_indiv);

	while (current_idv != NULL) {

		if ((current_idv->stage == VAGRANT)
			& (current_idv->age >= current_idv->age_settle)) {

			if (current_idv->sex == MALE) {
				g_ptr_array_add(array_vagrant_males, current_idv);
			}

			if (current_idv->sex == FEMALE) {
				g_ptr_array_add(array_vagrant_females, current_idv);
			}

		}

		current_idv = current_idv->next;

	}

	GPtrArray *array_packs_no_male = g_ptr_array_sized_new(pop->number_packs);
	GPtrArray *array_packs_no_female = g_ptr_array_sized_new(pop->number_packs);

	t_pack *current_pack = pop->all_packs;

	while (current_pack != NULL) {

		if (current_pack->alphaM == NULL) {
			g_ptr_array_add(array_packs_no_male, current_pack);
		}

		if (current_pack->alphaF == NULL) {
			g_ptr_array_add(array_packs_no_female, current_pack);
		}

		current_pack = current_pack->next;

	}

	if ((array_packs_no_male->len > 0) & (array_vagrant_males->len > 0)) {

		int replacement_males = fmin2(array_vagrant_males->len, array_packs_no_male->len);

		g_ptr_array_shuffle(array_packs_no_male);
		g_ptr_array_shuffle(array_vagrant_males);

		for (int i = 0; i < replacement_males; i++) {

			current_pack = (t_pack*)g_ptr_array_index(array_packs_no_male, i);
			current_idv = (t_individual*)g_ptr_array_index(array_vagrant_males, i);

			individual_joins_pack(current_idv, ALPHA, current_pack);

			pop->history_indiv[current_idv->unique-1][SETTLED] = current_idv->age;

		}

	}

	g_ptr_array_free(array_packs_no_male);

	if ((array_packs_no_female->len > 0) & (array_vagrant_females->len > 0)) {

		int replacement_females = fmin2(array_vagrant_females->len, array_packs_no_female->len);

		g_ptr_array_shuffle(array_packs_no_female);
		g_ptr_array_shuffle(array_vagrant_females);

		for (int i = 0; i < replacement_females; i++) {

			current_pack = (t_pack*)g_ptr_array_index(array_packs_no_female, i);
			current_idv = (t_individual*)g_ptr_array_index(array_vagrant_females, i);

			individual_joins_pack(current_idv, ALPHA, current_pack);

			pop->history_indiv[current_idv->unique-1][SETTLED] = current_idv->age;

		}

	}

	g_ptr_array_free(array_packs_no_female);

    g_ptr_array_empty(array_vagrant_males);
    g_ptr_array_empty(array_vagrant_females);

	current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		if ((current_idv->stage == VAGRANT)
			& (current_idv->age >= current_idv->age_settle) ) {

			if (current_idv->sex == MALE) {
				g_ptr_array_add(array_vagrant_males, current_idv);
			}

			if (current_idv->sex == FEMALE) {
				g_ptr_array_add(array_vagrant_females, current_idv);
			}

		}

		current_idv = current_idv->next;

	}

	if ((array_vagrant_males->len > 0) & (array_vagrant_females->len > 0)) {

		int pairing_trials = fmin2(array_vagrant_males->len, array_vagrant_females->len);

		mate_M = NULL;
		mate_F = NULL;

		g_ptr_array_shuffle(array_vagrant_males);
		g_ptr_array_shuffle(array_vagrant_females);

		for (int i = 0; i < pairing_trials; i++) {

			mate_M = (t_individual*)g_ptr_array_index(array_vagrant_males, i);
			mate_F = (t_individual*)g_ptr_array_index(array_vagrant_females, i);

			t_pack *a_pack = create_pack_empty(pop);
			individual_joins_pack(mate_F, ALPHA, a_pack);
			individual_joins_pack(mate_M, ALPHA, a_pack);

			pop->history_indiv[mate_F->unique-1][SETTLED] = mate_F->age;
			pop->history_indiv[mate_M->unique-1][SETTLED] = mate_M->age;

		}

	}

    g_ptr_array_empty(array_vagrant_females);

	current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		if ((current_idv->stage == VAGRANT)
			& (current_idv->age >= current_idv->age_settle)) {

			if (current_idv->sex == FEMALE) {
				g_ptr_array_add(array_vagrant_females, current_idv);
			}

		}

		current_idv = current_idv->next;

	}

    for (int i = 0; i < array_vagrant_females->len; i++) {

        t_pack *a_pack = create_pack_empty(pop);

		mate_F = (t_individual*)g_ptr_array_index(array_vagrant_females, i);

		individual_joins_pack(mate_F, ALPHA, a_pack);

		pop->history_indiv[mate_F->unique-1][SETTLED] = mate_F->age;

    }

	g_ptr_array_free(array_vagrant_males);
	g_ptr_array_free(array_vagrant_females);

}

/*******************************************************************************
 Hunting individuals
 *******************************************************************************/

void hunt_individuals(t_population *pop, int *quota_month) {

	GPtrArray *array_hunted_individuals = g_ptr_array_sized_new(pop->number_indiv);
	GPtrArray *array_hunted_packs = g_ptr_array_sized_new(pop->number_packs);

    t_individual *current_idv;
	t_pack *current_pack;

    int pos;
	int quota_effective;

	if (quota_month[0] > 0) {

		current_pack = pop->all_packs;

		while (current_pack != NULL) {

			if ( (current_pack->alphaF != NULL)
				& (current_pack->alphaM != NULL) ) {
				g_ptr_array_add(array_hunted_packs, current_pack);
			}
			current_pack = current_pack->next;

		}

		quota_effective = fmin2(quota_month[0], array_hunted_packs->len);

		while (quota_effective > 0) {

			pos = (int)runif(0, array_hunted_packs->len-1);
			current_pack = (t_pack*)g_ptr_array_index(array_hunted_packs, pos);

			current_pack->alphaF->alive = 0;
			current_pack->alphaM->alive = 0;

			g_ptr_array_remove_index_fast(array_hunted_packs, pos);
			quota_effective--;

		}

	}

	if (quota_month[1] > 0) {

		current_idv = pop->all_indiv;

		while (current_idv != NULL) {

			if ( (current_idv->stage == ALPHA)
                & (current_idv->alive == 1)) {
				g_ptr_array_add(array_hunted_individuals, current_idv);
			}
			current_idv = current_idv->next;

		}

		quota_effective = fmin2(quota_month[1], array_hunted_individuals->len);

		while (quota_effective > 0) {

			pos = (int)runif(0, array_hunted_individuals->len-1);
			current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
			current_idv->alive = 0;

			g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
			quota_effective--;

		}

	}

    g_ptr_array_empty(array_hunted_individuals);

	if (quota_month[2] > 0) {

		current_idv = pop->all_indiv;

		while (current_idv != NULL) {

			if ( (current_idv->stage == VAGRANT)
                & (current_idv->alive == 1)) {
				g_ptr_array_add(array_hunted_individuals, current_idv);
			}
			current_idv = current_idv->next;

		}

		quota_effective = fmin2(quota_month[2], array_hunted_individuals->len);

		while (quota_effective > 0) {

			pos = (int)runif(0, array_hunted_individuals->len-1);
			current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
			current_idv->alive = 0;

			g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
			quota_effective--;

		}

	}

    g_ptr_array_empty(array_hunted_individuals);

	if (quota_month[3] > 0) {

		current_idv = pop->all_indiv;
		while (current_idv != NULL) {

			if ( (current_idv->age > 5)
                & (current_idv->stage != ALPHA)
                & (current_idv->stage != VAGRANT)
                & (current_idv->alive == 1)) {
				g_ptr_array_add(array_hunted_individuals, current_idv);
			}

			current_idv = current_idv->next;

		}

		quota_effective = fmin2(quota_month[3], array_hunted_individuals->len);

		while (quota_effective > 0) {

			pos = (int)runif(0, array_hunted_individuals->len-1);
			current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);
			current_idv->alive = 0;

			g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
			quota_effective--;

		}

	}

	g_ptr_array_empty(array_hunted_individuals);

	if (quota_month[4] > 0) {

		current_idv = pop->all_indiv;

		while (current_idv != NULL) {

            if ( (current_idv->age > 5)
                & (current_idv->alive == 1) ) {

                g_ptr_array_add(array_hunted_individuals, current_idv);

            }

            current_idv = current_idv->next;

		}

		quota_effective = fmin2(quota_month[4], array_hunted_individuals->len);

		while (quota_effective > 0) {

			pos = (int)runif(0, array_hunted_individuals->len-1);
			current_idv = (t_individual*)g_ptr_array_index(array_hunted_individuals, pos);

			if (current_idv->alive == 1) {

				current_idv->alive = 0;

				g_ptr_array_remove_index_fast(array_hunted_individuals, pos);
				quota_effective--;
			}

		}

	}

	g_ptr_array_free(array_hunted_individuals);
	g_ptr_array_free(array_hunted_packs);

	remove_individuals(pop);
	remove_packs(pop);

}

/*******************************************************************************
 Monthly life cycle
 *******************************************************************************/


void cycle_month(t_population *pop) {

	survival_of_individuals(pop);

	remove_individuals(pop);

	dispersal_of_individuals(pop);

	remove_packs(pop);

	settle_in_packs(pop);

	ageing_of_individuals(pop);

    transition_of_individuals(pop);

}

/*******************************************************************************
 Yearly life cycle
 *******************************************************************************/

void cycle_year(t_population *pop, long i, long j, struct statistics *stats) {

	long k = 12*(j-1) + 1;

	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
	hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	reproduction_of_individuals(pop);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

	k++;
	cycle_month(pop);
    hunt_individuals(pop, R_quota[k-1]);
	do_statistics(pop, i, k, stats);

}

/*******************************************************************************
 Collect statistics
 *******************************************************************************/

void do_statistics(t_population *pop, long the_run, long the_month, struct statistics *stats) {

	int var = 0;
    double *stat = malloc(NUMBER_OF_STATS*sizeof(double));
	for (int i = 0; i < NUMBER_OF_STATS; i++) {
		stat[i] = 0.0;
	}

    t_individual* current_idv;
	t_pack *current_pack;

	stat[IDX_POPSIZE] = pop->number_indiv;

	current_idv = pop->all_indiv;

	while (current_idv != NULL) {

		if (current_idv->stage == ALPHA) {
            stat[IDX_ALPHAS]++;
        }

        if (current_idv->stage == ALPHA) {
			if((current_idv->pack->alphaF != NULL)
			   & (current_idv->pack->alphaM != NULL)
			   & (current_idv->pack->all_members->len == 2)) {
					stat[IDX_ALPHA_PAIRS]++;
			}
		}

		if (current_idv->stage == ALPHA) {
			if (current_idv->pack->all_members->len > 2) {
					stat[IDX_ALPHA_GROUPS]++;
				}
		}

		if (current_idv->stage == ALPHA) {
			if (current_idv->pack->all_members->len == 1) {
				stat[IDX_ALPHA_SINGLE]++;
			}
		}

		if (current_idv->stage == VAGRANT) {
            stat[IDX_VAGRANTS]++;
        }

		if (current_idv->stage == SUBADULT) {
            stat[IDX_SUBADULTS]++;
        }

		if (current_idv->stage == PUP) {
            stat[IDX_PUPS]++;
        }

		if (current_idv->sex == FEMALE) {
            stat[IDX_FEMALES]++;
        }

		if (current_idv->sex == MALE) {
            stat[IDX_MALES]++;
        }

		stat[IDX_AGE] += current_idv->age;

		current_idv = current_idv->next;
	}

	current_pack = pop->all_packs;

	while (current_pack != NULL) {

		if ( ((current_pack->alphaM != NULL) & (current_pack->alphaF != NULL))
            & (current_pack->all_members->len == 2) ) {
            stat[IDX_PAIRS]++;
        }

		if (current_pack->all_members->len > 2) {
			if (current_pack->alphaM != NULL || current_pack->alphaF != NULL) {
				stat[IDX_FAMILIES]++;
			}
        }

		if (current_pack->all_members->len > 2) {
			if (current_pack->alphaM != NULL || current_pack->alphaF != NULL) {
				stat[IDX_FAMILY_SIZE] += current_pack->all_members->len;
			}
		}

		var = 0;
		for (int i = 0; i < current_pack->all_members->len; i++) {
			current_idv = g_ptr_array_index(current_pack->all_members, i);
			if ( current_idv->stage == PUP ) {
				var = 1;
				break;
			}
		}
		stat[IDX_REPRODUCTIONS] = stat[IDX_REPRODUCTIONS] + var;

		current_pack = current_pack->next;
	}

	for (int i = 0; i < NUMBER_OF_STATS; i++) {
		stats->runs[the_run][the_month][i] = stat[i];
	}

	if (pop->number_indiv > 0) {
		stats->runs[the_run][the_month][IDX_AGE] = stats->runs[the_run][the_month][IDX_AGE] / pop->number_indiv / 12;
	}

	if (stats->runs[the_run][the_month][IDX_FAMILIES] > 0) {
		stats->runs[the_run][the_month][IDX_FAMILY_SIZE] = stats->runs[the_run][the_month][IDX_FAMILY_SIZE] / stats->runs[the_run][the_month][IDX_FAMILIES];

	} else {
        stats->runs[the_run][the_month][IDX_FAMILY_SIZE] = 0;
    }

	free(stat);

}

/*******************************************************************************
 Free population
 *******************************************************************************/

void free_population(t_population *pop) {

	t_individual *next_idv;

	while (pop->all_indiv != NULL) {
		next_idv = pop->all_indiv->next;
		free_individual(pop->all_indiv);
		pop->all_indiv = next_idv;
	}

	int k = ((int)pop->number_indiv_history/MAX_INDIV + 1)*MAX_INDIV;
	for (long i = 0; i < k; i++) {
		free(pop->history_indiv[i]);
	}
	free(pop->history_indiv);

	t_pack *next_pack;

	while (pop->all_packs != NULL) {
		next_pack = pop->all_packs->next;
		g_ptr_array_free(pop->all_packs->all_members);
		free(pop->all_packs);
		pop->all_packs = next_pack;
	}

}
