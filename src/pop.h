/* pop.h
 *
 * Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Guillaume Chapron.
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

#ifndef POP_H
#define POP_H

#include "tools.h"

#define MAX_INDIV 10000

#define NUMBER_OF_EVENTS 5
#define NUMBER_OF_STATS 15

#define IDX_POPSIZE 0
#define IDX_ALPHAS 1
#define IDX_ALPHA_PAIRS 2
#define IDX_ALPHA_GROUPS 3
#define IDX_ALPHA_SINGLE 4
#define IDX_VAGRANTS 5
#define IDX_SUBADULTS 6
#define IDX_PUPS 7
#define IDX_FEMALES 8
#define IDX_MALES 9
#define IDX_AGE 10
#define IDX_PAIRS 11
#define IDX_FAMILIES 12
#define IDX_FAMILY_SIZE 13
#define IDX_REPRODUCTIONS 14

int R_initial_pack_number;
int R_initial_vagrant_number;
int **R_initial_population;

extern double R_survival_av_PUP;
double R_survival_av_SUBADULT;
double R_survival_av_VAGRANT;
double R_survival_av_ALPHA;

double R_survival_sd_PUP;
double R_survival_sd_SUBADULT;
double R_survival_sd_VAGRANT;
double R_survival_sd_ALPHA;

double R_litter_size_av;
double R_litter_size_sd;

double R_dispersing_weib_shape_av;
double R_dispersing_weib_scale_av;
double R_settling_weib_shape_av;
double R_settling_weib_scale_av;

double R_dispersing_weib_shape_sd;
double R_dispersing_weib_scale_sd;
double R_settling_weib_shape_sd;
double R_settling_weib_scale_sd;

double R_pair1breed_av;
double R_pair1breed_sd;

int **R_quota;

typedef struct t_individual t_individual;
typedef struct t_pack t_pack;
typedef struct t_population t_population;

struct t_individual {
	int unique;
	int alive;
	int sex;
	int age;
	int stage;
	int age_disperse;
	int age_settle;
	t_pack *pack;
	t_individual *previous;
	t_individual *next;
};

struct t_pack {
	int did_bred_ever;
	int did_bred_yearbefore;
	int just_bred;
	int together;
	t_individual *alphaF;
	t_individual *alphaM;
	GPtrArray *all_members;
	t_pack *previous;
	t_pack *next;
};

struct t_population {
	int number_indiv;
	int number_initial_indiv;
	int number_indiv_history;
	int number_packs;
	double **history_indiv;

	t_individual *all_indiv;
	t_pack *all_packs;

	double survival[4];
	double litter_size;
	double dispersing_weib_shape;
	double dispersing_weib_scale;
	double settling_weib_shape;
	double settling_weib_scale;
	double pair1breed;
};

void create_population(t_population *pop);

void set_constant_parameters(t_population *pop);
void set_deterministic_parameters(t_population *pop);
void set_stochastic_parameters(t_population *pop);

void cycle_year(t_population *pop, long i, long j, struct statistics *stats);
void do_statistics(t_population *pop, long seed, long year, struct statistics *stats);
void free_population(t_population *pop);

#endif
