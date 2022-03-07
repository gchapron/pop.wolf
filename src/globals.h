/* globals.h
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

#ifndef GLOBALS_H
#define GLOBALS_H

#include "tools.h"

extern long R_number_of_years;
extern long R_number_mc_runs;

extern long number_of_months;

extern struct statistics *stats;

extern int R_initial_pack_number;
extern int R_initial_vagrant_number;
extern int **R_initial_population;

extern double R_survival_av_PUP;
extern double R_survival_av_SUBADULT;
extern double R_survival_av_VAGRANT;
extern double R_survival_av_ALPHA;

extern double R_survival_sd_PUP;
extern double R_survival_sd_SUBADULT;
extern double R_survival_sd_VAGRANT;
extern double R_survival_sd_ALPHA;

extern double R_litter_size_av;
extern double R_litter_size_sd;

extern double R_dispersing_weib_shape_av;
extern double R_dispersing_weib_scale_av;
extern double R_settling_weib_shape_av;
extern double R_settling_weib_scale_av;

extern double R_dispersing_weib_shape_sd;
extern double R_dispersing_weib_scale_sd;
extern double R_settling_weib_shape_sd;
extern double R_settling_weib_scale_sd;

extern double R_pair1breed_av;
extern double R_pair1breed_sd;

extern int **R_quota;

#endif
