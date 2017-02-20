/* mc.h
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

#ifndef MC_H
#define MC_H

#include "pop.h"
#include "tools.h"

long R_number_of_years;
long R_number_mc_runs;

long number_of_months;

void mc_allocate_statistics(struct statistics *stats);
void mc_free_results(struct statistics *stats);
void monte_carlo(struct statistics *stats);

#endif
