/* tools.h
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

#ifndef TOOLS_H
#define TOOLS_H

struct statistics {
	long number_mc_indiv_ever;
	long number_mc_indiv_eversize;
	long counter;
	double ***runs;
	double **individuals;
	double **posteriors;
};

double beta_shape(double mu, double sigma);
double beta_rate(double mu, double sigma);
double gamma_shape(double mu, double sigma);
double gamma_rate(double mu, double sigma);

#define g_ptr_array_index(array, index_) ((array)->pdata)[index_]

typedef struct _GPtrArray GPtrArray;
typedef void* gpointer;

struct _GPtrArray {
	gpointer *pdata;
	int len;
	int alloc;
};

GPtrArray* g_ptr_array_sized_new(int reserved_size);
void g_ptr_array_add (GPtrArray *array, gpointer data);
void g_ptr_array_remove_index_fast(GPtrArray *array, int index);
void g_ptr_array_remove_fast(GPtrArray *array, gpointer data);
void g_ptr_array_empty(GPtrArray *array);
void g_ptr_array_free(GPtrArray *array);
void g_ptr_array_shuffle(GPtrArray *array);

#endif
