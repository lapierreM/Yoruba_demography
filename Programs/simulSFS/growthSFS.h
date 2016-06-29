//
//  growthSFS.h
//  simulSFS
//
/*
 Copyright (C) 2016  Lapierre/Marguerite
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public License
 as published by the Free Software Foundation; either version 2.1
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 for more information, please contact: marguerite.lapierre@mnhn.fr
 */

#ifndef growthSFS_h
#define growthSFS_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void Kingman(int nRep, int taillePop);
void conditionned_Kingman(double tFond, int nRep, int taillePop);
void sudden_Kingman(double tFond, int nRep, int taillePop);
void growth_Kingman(char croissance, double tFond, int nRep, int taillePop);
void sim_SFS(double *t_k, double **P_i_k, int taillePop);
void fill_P_i_k(double **P_i_k, int taillePop);

#endif /* growthSFS_h */
