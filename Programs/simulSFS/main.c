//
//  main.c
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

#include "growthSFS.h"


int main(int argc, const char * argv[]) {
    
    srand((unsigned int) time(NULL)); //initialization of the random numbers generator
    
    char croissance;
    croissance=argv[1][0];
    int taillePop;
    taillePop=atoi(argv[2]);
    double tFond;
    tFond=atof(argv[3]);
    int nRep;
    nRep=atoi(argv[4]);
    
    if (croissance=='k')
        Kingman(nRep,taillePop);
    else if (croissance=='c')
        conditionned_Kingman(tFond, nRep, taillePop);
    else if (croissance=='s')
        sudden_Kingman(tFond, nRep, taillePop);
    else
        growth_Kingman(croissance, tFond, nRep, taillePop);
    return 0;
}
