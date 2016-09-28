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

void usage(char *prog){
    
    fprintf(stderr, "usage is '%s [options] nsamp nrep\n", prog);
    fprintf(stderr, "[options] are:\n");
    fprintf(stderr, "\t-x  #  : set seed to # (for random generator\n");
    fprintf(stderr, "\t-S xxx : read input demography from StairwayPlot output file xx\n");
    fprintf(stderr, "\t-c #.# : Kingman with Tmrca < #.#\n");
    fprintf(stderr, "\t-e #.# : Exp with rate #.#\n");
    fprintf(stderr, "\t-l #.# : Linear with foundation time #.#\n");
    fprintf(stderr, "\t-s #.# : Sudden growth at time #.#\n");
    
}


int main(int argc,  char * argv[]) {
    
    extern char *optarg;
    extern int optind;
    int o;
    
    int opt_demo='k';
    
    long seed=time(NULL);
    char *filename_SLP=NULL;
    int verbose=0;
    double tFond=0;
    
    if(argc < 3)
        usage(argv[0]),exit(1);
    
    while ( (o=getopt(argc, argv, "x:S:vc:e:s:l:")) != -1  ){
    
        switch(o){
        
            case 'x':
                seed=atol(optarg);
                break;
                
            case 'S':
                filename_SLP=optarg;
                opt_demo='S';
                break;
                
            case 'v':
                verbose=1;
                break;
                
            case 'c':
            case 'e':
            case 's':
            case 'l':
                opt_demo=(char)o;
                tFond=atof(optarg);
                break;
       }
    }
    
    srand((unsigned int) seed); //initialization of the random numbers generator
    
    int nsamp=atoi(argv[optind]);
    int nRep=atoi(argv[optind+1]);
    
    switch(opt_demo){
            
            case 'k':
                Kingman(nRep,nsamp);
                break;
            
            case 'e':
                growth_Kingman('e', tFond, nRep, nsamp);
                break;
            
            case 'c':
                conditionned_Kingman(tFond, nRep, nsamp);
                break;
            
            case 'l':
                growth_Kingman('l', tFond, nRep, nsamp);
                break;
            
            case 's':
                sudden_Kingman(tFond, nRep, nsamp);
                break;
            
            case 'S':
                stairway_demo(filename_SLP, nRep, nsamp);
                break;
            
    }

    return 0;
}
