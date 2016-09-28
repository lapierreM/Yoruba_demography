//
//  growthSFS.c
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

/*
 * Function: read_StairwayPlot
 * ----------------------------
 *      Reads the output summary file of the stairway plot method to get times (in years) and population sizes (Ne)
 *
 *      filename: output summary file of the stairway plot method
 *      Time: empty array to fill with times in years
 *      Size: empty array to fill with population sizes
 *
 *      returns: int, sample size of the stairway plot simulations
 */
int read_StairwayPlot(char * filename, double **Time, double **Size){
    
    
    FILE *f;
    int i, u=0;
    int samplesize;

    f = fopen(filename, "r");
    if(f==NULL)fprintf(stderr, "Cannot open file %s. Please check the file, bye\n", filename),exit(1);
    
    
    while( fgetc(f) != '\t' );
    fscanf(f,"%d",&samplesize);
    
    
    
    (*Time)=(double *)malloc( (size_t) sizeof(double)*(samplesize-1) );
    (*Size)=(double *)malloc( (size_t) sizeof(double)*(samplesize-1) );
    if( *Time ==NULL || *Size==NULL)
        fprintf(stderr, "read_SkyLinePlot: cannot malloc arrays, sorry, bye\n"),exit(3);
    
    
    for(i=0;i<6;i++)
        while( fgetc(f) != '\n' );
    
    
    for(u=0;u<(samplesize-1);u++)
    {        
        while( fgetc(f) != '\n' );

        for(i=0;i<5;i++)
            while( fgetc(f) != '\t' );
        
        fscanf(f,"%lf%lf",(*Time)+u,(*Size)+u);

        while( fgetc(f) != '\n' );
        
    }
    
    fclose(f);
    
    return samplesize;
    
}

/*
 * Function: stairway_demo
 * ----------------------------
 *      Computes the waiting times under a demography predicted by the stairway plot method
 *
 *      filename_SLP: output summary file of the stairway plot method
 *      nRep: number of repetitions for computing the mean waiting times
 *      nsamp: sample size
 *
 *      returns: void
 */
void stairway_demo(char * filename_SwP, int nRep, int nsamp)
{
    
    double *Time, *Size;

    int samplesize;
    
    samplesize=read_StairwayPlot(filename_SwP, &Time, &Size);
    
    
    int i=0;
    int j=0;
    int k=0;

    int i_max=samplesize-1;
    
    /*
    //scale times and sizes by a factor 0.1 because the genome length was divided by 10 for the stairway plot inference
    
    for (i=0 ; i<i_max ; i++)
    {
        Time[i]*=0.1;
        Size[i]*=0.1;
    }
    
    //convert times from years to coalescent units
    
    for (i=0 ; i<i_max ; i++)
    {
        Time[i]= ( Time[i] / 24.0 ) / ( Size[0] );
        //printf("Time %d : %f\n",i,Time[i]);
    }
    
    //normalize sizes so that N(0)=1
    
    double size_zero=Size[0];
    for (i=0 ; i<i_max ; i++)
    {
        Size[i] = Size[i]/size_zero;
        //printf("Size %d : %f\n",i,Size[i]);
    }
    */
    
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k_star[i]=0;
    
    double* t_k=NULL;
    t_k=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in the stairway demography
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k[i]=0;
    
    double tmp=0;

    double tfin=0;
    
    for (j=0 ; j<nRep ; j++)
    {
        i=0;
        tfin=0;
        
        for (k=nsamp ; k>1 ; k--)
        {
            int r=0;
            while (r==0)
                r=rand();
            
            t_k_star[k-2]=-( 2*log( r/(double)RAND_MAX ) )/(k*(k-1.0));
            //printf("t_k_star %d : %f\n",k,t_k_star[k-2]);
        }
        
        for (k=nsamp ; k>1 ; k--)
        {
            tmp = t_k_star[k-2] * Size[i];
            
            while (i<(i_max-1) && ( Time[i] - tfin ) < tmp )
            {
                tmp -= ( Time[i] - tfin );
                
                t_k[k-2] += ( Time[i] - tfin );
                
                tmp *= ( Size[i+1] / Size[i] );
                
                tfin = Time[i];
                
                i++;
            }
            
            t_k[k-2] += tmp;
            
            tfin += tmp;
            
        }
    }
    
    free(t_k_star);
    free(Time);
    free(Size);

    for (k=nsamp ; k>1 ; k--)
    {
        t_k[k-2]/=nRep;
        //printf("t_k %d :%f\n",k,t_k[k-2]);
    }
    
    double tmrca=0;
    
    for (k=2 ; k<=nsamp ; k++)
    {
        tmrca+=t_k[k-2];
    }
    
    //printf("tmrca : %f\n",tmrca);

    
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(nsamp-1)*sizeof(double*));
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(nsamp-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,nsamp);
    
    sim_SFS(t_k,P_i_k,nsamp);
    
    free(t_k);
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        free(P_i_k[k]);
    }
    free(P_i_k);
    
}

/*
 * Function: growth_Kingman
 * ----------------------------
 *      Computes the waiting times under linear or exponential growth
 *
 *      croissance: the type of growth ('e' for exponential and 'l' for linear)
 *      tFond: the time parameter to characterize growth
 *      nRep: number of repetitions for computing the mean waiting times
 *      nsamp: sample size
 *      
 *      returns: void
 */
void growth_Kingman(char croissance, double tFond, int nRep, int nsamp)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k_star[i]=0;

    double* t_k=NULL;
    t_k=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in exponential or linear coalescent
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k[i]=0;
    
    double tmp=0;
    double sum_v=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        sum_v=0;
        
        for (k=nsamp ; k>1 ; k--)
        {
            int r=0;
            while (r==0)
                r=rand();
            
            t_k_star[k-2]=-( 2*log( r/(double)RAND_MAX ) )/(k*(k-1.0));
        }
        
        for (k=nsamp ; k>1 ; k--)
        {
            switch (croissance)
            {
                case 'e':
                    tmp=(1/tFond)*log(1+(tFond*t_k_star[k-2]*exp(-tFond*sum_v)));
                    t_k[k-2]+=tmp;
                    break;
                case 'l':
                    tmp=(1-exp(-t_k_star[k-2]/tFond))*(tFond-sum_v);
                    t_k[k-2]+=tmp;
                    break;
            }
            sum_v+=tmp;
        }
    }
    free(t_k_star);
    
    for (k=2 ; k<=nsamp ; k++)
    {
        t_k[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=nsamp ; k++)
    {
        length_tree+=k*t_k[k-2];
        tmrca+=t_k[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    
    double length_tree_without_singletons=0;
    
    //for k=2
    length_tree_without_singletons+=2*t_k[0]*(1-((double)2/(double)(nsamp-1)));
                                                
    for (k=3 ; k<=nsamp ; k++)
    {
        length_tree_without_singletons+=k*t_k[k-2]*(1-((double)(k-1)/(double)(nsamp-1)));
    }
    
    printf("Length tree without singletons : %f\n",length_tree_without_singletons);
    */
    
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(nsamp-1)*sizeof(double*));
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(nsamp-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,nsamp);
    
    sim_SFS(t_k,P_i_k,nsamp);
    
    free(t_k);
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        free(P_i_k[k]);
    }
    free(P_i_k);
}

/*
 * Function: conditionned_Kingman
 * ----------------------------
 *      Computes the waiting times under the conditioned model (Kingman that conditions all lineage to have coalesce before 't')
 *
 *      tFond: the time parameter to characterize the conditionment
 *      nRep: number of repetitions for computing the mean waiting times
 *      nsamp: sample size
 *
 *      returns: void
 */
void conditionned_Kingman(double tFond, int nRep, int nsamp)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent for simulations <tFond
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k_star[i]=0;
    
    double* tmp_t_k_star=NULL;
    tmp_t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent for current simulation
    for (i=0 ; i<(nsamp-1) ; i++)
        tmp_t_k_star[i]=0;
    
    double sum_t_k=0;
    double t=0;
    
    i=0;
    while (i<nRep)
    {
        for (k=2 ; k<=nsamp && sum_t_k<tFond ; k++)
        {
            t=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
            tmp_t_k_star[k-2]=t;
            sum_t_k+=t;
        }
        
        if (sum_t_k<tFond) // if we keep the current simulation
        {
            for (k=2 ; k<=nsamp ; k++)
            {
                t_k_star[k-2]+=tmp_t_k_star[k-2];
            }
            i++;
        }
        sum_t_k=0;
    }
    free(tmp_t_k_star);

    for (k=2 ; k<=nsamp ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=nsamp ; k++)
    {
        length_tree+=k*t_k_star[k-2];
        tmrca+=t_k_star[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    */
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(nsamp-1)*sizeof(double*));
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(nsamp-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,nsamp);
    
    sim_SFS(t_k_star,P_i_k,nsamp);
    
    free(t_k_star);
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        free(P_i_k[k]);
    }
    free(P_i_k);
}

/*
 * Function: Kingman
 * ----------------------------
 *      Computes the waiting times under the standard Kingman coalescent
 *
 *      nRep: number of repetitions for computing the mean waiting times
 *      nsamp: sample size
 *
 *      returns: void
 */
void Kingman(int nRep, int nsamp)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k_star[i]=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        for (k=nsamp ; k>1 ; k--)
        {
            t_k_star[k-2]+=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
        }
    }
    
    for (k=2 ; k<=nsamp ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(nsamp-1)*sizeof(double*));
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(nsamp-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,nsamp);
    
    /* affichage des P_i_k
    for (i=2 ; i<=nsamp ; i++)
    {
        for (k=1 ; k<=(nsamp-i+1) ; k++)
        {
            printf("%f\t",P_i_k[k-1][i-2]);
        }
        printf("\n");
    }
    */
    
    sim_SFS(t_k_star,P_i_k,nsamp);
    
    free(t_k_star);
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        free(P_i_k[k]);
    }
    free(P_i_k);
}

/*
 * Function: sudden_Kingman
 * ----------------------------
 *      Computes the waiting times under the Sudden model (Kingman that started at time t (it was 1 individual before that))
 *
 *      tFond: the time parameter to characterize the time of sudden growth
 *      nRep: number of repetitions for computing the mean waiting times
 *      nsamp: sample size
 *
 *      returns: void
 */
void sudden_Kingman(double tFond, int nRep, int nsamp)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(nsamp-1)*sizeof(double)); //waiting times in standard coalescent, O if sum>tFond
    for (i=0 ; i<(nsamp-1) ; i++)
        t_k_star[i]=0;
    
    double sum_t_k=0;
    double t=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        for (k=nsamp ; k>1 && sum_t_k<tFond ; k--)
        {
            t=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
            sum_t_k+=t;
            
            if (sum_t_k<tFond)
            {
                t_k_star[k-2]+=t;
            }
            else
            {
                t_k_star[k-2]+=tFond-sum_t_k+t;
                break;
            }
        }
        sum_t_k=0;
    }
    
    for (k=2 ; k<=nsamp ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=nsamp ; k++)
    {
        length_tree+=k*t_k_star[k-2];
        tmrca+=t_k_star[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    */
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(nsamp-1)*sizeof(double*));
    
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(nsamp-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,nsamp);
    
    sim_SFS(t_k_star,P_i_k,nsamp);
    
    free(t_k_star);
    for (k=0 ; k<(nsamp-1) ; k++)
    {
        free(P_i_k[k]);
    }
    free(P_i_k);
}

/*
 * Function: fill_P_i_k
 * ----------------------------
 *      Fills the P_i_k matrix (probability that a random line at state i gives k descendants at state n (see Fu 1995)
 *
 *      P_i_k: matrix to be filled
 *      nsamp: sample size
 *
 *      returns: void
 */
void fill_P_i_k(double **P_i_k, int nsamp)
{
    int k=0;
    int i=0;
    
    for (k=1 ; k<nsamp ; k++)
    {
        P_i_k[k-1][0]=1/(double)(nsamp-1);
        for (i=2 ; i<(nsamp-k+1) ; i++)
        {
            P_i_k[k-1][i-1]=P_i_k[k-1][i-2]*((double)i*(double)(nsamp-k-i+1))/((double)(i-1)*(double)(nsamp-i));
        }
    }
}

/*
 * Function: sim_SFS
 * ----------------------------
 *      Computes the SFS with given waiting times and P_i_k matrix
 *
 *      t_k: list of waiting times
 *      P_i_k: matrix of probabilities P_i_k
 *      nsamp: sample size
 *
 *      returns: void
 */
void sim_SFS(double *t_k, double **P_i_k, int nsamp)
{
    int i=0;
    int k=0;

    double* SFS=NULL;
    SFS=(double*)malloc((size_t)(nsamp-1)*sizeof(double));
    for (i=0 ; i<(nsamp-1) ; i++)
        SFS[i]=0;
    
    //double sumSFS=0;
    
    for (k=1 ; k<nsamp ; k++)
    {
        for (i=2 ; i<(nsamp-k+2) ; i++)
        {
            SFS[k-1]+=(double)i*t_k[i-2]*P_i_k[k-1][i-2];
        }
        //sumSFS+=SFS[k-1];
    }
    
    for (k=1 ; k<nsamp ; k++)
    {
        //SFS[k-1]/=sumSFS;
        printf("%f\n",SFS[k-1]);
    }
    
    free(SFS);
}

