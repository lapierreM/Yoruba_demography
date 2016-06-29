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
 * Function: growth_Kingman
 * ----------------------------
 *      Computes the waiting times under linear or exponential growth
 *
 *      croissance: the type of growth ('e' for exponential and 'l' for linear)
 *      tFond: the time parameter to characterize growth
 *      nRep: number of repetitions for computing the mean waiting times
 *      taillePop: sample size
 *      
 *      returns: void
 */

void growth_Kingman(char croissance, double tFond, int nRep, int taillePop)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in standard coalescent
    for (i=0 ; i<(taillePop-1) ; i++)
        t_k_star[i]=0;

    double* t_k=NULL;
    t_k=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in exponential or linear coalescent
    for (i=0 ; i<(taillePop-1) ; i++)
        t_k[i]=0;
    
    double tmp=0;
    double sum_v=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        sum_v=0;
        
        for (k=taillePop ; k>1 ; k--)
        {
            t_k_star[k-2]=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
        }
        
        for (k=taillePop ; k>1 ; k--)
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
    
    for (k=2 ; k<=taillePop ; k++)
    {
        t_k[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=taillePop ; k++)
    {
        length_tree+=k*t_k[k-2];
        tmrca+=t_k[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    
    double length_tree_without_singletons=0;
    
    //for k=2
    length_tree_without_singletons+=2*t_k[0]*(1-((double)2/(double)(taillePop-1)));
                                                
    for (k=3 ; k<=taillePop ; k++)
    {
        length_tree_without_singletons+=k*t_k[k-2]*(1-((double)(k-1)/(double)(taillePop-1)));
    }
    
    printf("Length tree without singletons : %f\n",length_tree_without_singletons);
    */
    
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(taillePop-1)*sizeof(double*));
    
    for (k=0 ; k<(taillePop-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(taillePop-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,taillePop);
    
    sim_SFS(t_k,P_i_k,taillePop);
    
    free(t_k);
    for (k=0 ; k<(taillePop-1) ; k++)
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
 *      taillePop: sample size
 *
 *      returns: void
 */
void conditionned_Kingman(double tFond, int nRep, int taillePop)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in standard coalescent for simulations <tFond
    for (i=0 ; i<(taillePop-1) ; i++)
        t_k_star[i]=0;
    
    double* tmp_t_k_star=NULL;
    tmp_t_k_star=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in standard coalescent for current simulation
    for (i=0 ; i<(taillePop-1) ; i++)
        tmp_t_k_star[i]=0;
    
    double sum_t_k=0;
    double t=0;
    
    i=0;
    while (i<nRep)
    {
        for (k=2 ; k<=taillePop && sum_t_k<tFond ; k++)
        {
            t=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
            tmp_t_k_star[k-2]=t;
            sum_t_k+=t;
        }
        
        if (sum_t_k<tFond) // if we keep the current simulation
        {
            for (k=2 ; k<=taillePop ; k++)
            {
                t_k_star[k-2]+=tmp_t_k_star[k-2];
            }
            i++;
        }
        sum_t_k=0;
    }
    free(tmp_t_k_star);

    for (k=2 ; k<=taillePop ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=taillePop ; k++)
    {
        length_tree+=k*t_k_star[k-2];
        tmrca+=t_k_star[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    */
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(taillePop-1)*sizeof(double*));
    
    for (k=0 ; k<(taillePop-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(taillePop-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,taillePop);
    
    sim_SFS(t_k_star,P_i_k,taillePop);
    
    free(t_k_star);
    for (k=0 ; k<(taillePop-1) ; k++)
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
 *      taillePop: sample size
 *
 *      returns: void
 */
void Kingman(int nRep, int taillePop)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in standard coalescent
    for (i=0 ; i<(taillePop-1) ; i++)
        t_k_star[i]=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        for (k=taillePop ; k>1 ; k--)
        {
            t_k_star[k-2]+=-(2*log((double)(rand()+1)/(double)((unsigned)RAND_MAX + 1)))/(k*(k-1));
        }
    }
    
    for (k=2 ; k<=taillePop ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(taillePop-1)*sizeof(double*));
    
    for (k=0 ; k<(taillePop-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(taillePop-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,taillePop);
    
    /* affichage des P_i_k
    for (i=2 ; i<=taillePop ; i++)
    {
        for (k=1 ; k<=(taillePop-i+1) ; k++)
        {
            printf("%f\t",P_i_k[k-1][i-2]);
        }
        printf("\n");
    }
    */
    
    sim_SFS(t_k_star,P_i_k,taillePop);
    
    free(t_k_star);
    for (k=0 ; k<(taillePop-1) ; k++)
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
 *      taillePop: sample size
 *
 *      returns: void
 */
void sudden_Kingman(double tFond, int nRep, int taillePop)
{
    int i=0;
    int k=0;
    
    double* t_k_star=NULL;
    t_k_star=(double *)malloc((size_t)(taillePop-1)*sizeof(double)); //waiting times in standard coalescent, O if sum>tFond
    for (i=0 ; i<(taillePop-1) ; i++)
        t_k_star[i]=0;
    
    double sum_t_k=0;
    double t=0;
    
    for (i=0 ; i<nRep ; i++)
    {
        for (k=taillePop ; k>1 && sum_t_k<tFond ; k--)
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
    
    for (k=2 ; k<=taillePop ; k++)
    {
        t_k_star[k-2]/=nRep;
    }
    /*
    double length_tree=0;
    double tmrca=0;
    
    for (k=2 ; k<=taillePop ; k++)
    {
        length_tree+=k*t_k_star[k-2];
        tmrca+=t_k_star[k-2];
    }
    
    printf("Length tree : %f\n",length_tree);
    printf("tmrca : %f\n",tmrca);
    */
    double** P_i_k=NULL;
    P_i_k=(double**)malloc((size_t)(taillePop-1)*sizeof(double*));
    
    for (k=0 ; k<(taillePop-1) ; k++)
    {
        P_i_k[k]=(double*)malloc((size_t)(taillePop-k-1)*sizeof(double));
    }
    
    fill_P_i_k(P_i_k,taillePop);
    
    sim_SFS(t_k_star,P_i_k,taillePop);
    
    free(t_k_star);
    for (k=0 ; k<(taillePop-1) ; k++)
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
 *      taillePop: sample size
 *
 *      returns: void
 */
void fill_P_i_k(double **P_i_k, int taillePop)
{
    int k=0;
    int i=0;
    
    for (k=1 ; k<taillePop ; k++)
    {
        P_i_k[k-1][0]=1/(double)(taillePop-1);
        for (i=2 ; i<(taillePop-k+1) ; i++)
        {
            P_i_k[k-1][i-1]=P_i_k[k-1][i-2]*((double)i*(double)(taillePop-k-i+1))/((double)(i-1)*(double)(taillePop-i));
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
 *      taillePop: sample size
 *
 *      returns: void
 */
void sim_SFS(double *t_k, double **P_i_k, int taillePop)
{
    int i=0;
    int k=0;
    
    double* SFS=NULL;
    SFS=(double*)malloc((size_t)(taillePop-1)*sizeof(double));
    for (i=0 ; i<(taillePop-1) ; i++)
        SFS[i]=0;
    
    //double sumSFS=0;
    
    for (k=1 ; k<taillePop ; k++)
    {
        for (i=2 ; i<(taillePop-k+2) ; i++)
        {
            SFS[k-1]+=(double)i*t_k[i-2]*P_i_k[k-1][i-2];
        }
        //sumSFS+=SFS[k-1];
    }
    
    for (k=1 ; k<taillePop ; k++)
    {
        //SFS[k-1]/=sumSFS;
        printf("%f\n",SFS[k-1]);
    }
    free(SFS);
}

