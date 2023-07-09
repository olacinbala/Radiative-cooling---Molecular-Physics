#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

int N_etats_1(int xN, int xN1);
void mode(double *M,double dE, int m);
void calcul_fichier_dos(double **T,double *M, int N,int m );
void calcul_dos(double **T,double *M, int N,int m );



int N_etats_1(int xN, int xN1)
{
    if(xN%xN1==0 && xN>0) // xN is a strictly positive integer.
    {
        return 1;
    }
    else
    {
        return 0;
    }
}



void mode(double *M,double dE, int m) // On pose le nombre et les valeurs des modes, on construit ainsi un tableau.
{
    int i=0;
    for(i=0; i<m; i++)
    {
        printf("Quelle est la valeur du mode %d ? ",i+1);
        scanf("%lf", &M[i]);
        M[i]= ceil(M[i]/dE);
        printf("%f\n",M[i]);
    }
}


void calcul_dos(double **T,double *M, int N,int m ) //To calculate the number of states for each internal energy (energy bin= 1 cm-1)
{
    int i,j,Z;

   for(i=0; i<(N+1); i++) // On initialise la matrice T
    {
        for(j=0; j<m; j++)
        {
            T[j][i]=-1;
        }
    }

    for(i=0; i<N+1; i++)
    {
        T[0][i]= N_etats_1(i,M[0]);
        for(j=1; j<m; j++)
        {

            Z=i-M[j];

            if (Z<0)
            {
                T[j][i]= T[j-1][i];
            }
            else
            {
                if(Z==0)
                {
                    T[j][i]= 1+T[j-1][i];
                }
                else
                {
                    T[j][i]=T[j][Z]+T[j-1][i];
                }
            }
        }

    }

}

double Zpe(double *M,int m)
{
    int i=0;
    double a=0;
    for(i=0;i<m;i++)
    {
        a=a+M[i];
    }
    a=a/2;

    return a;
}


int valeur(double *E, int k)
{
    if((E[k+1]-floor(E[k+1]))<=0.5)
    {
        return floor(E[k+1]);
    }
    else
    {
        return ceil(E[k+1]);
    }
}


int main()
{

    double Ener,dE;
    Ener=50000; // Internal energy in cm-1.
    dE=1; // energy bin

    int N=ceil(Ener/dE);

    printf("N=%d\n",N);


    N=N+300;

    int m=0,i=0,j=0,k=0,q=0;

    int n_etats_excite=30; // The maximum quantum number n to take into account. This number gives also the number of excited electronic states taken into account in this computation.
    int ntraj=1000; // Number of trajectory for the Monte Carlo simulation

    m=174; // number of vibrational degrees of freedom for the C60.


    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions related to vibrational radiative emission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    double *M=NULL; // To contain data from the input file of normal vibrational modes.
    M=malloc(m*sizeof(double));
    FILE *fichier_modes=NULL;
    fichier_modes= fopen("c60_modes_scale.txt","r");
    if(fichier_modes!=NULL)
    {
        for(i=0;i<m;i++)
        {
            fscanf(fichier_modes,"%lf ",&M[i]);
        }
    }
    fclose(fichier_modes);
    
    double *Ein=NULL; // To contain Einstein coefficients (vibrational emission)
    Ein=malloc(m*sizeof(double));

    FILE *fichier_IR_intensity=NULL;
    fichier_IR_intensity=fopen("c60_IRintensity.txt","r");

    if(fichier_IR_intensity!=NULL)
    {
        for(i=0;i<m;i++)
        {
            fscanf(fichier_IR_intensity,"%lf ",&Ein[i]);
            Ein[i]=1.2512E-7*Ein[i]*M[i]*M[i];
        }
    }
    fclose(fichier_IR_intensity);
    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of the vibrational density of states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    double **T=NULL; // To contain the vibrational density of states as a function of the vibrational energy.
    T=malloc(m*sizeof(double *));
    for(i=0;i<m;i++)
    {
        T[i]=malloc((N+1)*sizeof(double));
    }
    calcul_dos(T,M,N,m);  // Computation of the vibrational density of states as a function of the vibrational energy (internal energy).

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions related to electronic radiative emission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    double *E_exc=NULL; // To contain data from the input file of electronic oscillator strength (no unit)
    E_exc=malloc((n_etats_excite+1)*sizeof(double));  // First n_etats _excite excited electronic states.
    FILE *fichier_level=NULL;
    fichier_level= fopen("c60_elec_level_scale.txt","r");
    E_exc[0]=0;
    if(fichier_level!=NULL)
    {
        for(i=1;i<=n_etats_excite;i++)
        {
            fscanf(fichier_level,"%lf ",&E_exc[i]);
            E_exc[i]= 8065.6*E_exc[i]; // Convertion of energies (eV) into cm-1.
        }
    }
    fclose(fichier_level);


    double *I_exc=NULL; // To contain data from the input file of vibrational oscillator strength (in km/mol)
    I_exc=malloc((n_etats_excite+1)*sizeof(double));
    FILE *fichier_fluo_rate=NULL;
    fichier_fluo_rate= fopen("c60_fluo_rate.txt","r");
    I_exc[0]=0;
    if(fichier_fluo_rate!=NULL)
    {
        for(i=1;i<=n_etats_excite;i++)
        {
            fscanf(fichier_fluo_rate,"%lf ",&I_exc[i]);
        }
    }
    fclose(fichier_fluo_rate);

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    
    

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions of M_bis and T_bis for computing the radiative emission rates (Beyer-Swinehart algorithm) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double *M_bis=NULL; // To contain vibrational normal modes, except one.
    M_bis=malloc((m-1)*sizeof(double));
    
    double **T_bis=NULL; // To contain vibrational density of states without one vibrational normal mode..
    T_bis= malloc((m-1)*sizeof(double *));

    for(k=0; k<(m-1); k++)
    {
        T_bis[k]= malloc((N+1)*sizeof(double));
    }

    double ***A=NULL; // To contain emission rates.
    A=malloc((N+1)*sizeof(double **));


    for(i=0;i<=N;i++)
    {
        A[i]=malloc((m+n_etats_excite)*sizeof(double *));
        for(j=0;j<m+n_etats_excite;j++)
        {
            A[i][j]=malloc((n_etats_excite+1)*sizeof(double));
        }
    }


    for(j=0;j<m+n_etats_excite;j++)
    {
        for(i=0;i<=N;i++)
        {
            for(k=0;k<=n_etats_excite;k++)
            {
                A[i][j][k]=0;
            }
        }
    }
    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    

    /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of the vibrational de-excitation rates from the electronic ground state  (Beyer-Swinehart algorithm) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    int n=0,nmax=0,Z=0, r=0, x=0;

    printf("Calcul des taux de desexcitation IR...");
    for(j=0; j<m;j++)
    {
        for(k=0;k<m-1;k++)
            {
                if(k<j)
                {
                    M_bis[k]=M[k];
                }
                else
                {
                    M_bis[k]=M[k+1];
                }
            }

        for(k=0; k<(m-1); k++)
        {
            for(r=0;r<=N;r++)
            {
               T_bis[k][r]= 0;
            }
        }
        calcul_dos(T_bis,M_bis,N,m-1);

        for(i=0;i<=N;i++)
        {
            if(T[m-1][i]==0) // In this case, the probability to occupy upper vibrational states is zero.
            {
                A[i][j][0]=0;
            }
            else
            {
                nmax= floor(i/M[j]);

                for(n=0; n<=nmax; n++)
                {
                    Z=i-n*M[j];
                    if(Z>0)
                        {
                            A[i][j][0]= A[i][j][0]+(T_bis[m-2][Z]/T[m-1][i])*n;
                        }
                    else
                        {
                            A[i][j][0]= A[i][j][0]+((T_bis[m-2][Z]+1)/T[m-1][i])*n;
                        }
                }
                A[i][j][0]=fabs(A[i][j][0]*Ein[j]);
            }
            x++;
        }
    }
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of the Poincaré fluorescence rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int l=0;

    printf("calcul des taux effectif de fluorescence de Poincare...");

    for(l=1;l<=n_etats_excite;l++)
    {
        if(N<E_exc[l])
        {
        }
        else
        {
            for(i=0;i<=N;i++)
            {
                if(i<E_exc[l])
                {

                }
                else
                {
                    if(T[m-1][i]==0)
                    {
                        // In that case, we know that A[i][m]=0.
                    }
                    else
                    {
                        A[i][m-1+l][0]= 0.665*E_exc[l]*E_exc[l]*I_exc[l]*T[m-1][i-valeur(E_exc,l-1)]/T[m-1][i]; //To use if Poincaré fluorescence is taken into account in the relaxation.
                       // A[i][m-1+l][0]= 0; //To use if Poincaré fluorescence is NOT taken into account in the relaxation.
                    }
                }

            }
        }
    }
   ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






   ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of the de-excitation rates from the excited electronic state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    for(l=1;l<=n_etats_excite;l++)
    {
        for(i=0;i<=N;i++)
            {

                for(j=0;j<m;j++)
                {
                    if(T[m-1][i]==0 || i<valeur(E_exc,l-1))
                    {
                    }
                    else
                    {
                        A[i][j][l]=A[i][j][0]*T[m-1][i-valeur(E_exc,l-1)]/T[m-1][i];
                        if(l==1)
                        {
                            
                        }
                    }
                }

            }
    }

    free(Ein);
    Ein=NULL;
    free(M_bis);
    M_bis=NULL;
    free(T);
    T=NULL;
    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N=N-300;
    int taille_max= ceil(Ener/M[0]);
    taille_max=1000;    // This integer will define the size of the following tables used for individual Monte Carlo trajectories. 
                    // Usually, one Monte Carlo trajectorycontains less than 1000 steps.
                    
                    
 //   printf("taille_max= %d\n",taille_max);



    double *E; // To contain the energies at each step of one Monte Carlo trajectory simulating the radiative relaxation.
    double *t; // To contain the time at each step of one Monte Carlo trajectory simulating the radiative relaxation.
    E=NULL;
    t=NULL;
    t=malloc(taille_max*sizeof(double));
    E=malloc(taille_max*sizeof(double));

    E[0]=N; // Initial conditions.
    t[0]=0;

    for(i=1;i<taille_max;i++)
    {
        E[i]=0;
        t[i]=0;
    }



    int ***IR_par_modes; // To contain the vibrational emission spectrum resulting from all Monte Carlo trajectories.
    IR_par_modes=NULL;
    IR_par_modes=malloc(taille_max*sizeof(double **)); // This dimension is used for following relaxation kinetics.
    for(i=0;i<taille_max;i++)
    {
        IR_par_modes[i]=malloc(m*sizeof(int *)); // This dimension is used for following the emission from each vibrational normal modes in the GROUND electronic states.
        for(j=0;j<m;j++)
        {
          IR_par_modes[i][j]=malloc((n_etats_excite+1)*sizeof(int)); // This dimension is used for following the emission from each vibrational normal modes in EXCITED electronic states.
        }

    }


    for(k=0;k<=n_etats_excite;k++) // Initialisation of the table IR_par_modes.
    {
        for(j=0;j<m;j++)
        {
            for(i=0;i<taille_max;i++)
            {
                IR_par_modes[i][j][k]=0;
            }
        }
    }


    double *Somme_E;
    double *Somme_t;
    Somme_t=NULL;
    Somme_E=NULL;
    Somme_E=malloc(taille_max*sizeof(double));
    Somme_t=malloc(taille_max*sizeof(double));
    for(i=0;i<taille_max;i++)
    {
        Somme_E[i]=0;
        Somme_t[i]=0;
    }

    int **Signal_IR_states;
    Signal_IR_states=NULL;
    Signal_IR_states=malloc(taille_max*sizeof(int *));
    for(i=0;i<taille_max;i++)
    {
        Signal_IR_states[i]=malloc((n_etats_excite+1)*sizeof(int));
    }
    for(i=0;i<taille_max;i++)
    {
        for(l=0;l<=n_etats_excite;l++)
        {
            Signal_IR_states[i][l]=0;
        }
    }

    srand(time(NULL));
    double a=0, a1=0,a2=0,p=0;
    int je; 

    int En=N;
    int c=0;
    r=0;
    double SA_kmc=0;
    long long int SA_kmc_entier=0;
    long long int nombre_hasard=0;


    int *C;
    C=malloc(taille_max*sizeof(int));

    int **C_UV;
    C_UV=malloc(taille_max*sizeof(int *));
    for(i=0;i<taille_max;i++)
    {
        C_UV[i]=malloc((n_etats_excite+1)*sizeof(int));
    }

    int **C_IR;
    C_IR=malloc(taille_max*sizeof(int *));
    for(i=0;i<taille_max;i++)
    {
        C_IR[i]=malloc((n_etats_excite+1)*sizeof(int));
    }

    int d=0;

    int *nombre_poincare; // To count the number of emitted photon from EACH electronic states over all the relaxation mechanism.
    nombre_poincare=malloc((n_etats_excite+1)*sizeof(double));
    for(i=0;i<=n_etats_excite;i++)
    {
        nombre_poincare[i]=0;
    }




    /// UV_En
    double **UV_En=NULL; //To contain the internal energy and for simulating the Poincaré fluorescence kinetics.
    UV_En=malloc(taille_max*sizeof(double *));
    for(i=0;i<taille_max;i++)
    {
        UV_En[i]=malloc((n_etats_excite+1)*sizeof(double));
    }


    double **Somme_UV_En=NULL; // To sum the internal energy at each step for all kinetic Monte Carlo (kMC) trajectories. This will allow to simulate the mean kMC trajectory.
    Somme_UV_En=malloc(taille_max*sizeof(double *));

    for(i=0;i<taille_max;i++)
    {
        Somme_UV_En[i]=malloc((n_etats_excite+1)*sizeof(double));
    }


    printf("UV_EN OK\n");

    /// UV_t 
    double **UV_t=NULL; // To contain the time and for simulating the Poincaré fluorescence kinetics.
    UV_t=malloc(taille_max*sizeof(double *));
    for(i=0;i<taille_max;i++)
    {
        UV_t[i]=malloc((n_etats_excite+1)*sizeof(double));
    }



    double **Somme_UV_t=NULL; // To sum the time at each step for all kinetic Monte Carlo (kMC) trajectories. This will allow to simulate the mean kMC trajectory.
    Somme_UV_t=malloc(taille_max*sizeof(double *));

    for(i=0;i<taille_max;i++)
    {
        Somme_UV_t[i]=malloc((n_etats_excite+1)*sizeof(double));
    }


    printf("UV_t OK\n");

    /// UV_compteur 
    int **UV_compteur=NULL; // To count the number of emitted photons (via Poincaré fluorescence) at each Monte Carlo step from EACH excited electronic states.
    UV_compteur=malloc(taille_max*sizeof(int *)); 
    for(i=0;i<taille_max;i++)
    {
        UV_compteur[i]=malloc((n_etats_excite+1)*sizeof(int));
    }



    int **Somme_UV_compteur=NULL; // To sum the number of emitted photons (via Poincaré fluorescence) at each Monte Carlo step from ALL excited electronic states.
    Somme_UV_compteur=malloc(taille_max*sizeof(int *));

    for(i=0;i<taille_max;i++)
    {
        Somme_UV_compteur[i]=malloc((n_etats_excite+1)*sizeof(int));
    }


    printf("UV_comtpeur OK\n");

    /// IR_signal_niveaux_En

    double **IR_signal_niveaux_En=NULL; // To contain the internal energy at each step and for simulating the vibrational emission kinetics.
    IR_signal_niveaux_En=malloc(taille_max*sizeof(double *));
    for(i=0;i<taille_max;i++)
    {
        IR_signal_niveaux_En[i]=malloc((n_etats_excite+1)*sizeof(double));
    }


    double **Somme_IR_signal_niveaux_En=NULL; // To sum the internal energy at each step for all kinetic Monte Carlo (kMC) trajectories. This will allow to simulate the mean kMC trajectory.
    Somme_IR_signal_niveaux_En =malloc(taille_max*sizeof(double *));

    for(i=0;i<taille_max;i++)
    {
        Somme_IR_signal_niveaux_En[i]=malloc((n_etats_excite+1)*sizeof(double));
    }

    printf("IR_signalniveaux_En OK\n");
    /// IR_signal_niveaux_t

    double **IR_signal_niveaux_t=NULL; // To contain the time at each step and for simulating the vibrational emission kinetics.
    IR_signal_niveaux_t=malloc(taille_max*sizeof(double *));
    for(i=0;i<taille_max;i++)
    {
        IR_signal_niveaux_t[i]=malloc((n_etats_excite+1)*sizeof(double));
    }


    double **Somme_IR_signal_niveaux_t=NULL; // To sum the time at each step for all kinetic Monte Carlo (kMC) trajectories. This will allow to simulate the mean kMC trajectory.
    Somme_IR_signal_niveaux_t =malloc(taille_max*sizeof(double *));

    for(i=0;i<taille_max;i++)
    {
        Somme_IR_signal_niveaux_t[i]=malloc((n_etats_excite+1)*sizeof(double));
    }
    printf("IR_signalniveaux_t OK\n");

    /// IR_signal_niveaux_compteur

    int **IR_signal_niveaux_compteur=NULL; // To count the number of emitted photons (via vibrational emission) at each Monte Carlo step from EACH excited electronic states.
    IR_signal_niveaux_compteur=malloc(taille_max*sizeof(int *)); 
    for(i=0;i<taille_max;i++)
    {
        IR_signal_niveaux_compteur[i]=malloc((n_etats_excite+1)*sizeof(int));
    }


    int **Somme_IR_signal_niveaux_compteur=NULL; // To sum the number of emitted photons (via vibrational emission) at each Monte Carlo step from ALL excited electronic states.
    Somme_IR_signal_niveaux_compteur =malloc(taille_max*sizeof(int *));

    for(i=0;i<taille_max;i++)
    {
        Somme_IR_signal_niveaux_compteur[i]=malloc((n_etats_excite+1)*sizeof(int));
    }

    printf("IR_signalniveaux_compteur OK\n");


    for(i=0;i<taille_max;i++) // Initialisation of all above-defined tables. 
    {
        for(j=0;j<=n_etats_excite;j++)
        {
            Somme_UV_compteur[i][j]=0;
            UV_compteur[i][j]=0;
            Somme_UV_t[i][j]=0;
            UV_t[i][j]=0;
            Somme_UV_En[i][j]=0;
            UV_En[i][j]=0;
            C[i]=0;
            C_UV[i][j]=0;
            IR_signal_niveaux_compteur[i][j]=0;
            IR_signal_niveaux_En[i][j]=0;
            IR_signal_niveaux_t[i][j]=0;
            Somme_IR_signal_niveaux_compteur[i][j]=0;
            Somme_IR_signal_niveaux_En[i][j]=0;
            Somme_IR_signal_niveaux_t[i][j]=0;
        }
    }



    i=0;
    int MC=0;
    double MC_reel=0;
    FILE *fichier3=NULL;
    fichier3= fopen("probleme.txt","w");

    printf("calcul des emissions infrarouge...\n");
    
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kinetic Monte Carlo simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for(k=1;k<=ntraj;k++)
    {
        En=N;
        i=0;


        for(j=0;j<taille_max;j++) // Set to zero for each individual Monte Carlo trajectory.
        {
            for(q=0;q<=n_etats_excite;q++)
            {
                UV_En[j][q]=0;
                UV_t[j][q]=0;
                UV_compteur[j][q]=0;

                IR_signal_niveaux_En[j][q]=0;
                IR_signal_niveaux_t[j][q]=0;
                IR_signal_niveaux_compteur[j][q]=0;
            }
        }


        do
        {

            SA_kmc=0;
            for(l=0;l<=n_etats_excite;l++)
            {
                for(r=0;r<m+n_etats_excite;r++)
                {
                    SA_kmc= SA_kmc+A[En][r][l]; // The random number for the Monte Carlo simulation is chosen between 0 and SA_kmc.
                }
            }

            
            ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% To determine the random number a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(SA_kmc==0)
            {
                break;
            }



            MC= (int) (1000*SA_kmc);
            if(MC>1)
            {
                do
                {
                    a1=rand()%MC;
                }
                while(a1==0);
            }
            else
            {
                a1=0;
            }
            a1=a1/1000;





            MC_reel=MC;
            MC= (int) ((SA_kmc-MC_reel/1000)*1000000);

            if(MC==0)
            {
                a2=0;
            }
            else
            {
               a2=rand()%MC;
               a2=a2/1000000;
            }


            /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            a=a1+a2;

            p=0;
            j=-1;
            c=-1;
            l=0;
            do // This loop allows to randomly select the relaxation pathways (vibrational normal mode or Poincaré fluorescence from an excited electronic states).
            {
                if(j==m-1+n_etats_excite)
                    {
                        l++;
                        j=-1;
                    }
                j++;
                p=p+A[En][j][l];

                c++;
            }
            while(a>p && c<((m-1+n_etats_excite)*(n_etats_excite+1)));



                if(j>=m && l==0) // Poincaré fluorescence
                {
                    E[i+1]=E[i]-E_exc[j-m+1]; // E_exc[j-m+1] is the energy of the emitted photon.
                    t[i+1]=1/SA_kmc; // Mean temporal in kinetic Monte Carlo (see the theory for more details)

                    UV_En[i][j-m+1]=E[i];
                    UV_t[i][j-m+1]=1/SA_kmc;
                    UV_compteur[i][j-m+1]++;
                    C_UV[i][j-m+1]++;

                    nombre_poincare[j-m+1]++; //
                    En=E[i+1];

                }
                else // Vibrational emission
                {

                    je=j;

                    IR_par_modes[i][je][l]++;

                    E[i+1]=E[i]-M[je]; // M[je] is the energy of the emitted photon.
                    t[i+1]=1/SA_kmc;

                    IR_signal_niveaux_En[i][l]=E[i];

                    IR_signal_niveaux_t[i][l]=1/SA_kmc;
                    IR_signal_niveaux_compteur[i][l]++;
                    C_IR[i][l]++;

                }
                En=valeur(E,i);
                i++;
        }
        while(En>=1000);


        for(d=1;d<=i;d++)
        {
            C[d]=C[d]+1; // To count the number of trajectory for each step "d". Useful for determining the mean kinetic Monte Carlo trajectory.
        }

        for(j=0;j<=i;j++)
        {
            Somme_E[j]=E[j]+Somme_E[j];
            Somme_t[j]=t[j]+Somme_t[j];

            for(q=0;q<=n_etats_excite;q++)
            {
               Somme_UV_En[j][q]=UV_En[j][q]+Somme_UV_En[j][q];
               Somme_UV_t[j][q]=UV_t[j][q]+Somme_UV_t[j][q];
               Somme_UV_compteur[j][q]=UV_compteur[j][q]+Somme_UV_compteur[j][q];

               Somme_IR_signal_niveaux_En[j][q]=IR_signal_niveaux_En[j][q] + Somme_IR_signal_niveaux_En[j][q];
               Somme_IR_signal_niveaux_t[j][q]=IR_signal_niveaux_t[j][q] + Somme_IR_signal_niveaux_t[j][q];
               Somme_IR_signal_niveaux_compteur[j][q]=IR_signal_niveaux_compteur[j][q] + Somme_IR_signal_niveaux_compteur[j][q];
            }
        }
    }

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E[0]=N;

    t[0]=0;
    for(q=0;q<=n_etats_excite;q++)
    {
        UV_En[0][q]=N;
        Somme_IR_signal_niveaux_En[0][q]=N;
        UV_t[0][q]=0;
        Somme_IR_signal_niveaux_t[0][q]=0;

    }

    ///%%%%%%%%%%%%%%%%%%%%% Mean kinetic Monte Carlo trajectory %%%%%%%%%%%%%%%%%%%%%%
    
    for(i=1;i<taille_max;i++)
    {
        E[i]=Somme_E[i]/C[i];
        t[i]=Somme_t[i]/C[i];

        for(q=0;q<=n_etats_excite;q++)
        {
            UV_En[i][q]=Somme_UV_En[i][q]/C_UV[i][q];
            UV_t[i][q]=Somme_UV_t[i][q]/C_UV[i][q];

            IR_signal_niveaux_En[i][q]=Somme_IR_signal_niveaux_En[i][q]/C_IR[i][q];
            IR_signal_niveaux_t[i][q]= Somme_IR_signal_niveaux_t[i][q]/C_IR[i][q];
        }
    }


    for(i=1;i<taille_max;i++)
    {
        t[i]=t[i-1]+t[i];
        for(q=0;q<=n_etats_excite;q++)
        {
            UV_t[i][q]=UV_t[i-1][q]+UV_t[i][q];
            IR_signal_niveaux_t[i][q]=IR_signal_niveaux_t[i-1][q]+IR_signal_niveaux_t[i][q];
        }
    }
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    FILE *fichier_decroissance= NULL;
    fichier_decroissance= fopen("decroissance_E.txt", "w");

    i=0;
    do
    {
        fprintf(fichier_decroissance,"%f    %f\n",t[i],E[i]);
        i++;
    }
    while(E[i]>=1000);
    fclose(fichier_decroissance);

    
    
    /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Emission from each vibrational normal modes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FILE *fichier_emission=NULL;
    fichier_emission=fopen("emissionIR_par_modes.txt","w"); // Number of emitted photon for each vibrational normal modes (at the ground electronic states) at each Monte Carlo step.


    fprintf(fichier_emission,"time     ");
    for(j=0;j<m;j++)
    {
        fprintf(fichier_emission,"%f        ",M[j]);
    }
    fprintf(fichier_emission,"\n");

    for(i=0;i<taille_max;i++)
    {
       fprintf(fichier_emission,"%f            ",t[i]);
       for(j=0;j<m;j++)
        {
                fprintf(fichier_emission,"%d                ",IR_par_modes[i][j][0]);
        }
        fprintf(fichier_emission,"\n");
    }
    fprintf(fichier_emission,"\n\n\n");

    fprintf(fichier_emission,"nombre de trajectoires= %d\n",k-1);
    int total=0;
    k=0;
    for(l=0;l<=n_etats_excite;l++)
    {
        total=0;
        for(j=0;j<m;j++)
        {
            for(i=0;i<taille_max;i++)
            {
                total=IR_par_modes[i][j][l]+total;
            }
        }
        fprintf(fichier_emission,"nombres de photons IR émis depuis l'état %d= %d\n",l,total);
        k=total+k;
    }

    fprintf(fichier_emission,"nombres total de photons IR émis= %d\n",k);


    FILE *fichier_emissionIRglobal=NULL;
    fichier_emissionIRglobal=fopen("emissionIR_par_modes_global.txt","w");

    for(j=0;j<m;j++)
    {
        k=0;
        fprintf(fichier_emissionIRglobal,"%f        ",10000/M[j]);
        for(l=0;l<=n_etats_excite;l++)
        {
           for(i=0;i<100;i++)
           {
                k=k+IR_par_modes[i][j][l];
           }
        }
        fprintf(fichier_emissionIRglobal,"%d        ",k);
        fprintf(fichier_emissionIRglobal,"\n");
    }

    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Emission via Poincaré fluorescence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FILE *fichier_spectreUV=NULL;
    fichier_spectreUV= fopen("spectreUV.txt","w");
    for(i=1;i<=n_etats_excite;i++)
    {
        fprintf(fichier_spectreUV,"%f        %d\n",10000/E_exc[i],nombre_poincare[i]);
    }

    total=0;
    for(i=0;i<=n_etats_excite;i++)
    {
        total=nombre_poincare[i]+total;
    }
    fprintf(fichier_emission,"nombre de photons UV poincaré émis:%d\n",total);


    FILE *fichier_emissionUV= NULL;
    fichier_emissionUV= fopen("emissionUV.txt", "w");
    i=0;
    do
    {
        for(j=1;j<=n_etats_excite;j++)
        {
            fprintf(fichier_emissionUV,"%f     %f      %d               ",UV_t[i][j],UV_En[i][j],C_UV[i][j]);
        }
        fprintf(fichier_emissionUV,"\n");
        i++;
    }
    while(UV_En[i][1]>=1000);








    /// emissionIR_niveaux_elec.txt
    i=0;
    FILE *fichier_emissionIR=NULL;
    fichier_emissionIR= fopen("emissionIR_niveaux_elec.txt","w");
    do
    {
        for(j=0;j<=n_etats_excite;j++)
        {
            fprintf(fichier_emissionIR,"%f     %f      %d               ",IR_signal_niveaux_t[i][j],IR_signal_niveaux_En[i][j],C_IR[i][j]);
        }
        fprintf(fichier_emissionIR,"\n");
        i++;
    }
    while(IR_signal_niveaux_En[i][0]>=1000);

    fclose(fichier_emission);
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return 0;
}


