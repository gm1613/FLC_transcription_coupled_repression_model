//
//  20230117_autopathway_flc.cpp
//  
//
//  Created by Govind Menon (JIC) on 17/01/2023.
//
//
// clang++ -Wall -pedantic 20230117_autopathway_flc.cpp -o 20230117_autopathway_flc  -lgsl -lgslcblas
//  ./20230117_autopathway_flc 5 5 1 0.5 (3/0/-1) typeofsimulation

//  Simulations for combined model of Transcription-coupled repression through proximal termination and Polycomb silencing at FLC, Menon et al., 2023

//#include <omp.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>



// SETUP FUNCTION THAT GETS A SEED FOR THE RNG FROM THE SYSTEM CLOCK
unsigned long int random_seed()
{
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
}



// MAIN WITH OUTPUTS WRITTEN DIRECTLY TO TEXT FILES, TAKES SIX USER INPUTS: (1) NUMBER OF SIMULATIONS, (2) DURATION OF INDIVIDUAL SIMULATION (NUMBER OF CELL CYCLES), (3) FCA EFFECTIVENESS PARAMETER, (4) FLD MEDIATED DEMETHYLATION PROBABILITY (PER H3 HISTONE PER PROXIMAL TERMINATION EVENT), (5) INITIAL STATE OF H3 HISTONES (UNIFORM ACROSS LOCUS), (6) NAME TAG FOR OUTPUT FILES, INDICATING TYPE OF SIMULATION (Ler, fca-3, fca-1 etc.)
int main(int argc, char *argv[])
{
    
    // Name tag for output files
    std::string out_file_tag = "20230117_autopathway_flc_";
    
    // SET UP TIMER
    clock_t start,end;
    double cpu_time_used;
    int nthreads;
    int tid;
    
    int flag=1;
    start=clock();
    
    //SETTING USER INPUT PARAMETERS
    int nsims = atoi(argv[1]);
    int ncycles = atoi(argv[2]);
    float fca_parameter = atof(argv[3]); //
    float pdem_fld = atof(argv[4]);
    //  INITIAL CONDITIONS
    int init_H3state = atoi(argv[5]);
    // OUTPUT FILENAME TAG INDICATING TYPE OF SIMULATION
    std::string out_file_tag_plus = argv[6];

    
    
    // SYSTEM PARAMETERS
    // SYSTEM SIZE
    const int N = 2*30; //Number of H3 histones at FLC
    const int nrxn = 4*N+2; // total number of reactions: K4me +,- K27me +,-, at each H3 histone, Sense initiation, AS initiation
    
    // REACTION PARAMETERS
    // FREE PARAMETERS
    const float kme = 8e-6; //PRC2 mediated me2 to me3 methylation rate (hist[-1]second[-1])

    const float fmax = 5.3e-4; //maximum transcription frequency (second[-1])

    const float pdem = 25*0.004; //0.004 //demethylation probability per histone during transcription (hist[-1]transcription[-1])
    const float pex = 4e-3; // histone exchange proability per histone during transcription (hist[-1]transcription[-1])
    const float pinherit = 0.5; // H3 pair inheritance probability at replication
    const float PT27 = 1.f/3; // K27 me2/me3 repression threshold: Set PT27<=1 to incorporate threshold behaviour
    
    // DEPENDENT PARAMETERS
    
    const float kme01 = 9*kme;
    const float kme12 = 6*kme;
    const float gme01 = kme01/20; // Background (noisy) methylation of H3K27
    const float gme12 = kme12/20; // Background (noisy) methylation of H3K27
    const float gme23 = kme/20;   // Background (noisy) methylation of H3K27
    const float rho_me2 = 0.1;   // Parameter capturing H3K27me2 contribution to PRC2 read-write feedback
    const float fmin = fmax/60; // minimum transcription frequency (second[-1])
    const float gdem = fmin*pdem; //Background (noisy) demethylation (hist[-1]second[-1])
    
    const float gme01_nr = 7*kme01/20; // De-novo methylation of H3K27 in nucleation region
    const float gme12_nr = 7*kme12/20; // De-novo methylation of H3K27 in nucleation region
    const float gme23_nr = 7*kme/20; // De-novo methylation of H3K27 in nucleation region
    
    // PROXIMAL VS DISTAL DECISION PARAMETERS
    const float delta_proximal = 0.05; // basal probability of proximal termination of AS transcript
    const float mu_distal = 0.05; // basal probability of distal termination of AS transcript
    
    // Trans activation PARAMETERS
    
    const float alpha_sense = 1; // Trans-factor activation parameter for sense transcription
    const float alpha_antisense = 0.2; // Trans-factor activation parameter for antisense transcription
    
    // Polycomb activity parameter
    const float beta = 0; // Set to 0.7 (previously 0.55) as default; set to zero for simulations of analog module in isolation
    
    
    // Cell Cycle duration PARAMETERS
    const float cyclehr = 22; //cycle length in HOURS (based on Rahni,2019)
    const float cycletime = cyclehr*3600; //cell cycle length in seconds
    
    // Refractory periods for mutual exclusivity
    const float tau_S = 6*60; // refractory period after sense initiation event (6mins)
    const float tau_AS = 6*60; // refractory period after anti-sense initiation event (6mins)
    
    //Defining the nucleation region and looping region (nucleosomes 1-3 and 19-21)
    const int nr_start = 0;
    const int nr_end = 5;
    const int lr_start = 39;
    const int lr_end = 44;
    const int N_NR = 6+6;
 
    
    //Defining OFF criterion parameter
    const float chi = 0.25;
    
    //Contact parameters for clf

        const float xi_nr_nr = 0.8; //Unchanged in clf
        const float xi_nr_lr = 0.8; //Unchanged in clf
        const float xi_nr_body = 0.12; //This is set to a lower value (0.12) for clf, 0.2 otherwise
        const float xi_lr_nr = 0.8; //Unchanged in clf
        const float xi_lr_lr = 0.8;//Unchanged in clf
        const float xi_lr_body = 0.12; //This is set to a lower value (0.12) for clf, 0.2 otherwise
        const float xi_body_nr = 0.05; //This is set to a lower value (0.01) for clf, 0.5 otherwise (NR contribution TO body region)
        const float xi_body_lr = 0.05; //This is set to a lower value (0.01) for clf, 0.5 otherwise (LR contribution TO body region)
        const float xi_body_body = 0.01; //This is set always to 0.01

    //Contact parameters for WT
//        const float xi_nr_nr = 1; //Unchanged in clf
//        const float xi_nr_lr = 1; //Unchanged in clf
//        const float xi_nr_body = 0.2; //This is set to a lower value (0.12) for clf, 0.2 otherwise
//        const float xi_lr_nr = 1; //Unchanged in clf
//        const float xi_lr_lr = 1;//Unchanged in clf
//        const float xi_lr_body = 0.2; //This is set to a lower value (0.12) for clf, 0.2 otherwise
//        const float xi_body_nr = 0.5; //This is set to a lower value (0.01) for clf, 0.5 otherwise
//        const float xi_body_lr = 0.5; //This is set to a lower value (0.01) for clf, 0.5 otherwise
//        const float xi_body_body = 0.01; //This is set always to 0.01
    

    
    // WEIGHTING PARAMETER FOR COMPUTING EFFECTIVE H3K27ME3 COVERAGE THAT DETERMINES REPRESSION OF TRANSCRIPTION
    const float eta = 0.25; // Full coverage of the NR by itself contributes 0.25 (if all nucleosomes were equal, this would be 0.2)
    
    // FRACTION OF TOTAL CELL CYCLE WHERE SPREADING IS POSSIBLE
    const float spreading_phase = 0.5;
    
    
    // K4 PARAMETERS
    const int K4_switch = 1; // Used to switch on(1)/off(0) the K4me reactions
    const float lambda = 3.128/0.959; //Ratio of K27me3 half-life to K4me1 half-life (Zee et al., JBC, 2010)
    const float pk4me = K4_switch*100*kme/fmax;//Transcription mediated K4 methylation probability per histone during transcription (hist[-1]transcription[-1])
    const float gk4me = K4_switch*kme; //Background (noisy) methylation rate
    const float gk4dem = 0.4*K4_switch*((lambda-1)*( fmax*pex + (1-pinherit)/cycletime) + lambda*fmax*pdem + lambda*gdem)*1; //Background (noisy) demethylation rate
    
    

    
    

    
    //##########################################################################################################################################################################################//
    //WRITING PARAMETER VALUES TO OUTPUT FILE
    std::string parameter_values_file_name = out_file_tag + "parameters_" + out_file_tag_plus + ".txt";
    std::ofstream file0_;
    file0_.open(parameter_values_file_name);
    if(file0_.is_open())
    {
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "USER INPUT PARAMETERS" << "\n";
        file0_ << "nsims = "<< nsims << "\n";
        file0_ << "ncycles = "<< ncycles << "\n";
        file0_ << "fca_parameter = "<< fca_parameter << "\n";
        file0_ << "pdem_fld = "<< pdem_fld << "\n";
        file0_ << "init_H3state = "<< init_H3state << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "SYSTEM SIZE" << "\n";
        file0_ << "N_histones = "<< N << "\n";
        file0_ << "nrxn = "<< nrxn << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "REACTION PARAMETERS FOR H3K27ME3" << "\n";
        file0_ << "kme = "<< kme << "\n";

        file0_ << " pdem = "<< pdem << "\n";
        file0_ << " pex = "<< pex << "\n";
        file0_ << " pinherit = "<< pinherit << "\n";
        
        file0_ << " kme01 = "<< kme01 << "\n";
        file0_ << " kme12 = "<< kme12 << "\n";
        file0_ << " gme01 = "<< gme01 << "\n";
        file0_ << " gme12 = "<< gme12 << "\n";
        file0_ << " gme23 = "<< gme23 << "\n";
        file0_ << " rho_me2 = "<< rho_me2 << "\n";
        
        file0_ << " gdem = "<< gdem << "\n";
        file0_ << " gme01_nr = "<< gme01_nr << "\n";
        file0_ << " gme12_nr = "<< gme12_nr << "\n";
        file0_ << " gme23_nr = "<< gme23_nr << "\n";
        file0_ << " beta = "<< beta << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "DEFINITION OF H3K27ME3 NUCLEATION AND LOOPING REGIONS AND ASSOCIATED INTERACTION PARAMETERS" << "\n";
        file0_ << " nr_start = "<< nr_start << "\n";
        file0_ << " nr_end = "<< nr_end << "\n";
        file0_ << " lr_start = "<< lr_start << "\n";
        file0_ << " lr_end = "<< lr_end << "\n";
        file0_ << " xi_nr_nr = "<< xi_nr_nr << "\n";
        file0_ << " xi_nr_lr = "<< xi_nr_lr << "\n";
        file0_ << " xi_nr_lr = "<< xi_nr_body << "\n";
        file0_ << " xi_lr_nr = "<< xi_lr_nr << "\n";
        file0_ << " xi_lr_lr = "<< xi_lr_lr << "\n";
        file0_ << " xi_lr_body = "<< xi_lr_body << "\n";
        file0_ << " xi_body_nr = "<< xi_body_nr << "\n";
        file0_ << " xi_body_lr = "<< xi_body_lr << "\n";
        file0_ << " xi_body_body = "<< xi_body_body << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "WEIGHTING PARAMETER FOR COMPUTING EFFECTIVE H3K27ME3 COVERAGE THAT DETERMINES REPRESSION OF TRANSCRIPTION" << "\n";
        file0_ << " eta = "<< eta << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "FRACTION OF TOTAL CELL CYCLE WHERE SPREADING IS POSSIBLE" << "\n";
        file0_ << " spreading_phase = "<< spreading_phase << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "OFF CRITERION PARAMETER" << "\n";
        file0_ << " chi = "<< chi << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "REACTION PARAMETERS FOR TRANSCRIPTION" << "\n";
        file0_ << "fmax = "<< fmax << "\n";
        file0_ << " fmin = "<< fmin << "\n";
        file0_ << " alpha_sense = "<< alpha_sense << "\n";
        file0_ << " alpha_antisense = "<< alpha_antisense << "\n";
        file0_ << " PT27 = "<< PT27 << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "PROXIMAL AND DISTAL TERMINATION PARAMETERS" << "\n";
        file0_ << " delta_proximal = "<< delta_proximal << "\n";
        file0_ << " mu_distal = "<< mu_distal << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "REFRACTORY PERIODS FOR MUTUAL EXCLUSIVITY OF SENSE AND ANTISENSE TRANSCRIPTION" << "\n";
        file0_ << " tau_S = "<< tau_S << "\n";
        file0_ << " tau_AS = "<< tau_AS << "\n";
        file0_ << "\n";
        file0_ << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        file0_ << "REACTION PARAMETERS FOR H3K4 METHYLATION" << "\n";
        file0_ << " K4_switch = "<< K4_switch << "\n";
        file0_ << " lambda = "<< lambda << "\n";
        file0_ << " pk4me = "<< pk4me << "\n";
        file0_ << " gk4me = "<< gk4me << "\n";
        file0_ << " gk4dem = "<< gk4dem << "\n";

        file0_.close();
    }
    
    
    
    // SIMULATION and ANALYSIS PARAMETERS
    
    double tmax = (ncycles*cycletime); //TIME LIMIT FOR SIMULATION
    float avg_start_cycle = 0; // Averaging starts at end of this cycle
    float avg_start_time = avg_start_cycle*cycletime; // Start time for averaging
    
    
    //#######################################################################################################################################################################################//
    // DEFINE ARRAYS CONTAINING TIMEPOINTS, AND CORRESPONDING SUMS, FOR AVERAGING
    const int ntimes = 50; //50 timepoints per cell cycle
    const int t_traj_size = ntimes*(ncycles-avg_start_cycle) + 1;
    float *t_traj = new float[t_traj_size];
    for(int i=0;i<t_traj_size;i++)
    {
        t_traj[i] = avg_start_time + i*(cycletime/ntimes);
    }
    
    
    int t_term_freq_size = (((ncycles-avg_start_cycle)*cyclehr)/1) + 1;
    float *t_term_freq = new float[t_term_freq_size];
    for(int i=0;i<t_term_freq_size;i++)
    {
        t_term_freq[i] = avg_start_time + i*3600;
    }
    
    
  
    

    tid = 1; // Set thread number (only important for parallel implementation using OpenMPI
                               
        
        // SET UP pseudo-random number generator
        
        gsl_rng *r=gsl_rng_alloc(gsl_rng_mt19937);
        unsigned long int seed_value = random_seed();
        gsl_rng_set(r, seed_value);
        
        
        
        // SET OUTPUT FILENAMES
        std::string type1a_file_name= out_file_tag + "type1a_tid" + std::to_string(tid) + "_" + out_file_tag_plus + ".txt";
        std::string type1b_file_name= out_file_tag + "type1b_tid" + std::to_string(tid) + "_" + out_file_tag_plus + ".txt";
        std::string type1c_file_name= out_file_tag + "type1c_tid" + std::to_string(tid) + "_" + out_file_tag_plus + ".txt";
        std::string type1d_file_name= out_file_tag + "type1d_tid" + std::to_string(tid) + "_" + out_file_tag_plus + ".txt";
        
        //  SYSTEM VARIABLES
        // TIME
        double t;
        double t_prev,reptime;
        
        
        
        
        // STATE VARIABLES
        int *cdomK27K4 = new int[N]; //ARRAY OF HISTONE H3 K27/K4 METHYLATION STATES: {3} for K27me1/2/3 and {-1} for K4me1 and {0} for no modification
        float PK27me3,PK27me3_nr,PK27me3_whole,PK27me3_rem, PK4me1; //FRACTIONAL COVERAGE OF MODIFICATIONS
        float P_proximal;
        float omega_t;
        
        // Analysis variables
        float *sumovertraj_cdomK27 = new float[N]; //ARRAY OF FINAL H3K27ME3 LEVELS AT EACH HISTONE, AVERAGED OVER TRAJECTORIES
        float *sumovertraj_cdomK4 = new float[N]; //ARRAY OF FINAL H3K4ME1 LEVELS AT EACH HISTONE, AVERAGED OVER TRAJECTORIES
//        float *sumovertraj_K27me_prpnsty = new float[N];
        
        for(int i=0;i<N;i++)
        {
            sumovertraj_cdomK27[i] = 0;
            sumovertraj_cdomK4[i] = 0;
//            sumovertraj_K27me_prpnsty[i] = 0;
        }
        
        
        float *sumovertraj_K27me3_nr = new float[t_traj_size];
        float *sumovertraj_K27me3_whole = new float[t_traj_size];
        float *sumovertraj_K4me1 = new float[t_traj_size];
        float *sumovertraj_OFF_state = new float[t_traj_size];
        
        for(int i=0;i<t_traj_size;i++)
        {
            sumovertraj_K27me3_nr[i] = 0;
            sumovertraj_K27me3_whole[i] = 0;
            sumovertraj_K4me1[i] = 0;
            sumovertraj_OFF_state[i] = 0;
            
        }

        
        float *sumovertraj_sense_prox_events = new float[t_term_freq_size];
        float *sumovertraj_sense_dist_events = new float[t_term_freq_size];
        float *sumovertraj_antisense_prox_events = new float[t_term_freq_size];
        float *sumovertraj_antisense_dist_events = new float[t_term_freq_size];
        
        for(int i=0;i<t_term_freq_size;i++)
        {
            sumovertraj_sense_prox_events[i] = 0;
            sumovertraj_sense_dist_events[i] = 0;
            sumovertraj_antisense_prox_events[i] = 0;
            sumovertraj_antisense_dist_events[i] = 0;
        }
        
        float *switchoff_time = new float[nsims];
        for(int i=0;i<nsims;i++)
        {
            switchoff_time[i] = 0;
        }

        
        
        // Other variables
        float Ei,nr_contribution,lr_contribution,body_contribution;
        float nr_contribution_corr,lr_contribution_corr,body_contribution_corr;
    float sense_prpnsty, antisense_prpnsty;
        int i_in_nr;
        int sense_flag, antisense_flag, switch_flag;
        int t_traj_counter, t_term_freq_counter;
        
        
        
        
        // SIMULATION VARIABLES
        //RANDOM NUMBERS USED IN ALGORITHM
        double rtime,rxn;
        int reploss;
        int prox_term, histexc;
        
        //REACTION PROPENSITIES
        double totalprop;  // sum of propensities
        double *prpnsty = new double[nrxn]; // array of propensities (K4me +,- K27me +,-, Sense initiation, AS initiation)
        double *cmprpnsty = new double[nrxn]; // array of cumulative propensities
        
        
        
             
        
        // LOOP OVER NUMBER OF SIMULATIONS (nsims with initial uniform me0 and nsims with initial uniform me3)
        
        for(int j=1;j<=1*nsims;j++)
        {
            
            //INITIAL TIME AND SYSTEM STATE (INITIALIZE AND RECORD)
            t=0; // START TIME IS ALWAYS SET TO ZERO, ASSUMED TO BE JUST AFTER A REPLICATION EVENT
            reptime =  gsl_rng_uniform(r)*cycletime;  // SET TIME FOR FIRST REPLICATION EVENT: 1*cycletime for synchronised trajectories, otherwise multiply by gsl_rng_uniform(r)
            // SET COUNTERS FOR RECORDING TIMEPOINTS
            t_traj_counter = 0;
            t_term_freq_counter = 0;
            
            
            // Initialize transcription flags to 0
            sense_flag = 0;
            antisense_flag = 0;
            
            // Switch indicator set to zero at start of each simulation
            switch_flag = 0;
            
            // Initialize and record system state and time
            PK27me3_nr = 0;
            PK27me3_rem = 0;
            PK27me3_whole = 0;
            PK4me1 = 0;
            
            
            for(int i=0;i<N;i++)
            {
                cdomK27K4[i]= init_H3state;  //cdomK27K4[i]=3 for uniform K27me3 initial state; cdomK27K4[i]=-1 for uniform K4me1 initial state
                i_in_nr = ((nr_start <= i and i<= nr_end) or (lr_start <= i and i<= lr_end));
                
                PK27me3_nr += i_in_nr*(cdomK27K4[i]==3);
                PK27me3_rem += (1 - i_in_nr)*(cdomK27K4[i]==3);
                PK27me3_whole += (cdomK27K4[i]==3);
                PK4me1 += (cdomK27K4[i]==-1);
                
                
            }
            
            PK27me3_nr = PK27me3_nr/N_NR;
            PK27me3_rem = PK27me3_rem/(N-N_NR);
            PK27me3_whole = PK27me3_whole/N;
            PK4me1 = PK4me1/N;
            
            

            
            // STARTING INDIVIDUAL SIMULATION
            while(t<tmax and t>=0)
            {
                // COMPUTE CELL CYCLE PHASE DEPENDENT FACTOR THAT CONTROLS SPREADING (CONTROLS INTERACTIONS BETWEEN NR+LR AND BODY REGION)
                if(t<(reptime-(1-(0.5*spreading_phase))*cycletime) or t>(reptime-(0.5*spreading_phase)*cycletime))
                {
                    omega_t = 1;
                }
                else
                {
                    omega_t = 0.1;
                }
                
                // COMPUTE H3K27 METHYLATION AND H3K4 methylation PROPENSITIES
                
                // FIRST COMPUTE CONTRIBUTIONs TO FEEDBACK FROM HISTONES in different regions
                nr_contribution = 0;
                lr_contribution = 0;
                body_contribution = 0;
                
                for(int i=0;i<N;i++)
                {
                    if(nr_start<=i and i<=nr_end)
                    {
                        nr_contribution += (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3));
                        
                    }
                    else if(lr_start<=i and i<=lr_end)
                    {
                        lr_contribution += (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3));
                    }
                    else
                    {
                        body_contribution += (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3));
                    }
                    
                }
                // LOOP OVER HISTONES to compute propensities for adding H3K27me
                for(int i=0;i<N;i++)
                {
                    // IF INSIDE NR
                    if(nr_start<=i and i<=nr_end)
                    {
                        nr_contribution_corr = nr_contribution - (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3)); //Correction to ensure self-feedback not allowed
                        Ei = xi_nr_nr*nr_contribution_corr + xi_nr_lr*lr_contribution + omega_t*xi_nr_body*body_contribution;
                        prpnsty[i] = beta*((cdomK27K4[i]==0)*(gme01_nr + kme01*Ei) + (cdomK27K4[i]==1)*(gme12_nr + kme12*Ei) + (cdomK27K4[i]==2)*(gme23_nr + kme*Ei));
                        
                        
                    }
                    // ELSE IF INSIDE LR
                    else if(lr_start<=i and i<=lr_end)
                    {
                        lr_contribution_corr = lr_contribution - (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3)); //Correction to ensure self-feedback not allowed
                        Ei = xi_lr_nr*nr_contribution + xi_lr_lr*lr_contribution_corr + omega_t*xi_lr_body*body_contribution;
                        prpnsty[i] = beta*((cdomK27K4[i]==0)*(gme01_nr + kme01*Ei) + (cdomK27K4[i]==1)*(gme12_nr + kme12*Ei) + (cdomK27K4[i]==2)*(gme23_nr + kme*Ei));
                        
                    }
                    // ELSE IN BODY
                    else
                    {
                        body_contribution_corr = body_contribution - (rho_me2*(cdomK27K4[i]==2)+(cdomK27K4[i]==3)); //Correction to ensure self-feedback not allowed
                        Ei = omega_t*xi_body_nr*nr_contribution + omega_t*xi_body_lr*lr_contribution + xi_body_body*body_contribution_corr;
                        prpnsty[i] = beta*((cdomK27K4[i]==0)*(gme01 + kme01*Ei) + (cdomK27K4[i]==1)*(gme12 + kme12*Ei) + (cdomK27K4[i]==2)*(gme23 + kme*Ei));
                        
                    }
                    
                    
                    
                }
                
                
                
                
                // COMPUTE H3K27 background demethylation AND H3K4 background methylation AND demethylation PROPENSITIES
                
                for(int i=0;i<N;i++)
                {
                    
                    prpnsty[i+1*N]= gdem*(cdomK27K4[i]>0);  // Background removal of K27me
                    prpnsty[i+2*N]= gk4me*(cdomK27K4[i]==0); // Background addition of K4me1
                    prpnsty[i+3*N]= gk4dem*(cdomK27K4[i]==-1);   // Background removal of K4me1
                    
                }
                
                
                
                // COMPUTE TRANSCRIPTION PROPENSITIES
                PK27me3 = eta*PK27me3_nr + (1-eta)*PK27me3_rem;
                
                if(PK27me3<=PT27)
                {
                    sense_prpnsty = alpha_sense * (fmax - (PK27me3/PT27)*(fmax-fmin));//Propensity for sense transcription event; SET TO ZERO IF CURRENTLY IN A REFREACTORY PERIOD AFTER AN ANTISENSE TRANSCRIPTION EVENT
                    antisense_prpnsty = alpha_antisense * (fmax - (PK27me3/PT27)*(fmax-fmin)) * (1-sense_flag); //Propensity for antisense transcription event; SET TO ZERO IF CURRENTLY IN A REFREACTORY PERIOD AFTER AN SENSE TRANSCRIPTION EVENT
                }
                else
                {
                    sense_prpnsty = alpha_sense * fmin;
                    antisense_prpnsty = alpha_antisense * fmin;
                }
                
                prpnsty[4*N] = sense_prpnsty * (1-antisense_flag);
                prpnsty[4*N+1] = antisense_prpnsty * (1-sense_flag);
                
                // COMPUTE PROBABILITY OF PROXIMAL TERMINATION BASED ON CURRENT H3K4me coverage
                
                P_proximal = delta_proximal + fca_parameter*(1-PK4me1)*(1-delta_proximal-mu_distal);
                
                //CHECK IF TRAJECTORY HAS "SWITCHED"
                if((PK27me3 > chi*PT27) and switch_flag == 0) // USING OUR DEFINITION OF SWITCHING
                {
                    // IF SWITCHED THEN RECORD CURRENT TIME IN SWITCHOFF_TIME ARRAY
                    switchoff_time[j-1] = t/cycletime;  // SWITCH-OFF TIME FOR THE jTH TRAJECTORY

                    switch_flag = 1;
                }
                

                
                // COMPUTE TOTAL PROPENSITY AND  CUMULATIVE PROPENSITY
                totalprop=0;
                for(int k=0;k<nrxn;k++)
                {
                   
                    totalprop+=prpnsty[k];
                    cmprpnsty[k]=totalprop;
                }
                
                rtime=gsl_rng_uniform(r);
                rxn=gsl_rng_uniform(r);
                
                // COMPUTE NEXT REACTION TIME
                t_prev = t;  // STORE PREVIOUS REACTION TIME BEFORE UPDATING
                t-=log(rtime)/totalprop; // UPDATE TIME TO NEXT REACTION TIME

                // ACCOUNT FOR REFRACTORY PERIOD BEFORE PROCEEDING
                if(sense_flag==1 and t>=(t_prev+tau_S)) //CHECK IF REFREACTORY PERIOD AFTER A SENSE TRANSCRIPTION EVENT HAS ELAPSED
                {
                    //RESET SENSE TRANSCRIPTION FLAG
                    sense_flag = 0;
                    //Interrupt Gillespie algorithm to set current time to END OF CURRENT REFRACTORY PERIOD
                    t=t_prev+tau_S;
                    //RE-COMPUTE PROPENSITIES, ALLOWING NON-ZERO PROPENSITY FOR BOTH S AND AS TRANSCRIPTION EVENTS
                    totalprop += antisense_prpnsty;
                    cmprpnsty[nrxn-1] = totalprop;
                    //RE-START GILLESPIE ALGORITHM BY COMPUTING NEXT REACTION TIME
                    rtime=gsl_rng_uniform(r);
                    t-=log(rtime)/totalprop;
                }
                
                if(antisense_flag==1 and t>=(t_prev+tau_AS)) //CHECK IF REFREACTORY PERIOD AFTER AN ANTISENSE TRANSCRIPTION EVENT HAS ELAPSED
                {
                    //RESET SENSE TRANSCRIPTION FLAG
                    antisense_flag = 0;
                    //Interrupt Gillespie algorithm to set current time to END OF CURRENT REFRACTORY PERIOD
                    t=t_prev+tau_AS;
                    //RE-COMPUTE PROPENSITIES, ALLOWING NON-ZERO PROPENSITY FOR BOTH S AND AS TRANSCRIPTION EVENTS
                    totalprop += sense_prpnsty;
                    cmprpnsty[nrxn-2] = cmprpnsty[nrxn-3] + sense_prpnsty;
                    cmprpnsty[nrxn-1] = totalprop;
                    //RE-START GILLESPIE ALGORITHM BY COMPUTING NEXT REACTION TIME
                    rtime=gsl_rng_uniform(r);
                    t-=log(rtime)/totalprop;
                    
                }
                
                
                //CHECK IF SIMULATION SHOULD TERMINATE IMMEDIATELY
                if(t>=tmax)
                {
                    
                    //UPDATE SUMS (CORRESPONDING TO REMAINING RECORDING TIMEPOINTS)AND TERMINATE
                    while(t_traj_counter<t_traj_size)
                    {
                        
                        //UPDATE SUMS OVER TRAJECTORIES
                        sumovertraj_K27me3_nr[t_traj_counter] += PK27me3_nr;
                        sumovertraj_K27me3_whole[t_traj_counter] += PK27me3_whole;
                        sumovertraj_K4me1[t_traj_counter] += PK4me1;
                        sumovertraj_OFF_state[t_traj_counter]+= (PK27me3 > chi*PT27);
                        //Increment traj_counter to move to next recording timepoint
                        t_traj_counter++;
                    }
                    //UPDATE SPATIAL PROFILE
                    for(int i=0;i<N;i++)
                    {
                        sumovertraj_cdomK27[i] += (cdomK27K4[i]==3);
                        sumovertraj_cdomK4[i] += (cdomK27K4[i]==-1);
                    }
                    break;
                    
                    
                }
                
                //ELSE CHECK IF REPLICATION SHOULD OCCUR FIRST
                else if(t>reptime)
                {
                    //Interrupt Gillespie algorithm to set current time to reptime
                    t=reptime;
                    while(t>t_traj[t_traj_counter])
                    {
                        
                        //UPDATE SUMS OVER TRAJECTORIES
                        sumovertraj_K27me3_nr[t_traj_counter] += PK27me3_nr;
                        sumovertraj_K27me3_whole[t_traj_counter] += PK27me3_whole;
                        sumovertraj_K4me1[t_traj_counter] += PK4me1;
                        sumovertraj_OFF_state[t_traj_counter]+= (PK27me3 > chi*PT27);
                        //Increment t_traj_counter to move to next recording timepoint
                        t_traj_counter++;
                        
                    }
                    //CARRY OUT REPLICATION
                   for(int i=0;i<N;i+=2)
                    {
                        
                        reploss = (gsl_rng_uniform(r)<pinherit);//Set reploss=0 if gsl_rng_uniform(r)>pinherit
                        
                        if(cdomK27K4[i]>0 || cdomK27K4[i+1]>0) // Pairs of H3 histones inherited with probability 0.5 only if H3K27 carries methylation
                        {
                            cdomK27K4[i]=reploss*cdomK27K4[i];
                            cdomK27K4[i+1]=reploss*cdomK27K4[i+1];
                        }
                        else                                  // No inheritance if H3K27me is absent
                        {
                            cdomK27K4[i]=0;
                            cdomK27K4[i+1]=0;
                        }
                        
                        
                    }
                    //UPDATE reptime (TIME FOR NEXT REPLICATION EVENT)
                    reptime+=1*cycletime;
                }
                
                //ELSE UPDATE SUMS AND CARRY OUT NEXT REACTION
                else
                {
                    // UPDATE SUMS
                    // IF current reaction time t is less or equal to the current recording time, no update
                    
                    // Else update sums corresponding to current recording time, then move to next recording timepoint and repeat until current recording time crosses current reaction time
                    while(t>t_traj[t_traj_counter])
                    {
                        
                        //UPDATE SUMS OVER TRAJECTORIES
                        sumovertraj_K27me3_nr[t_traj_counter] += PK27me3_nr;
                        sumovertraj_K27me3_whole[t_traj_counter] += PK27me3_whole;
                        sumovertraj_K4me1[t_traj_counter] += PK4me1;
                        sumovertraj_OFF_state[t_traj_counter]+= (PK27me3 > chi*PT27);
                        //Increment t_traj_counter to move to next recording timepoint
                        t_traj_counter++;
                        
                    }
                    
                    // NEXT REACTION: METHYLATION OR DEMETHYLATION OR TRANSCRIPTION (DECIDE PROXIMAL OR DISTAL)
                    for(int k=0;k<=(nrxn-1);k++)
                    {
                        if(rxn<=(cmprpnsty[k]/totalprop))
                        {   //std::cout << "Reaction "<< k << "\n";
                            if(k<N)
                            {cdomK27K4[k]+=1;}  //Adding K27me
                            else if(k<2*N)
                            {cdomK27K4[k-N]-=1;} //Removing K27me
                            else if(k<3*N)
                            {cdomK27K4[k-2*N]-=1;} //Adding K4me
                            else if(k<4*N)
                            {cdomK27K4[k-3*N]+=1;} //Removing K4me
                            else
                            {
                                //FIRST MAKE TERMINATION CHOICE
                                prox_term = (gsl_rng_uniform(r) < P_proximal);
                                if(k==4*N) // IF SENSE TRANSCRIPTION EVENT
                                {
                                    // Set flag used to ensure mutual exclusivity of S and AS
                                    sense_flag = 1;
                                    if(prox_term == 1) // FOR A PROXIMAL SENSE TERMINATION EVENT
                                    {
                                        // RECORD SENSE PROXIMAL EVENT
                                        if(t<=t_term_freq[t_term_freq_counter])
                                        {
                                            sumovertraj_sense_prox_events[t_term_freq_counter]++;
                                        }
                                        else
                                        {
                                            while(t>t_term_freq[t_term_freq_counter])
                                            {
                                                t_term_freq_counter++;
                                            }
                                            sumovertraj_sense_prox_events[t_term_freq_counter]++;
                                        }
                                        
                                        //CARRY OUT REACTIONS ASSOCIATED WITH EVENT: REMOVAL OF H3K4ME1 ACROSS WHOLE LOCUS
                                        for(int i=0;i<N;i++)
                                        {
                                            cdomK27K4[i]+=(cdomK27K4[i]==-1)*(gsl_rng_uniform(r)<pdem_fld);
                                        }
                                        
                                    } // END OF PROXIMAL SENSE TERMINATION EVENT
                                    else   // FOR A DISTAL SENSE TERMINATION EVENT
                                    {
                                        // RECORD SENSE DISTAL EVENT
                                        if(t<=t_term_freq[t_term_freq_counter])
                                        {
                                            sumovertraj_sense_dist_events[t_term_freq_counter]++;
                                        }
                                        else
                                        {
                                            while(t>t_term_freq[t_term_freq_counter])
                                            {
                                                t_term_freq_counter++;
                                            }
                                            sumovertraj_sense_dist_events[t_term_freq_counter]++;
                                        }
                                        
                                        //CARRY OUT ALL REACTIONS ASSOCIATED WITH EVENT: REMOVAL OF H3K27ME, ADDITION OF H3K4ME, HISTONE EXCHANGE (ALL THREE REACTIONS CAN OCCUR THOUGHOUT THE LOCUS)
                                        for(int i=0;i<N;i++)
                                        {
                                            cdomK27K4[i]-=1*(cdomK27K4[i]>0)*(gsl_rng_uniform(r)<pdem); //REMOVAL OF H3K27ME3
                                            cdomK27K4[i]-=1*(cdomK27K4[i]==0)*(gsl_rng_uniform(r)<pk4me); //ADDITION OF H3K4ME1
                                            
                                            // HISTONE EXCHANGE
                                            if( (i%2) == 0)
                                            {
                                                histexc = (gsl_rng_uniform(r)>(2*pex)); //Set histexc=0 if gsl_rng_uniform(r)<(2*pex), i.e. if the histone octamer is lost
                                                cdomK27K4[i] = histexc*cdomK27K4[i];
                                                cdomK27K4[i+1] = histexc*cdomK27K4[i+1];
                                            }
                                        }
                                        
                                        
                                    } // END OF DISTAL SENSE TERMINATION EVENT
                                    
                                    
                                } // END OF SENSE TRANSCRIPTION EVENT
                                else   //IF ANTISENSE TRANSCRIPTION EVENT
                                {
                                    // Set flag used to ensure mutual exclusivity of S and AS
                                    antisense_flag = 1;
                                    if(prox_term == 1) // FOR A PROXIMAL SENSE TERMINATION EVENT
                                    {
                                        // RECORD ANTISENSE PROXIMAL EVENT
                                        if(t<=t_term_freq[t_term_freq_counter])
                                        {
                                            sumovertraj_antisense_prox_events[t_term_freq_counter]++;
                                        }
                                        else
                                        {
                                            while(t>t_term_freq[t_term_freq_counter])
                                            {
                                                t_term_freq_counter++;
                                            }
                                            sumovertraj_antisense_prox_events[t_term_freq_counter]++;
                                        }
                                        
                                        //CARRY OUT REACTIONS ASSOCIATED WITH EVENT: REMOVAL OF H3K4ME1 ACROSS WHOLE LOCUS
                                        for(int i=0;i<N;i++)
                                        {
                                            cdomK27K4[i]+=(cdomK27K4[i]==-1)*(gsl_rng_uniform(r)<pdem_fld);
                                        }
                                        
                                    } // END OF PROXIMAL ANTISENSE TERMINATION EVENT
                                    else   // FOR A DISTAL ANTISENSE TERMINATION EVENT
                                    {
                                        // RECORD ANTISENSE DISTAL EVENT
                                        if(t<=t_term_freq[t_term_freq_counter])
                                        {
                                            sumovertraj_antisense_dist_events[t_term_freq_counter]++;
                                        }
                                        else
                                        {
                                            while(t>t_term_freq[t_term_freq_counter])
                                            {
                                                t_term_freq_counter++;
                                            }
                                            sumovertraj_antisense_dist_events[t_term_freq_counter]++;
                                        }
                                        
                                        //CARRY OUT ALL REACTIONS ASSOCIATED WITH EVENT: REMOVAL OF H3K27ME, ADDITION OF H3K4ME, HISTONE EXCHANGE (ALL THREE REACTIONS CAN OCCUR THOUGHOUT THE LOCUS)
                                        for(int i=0;i<N;i++)
                                        {
                                            cdomK27K4[i]-=1*(cdomK27K4[i]>0)*(gsl_rng_uniform(r)<pdem); //REMOVAL OF H3K27ME3
                                            cdomK27K4[i]-=1*(cdomK27K4[i]==0)*(gsl_rng_uniform(r)<pk4me); //ADDITION OF H3K4ME1
                                            
                                            // HISTONE EXCHANGE
                                            if( (i%2) == 0)
                                            {
                                                histexc = (gsl_rng_uniform(r)>(2*pex)); //Set histexc=0 if gsl_rng_uniform(r)<(2*pex), i.e. if the histone octamer is lost
                                                cdomK27K4[i] = histexc*cdomK27K4[i];
                                                cdomK27K4[i+1] = histexc*cdomK27K4[i+1];
                                            }
                                        }
                                        
                                        
                                    } // END OF DISTAL ANTISENSE TERMINATION EVENT
                                }
                            }
                            
                            break;
                        }
                        
                        
                        
                    }
                    
                                       

                }
                // Compute and update current system state measures (after carrying out a reaction)
                PK27me3_nr = 0;
                PK27me3_rem = 0;
                PK27me3_whole = 0;
                PK4me1 = 0;
                
                for(int i=0;i<N;i++)
                {
                    i_in_nr = ((nr_start <= i and i<= nr_end) or (lr_start <= i and i<= lr_end));
                    PK27me3_nr += i_in_nr*(cdomK27K4[i]==3);
                    PK27me3_rem += (1 - i_in_nr)*(cdomK27K4[i]==3);
                    PK27me3_whole += (cdomK27K4[i]==3);
                    PK4me1 += (cdomK27K4[i]==-1);
                    
                    
                }
                PK27me3_nr = PK27me3_nr/N_NR;
                PK27me3_rem = PK27me3_rem/(N-N_NR);
                PK27me3_whole = PK27me3_whole/N;
                PK4me1 = PK4me1/N;
                
                
 
                
                
            }
            if(t<0)
            {
                std::cout << t << " totalprop " << totalprop  << "\n";
                break;
            }
        }
    
    std::cout << "Thread number " << tid << " completed job \n";
                //EACH THREAD WRITES TO CORRESPONDING OUTPUT FILES
                
                // OUTPUT TYPE 1a
                std::ofstream file1_;
                file1_.open(type1a_file_name);
                if(file1_.is_open())
                { for(int i=0;i<t_traj_size;i++)
                {
                    file1_ << t_traj[i]/cycletime << " " << sumovertraj_K27me3_nr[i]/nsims << " "<< sumovertraj_K27me3_whole[i]/nsims << " "<< sumovertraj_K4me1[i]/nsims <<  " " << sumovertraj_OFF_state[i] << "\n";
                    
                }
                    file1_.close();
                }
                
                // OUTPUT TYPE 1b
                std::ofstream file2_;
                file2_.open(type1b_file_name);
                if(file2_.is_open())
                { for(int i=0;i<t_term_freq_size;i++)
                {
                    file2_ << t_term_freq[i]/cycletime << " "<< sumovertraj_sense_prox_events[i]/nsims << " "<< sumovertraj_sense_dist_events[i]/nsims << " " << sumovertraj_antisense_prox_events[i]/nsims << " "<< sumovertraj_antisense_dist_events[i]/nsims<< "\n";
                    
                }
                    file2_.close();
                }
                
                // OUTPUT TYPE 1c
                std::ofstream file3_;
                file3_.open(type1c_file_name);
                if(file3_.is_open())
                {
                for(int i=0;i<nsims;i++)
                {
                    file3_ << switchoff_time[i] << "\n";

                    
                }
                    file3_.close();
                }
                
                // OUTPUT TYPE 1d
                std::ofstream file4_;
                file4_.open(type1d_file_name);
                if(file4_.is_open())
                { for(int i=0;i<N;i+=1)
                {
                    file4_ << i+1 << " " << (sumovertraj_cdomK27[i])/nsims << " " << (sumovertraj_cdomK4[i])/nsims <<"\n";
                    
                }
                    file4_.close();
                }
                
            
       
                
                //FREE MEMORY
                gsl_rng_free(r);
                delete [] cdomK27K4 ;
                delete [] sumovertraj_cdomK27 ;
                delete [] sumovertraj_cdomK4 ;
                delete [] prpnsty;
                delete [] cmprpnsty;
                delete [] sumovertraj_K27me3_nr;
                delete [] sumovertraj_K27me3_whole;
                delete [] sumovertraj_K4me1;
                delete [] sumovertraj_OFF_state;
                delete [] sumovertraj_sense_prox_events;
                delete [] sumovertraj_sense_dist_events;
                delete [] sumovertraj_antisense_prox_events;
                delete [] sumovertraj_antisense_dist_events;
                delete [] switchoff_time;
            

    
    delete [] t_traj;
    delete [] t_term_freq;
    
    nthreads = 1;
    end=clock();
    cpu_time_used = ((double) (end - start)) / (CLOCKS_PER_SEC);
    
    
    std::cout << "Time used:" << cpu_time_used << " " << "nthreads:" << nthreads << "\n";

    
    
        flag=0;
        return flag;
        
    }
    
