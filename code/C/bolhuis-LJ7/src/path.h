/*------------------INLINE SUBSTITUTIONS-------------------------------------*/
#define TRUE 1
#define FALSE 0
#define CANCELLED -99

#define P_INDEX(ptr)  (ptr - pts)
#define F_COMP(a,b) (fabs( a-b ) < 1e-10)
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define MESSAGE(a) printf("Message:"#a"\n")
#define SIGN(a)   ( 2*(a>0) -1) 

#define PI           3.141592653589793
#define SQRT2        1.4142136
#define BIGNUM       1e99
#define EPS          1e-10
#define NVT          1

#define dprint(expr) printf(#expr " = %d\n",expr)
#define gprint(expr) printf(#expr " = %g\n",expr)
#define vprint(expr) printf(#expr " = ( %16.14g %16.14g %16.14g ) \n",expr.x, expr.y, expr.z)
#define vector_inp(a, b) (a.x * b.x + a.y * b.y + a.z * b.z )
#define vector_cross(a,b,h)   h.x = a.y * b.z - a.z * b.y; h.y = a.z * b.x - a.x * b.z; h.z = a.x * b.y - a.y * b.x;
#define vector_times(a,b,h)   h.x = a.x * b.x; h.y = a.y * b.y; h.z = a.z * b.z;
#define vector_divide(a,b,h)   h.x = a.x / b.x; h.y = a.y / b.y; h.z = a.z / b.z;
#define vector_plustimes(a,b,h)   h.x += a.x * b.x; h.y += a.y * b.y; h.z += a.z * b.z;
#define vector_mintimes(a,b,h)   h.x -= a.x * b.x; h.y -= a.y * b.y; h.z -= a.z * b.z;
#define scalar_times(a,b,h)  h.x = a.x *b; h.y =a.y * b; h.z = a.z * b;  
#define scalar_divide(a,b,h)  h.x = a.x /b; h.y =a.y / b; h.z = a.z / b;  
#define scalar_plustimes(a,b,h)  h.x += a.x *b; h.y +=a.y * b; h.z += a.z * b;  
#define scalar_mintimes(a,b,h)  h.x -= a.x *b; h.y -=a.y * b; h.z -= a.z * b;  
#define vector_add(a,b,h)  {h.x = a.x +b.x; h.y =a.y + b.y; h.z = a.z + b.z;  }
#define vector_minus(a,b,h) { h.x = a.x - b.x; h.y =a.y - b.y; h.z = a.z - b.z;  }
#define update_average(a,b) {a.now =b; a.sum += a.now; a.sumsq += a.now*a.now; a.n++;}
#define SWAP(a,b,c) {c=a;a=b;b=c;}

#define POTENERGY sys.block_stats.aver[0]
#define KIN_ENER  sys.block_stats.aver[1]
#define TOT_ENER  sys.block_stats.aver[2]
#define PX_MOM    sys.block_stats.aver[3]
#define PY_MOM    sys.block_stats.aver[4]
#define PZ_MOM    sys.block_stats.aver[5]
#define LX_MOM    sys.block_stats.aver[6]
#define LY_MOM    sys.block_stats.aver[7]
#define LZ_MOM    sys.block_stats.aver[8]
#define TEMP      sys.block_stats.aver[9]
#define TEMP_L    sys.block_stats.aver[10]
#define HAHB      sys.block_stats.aver[11]
#define HB        sys.block_stats.aver[12]
#define BLIP       sys.block_stats.aver[13]


/*#define FORWARD_ACC      sys.block_stats.mcacc[0]
#define BACKWARD_ACC     sys.block_stats.mcacc[1]
#define REPTLEFT_ACC     sys.block_stats.mcacc[2]
#define REPTRIGHT_ACC    sys.block_stats.mcacc[3]
#define DIFLEFT_ACC      sys.block_stats.mcacc[4]
#define DIFRIGHT_ACC     sys.block_stats.mcacc[5]
#define INREGIONB        sys.block_stats.mcacc[6]
*/


#define NSTAT  15
#define NACC    160
#define NSHOOTAV 30
#define MBIN  1000
#define MAXSTATES  30
#define MAXRUN 50
#define MAXPATH 10000
#define MAXWINDOW 50
#define MAXREPLICA 20
#define LATBIN 10
#define STATINTERVAL 40
#define NPART 20

/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct vector_type {

      double         x,y,z;

} vector;


typedef struct tensor_type {

      vector         x,y,z;

} tensor;



typedef struct particle_type *Pts_ptr;

typedef struct particle_type {

      vector        r                 ,
                    v                 ,
                    f                 ;
	
      double        dr2;

} Pts;

typedef struct slice_type {

       Pts          pts[NPART]      ;

       int          ha                ,
                    hb                ,
                    empty             ;

       double       xi,
                    s,
	            state_rc[MAXSTATES],
                    ordered_dr2[NPART];

} Slice;



typedef struct replica_type {

  //  Slice         slice[MAXPATH];
               
  
  int           pathlen,
                index, 
                swapindex,
                ntotal,
                avlen,
                navlen,
                type;  

  double           lambda, 
                   dos;
  //                   window[4][2];

  vector        string;

  FILE          *pathfilep[5];

} Replica;

typedef struct aver_type {

      char          name[100]              ;

      int           n                  ;

      double        now                ,
                    sum                ,
                    sumsq              ;

} Average;

typedef struct mcacc_type {

      char          name[100]              ;

      int           acc,try             ;

      double        ratio              ;

} Mcacc;



typedef struct stats_type {

      Average       aver[NSTAT];

      Average       shoot[NSHOOTAV];       


      Mcacc         mcacc[NACC];

  

} Stats;

typedef struct corr_type {
  
      double        flux0[MAXSTATES]       ,
                    flux1[MAXSTATES]       ,
                    free_en[MBIN]         ,
	            energy[MBIN]          ,
                    energy2[MBIN]         ,
                    kin_ener[MBIN]        ,
                    kin_ener2[MBIN]       ,
                    tot_ener[MBIN]        ,
	            tot_ener2[MBIN]       ,
	            crosshistAB[MAXSTATES][MAXREPLICA][MBIN],
	            crosshistBA[MAXSTATES][MAXREPLICA][MBIN],
     	            rate[MAXSTATES][MAXSTATES],
                    currentmat[LATBIN][LATBIN][2];       

       
      int           navb[MBIN]             ,
            	    route[MAXREPLICA] ,
            	    turnovers[MAXREPLICA] ,
	            tothist[MAXREPLICA][MBIN],
	            crosshistAB_L[MAXREPLICA][MBIN],
	            crosshistBA_L[MAXREPLICA][MBIN],
	            crosshistAB_H[MAXREPLICA][MBIN],
	            crosshistBA_H[MAXREPLICA][MBIN],
	            crosshistAB_block[MAXREPLICA][MBIN],
	            crosshistBA_block[MAXREPLICA][MBIN],
	            crosshistAB_L_block[MAXREPLICA][MBIN],
	            crosshistBA_L_block[MAXREPLICA][MBIN],
	            crosshistAB_H_block[MAXREPLICA][MBIN],
	            crosshistBA_H_block[MAXREPLICA][MBIN],
	            fehistAB[MAXSTATES][MAXREPLICA][MBIN],
	            fehistBA[MAXSTATES][MAXREPLICA][MBIN],
	            boundhistAB[MAXREPLICA][MBIN],
	            boundhistBA[MAXREPLICA][MBIN],
	            densmat[LATBIN][LATBIN],
	            ncurrentmat[LATBIN][LATBIN],
            	    pathtype[NSHOOTAV][4],
            	    type_mat[MAXSTATES][MAXREPLICA][4],
     	            mstis_mat[MAXSTATES][MAXSTATES],
                    nflux0[MAXSTATES]       ,
                    nflux1[MAXSTATES]       ,
      //	            swaprepind[2000000],
      //	            stateind[2000000],
	            fcount                ,
                    n                     ;
      double        free_ch               ;
} Scr;



typedef struct bivar_type {
      double        sig_r              ,
                    sig_v              ,
                    sigma              ,
                    c_rv               ,
                    s12os11            ,
                    sqrts11            ,
                    sqrtSos11          ,
                    norm               ,
                    norder             ,
                    pfmax              ,
                    pbmax              ,
                    pb[10]             ,
                    pf[10]             ;


} Bivar;

typedef struct langevin_type {
  
  
     double         gamma              ,
                    gamma_a            ,
                    gamma_b            ,
                    mu                 ,
	            c0                 ,
	            c1                 ,
	            c2                 ,
	            c3                 ,
	            c1min2             ;

} Langevin;

typedef struct state_type {

  vector pos;
  double min, 
         min2;
  double lambda[MAXREPLICA];
  
  
  Replica srep[MAXREPLICA]; ;

  Slice target;

} State;


typedef struct status_type {
  
      int           move               ,
                    itime              ,
                    accepted           ,
                    hAlast             ,
                    hBfirst            , 
	            critL1              ,
                    critB1              ,
                    critS1              ,
                    crit01              ,	
                    nbar               ;
 
      double        critL              ,
                    critB              ,
                    critS              ,
                    crit0              ,	          
                    en0                ;
} Status;

typedef struct path_typ {
      int           nslices            ,
                    ninter             ,
                    nneigh1            ,
                    nneigh2            ,
                    tagged1            ,
                    tagged2            ,
	            updated            ,
	            int_failure        ,	
     	            ntarget            ,
	            hbfirst            ,
                    quench             ,
	            endpoint_constrained,
                    which_array        ;

      double        energy             ,
                    action             ,
                    kin_ener           ,
	            tot_ener           ,
 	            hAx                ,
 	            hAy                ,
 	            hAr                ,
 	            hBx                ,
 	            hBy                ,
 	            hBr                ,
                    order_par          ,
                    ring_dist2         ,
	            final_width        ,
   	            divide             ,
                    temp               ,
	            targetdif          ,
                    upper_moi_1        ,
                    lower_moi_1        ,
                    upper_moi_2        ,
                    lower_moi_2        ,
  	            crit[5]            ,
                    hAmax[5]           ,
                    hBmin[5]           ;

      vector        region             ,
	            initial            ,
                    final              ;

      
      Status        status              ;

} Path;  





typedef struct parameter_type {

      int           ncycle1            ,
                    ncycle2            ,
                    npart              ,
	            nshoot             ,
                    nreverse           ,
                    nrept              ,
	            ndif               ,
	            nswap              ,
	            nswapstates        ,
                    nacc               ,
                    nfix               ,
                    maxtrial           ,
                    maxlength          ,
	            nreplica           ,
	            current_replica    ,
                    block_no           ,

	            initial_path       ,
                    sim_type           ,
                    start_type         ,
                    ensemble_type      ,
       	            freq_check         ,
                    freq_graphics      ,
                    freq_config        ,
	            freq_corr          ,
	            freq_dump          ,
                    graphics           ,
      //	            ninter             ,
	            hAmin[5]           ,
	            hAmax[5]           ,
	            hBmin[5]           ,
	            hBmax[5]           ,
	            ncrit              , 
	            nstates            ,
	            initial_state       ,
	            vers               ;

      double        volume             ,
                    energy             ,
                    fixed_energy       ,
                    beta               ,
                    gamma              ,
                    delta_t            ,
                    temp               ,
	            sqrttemp           ,
                    dt                 ,
                    bDdt               ,
                    maxdist            ,
                    scalefactor        ,
                    potw               ,
	            saddlex            ,
	            min_stable         ,
                    saddley            ,
 	            dvmax              ,
	            sigma              , 
                    trial              ,
	            initx              ,
	            inity              ,
  	            lambda[MAXREPLICA]  ,
                    **norm             ;

  //      vector        states;

      vector        boxl               ;


      Stats         final_stats        ,
                    block_stats        ;

      
      FILE          *filep             ;

} Sys;

                    
/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED VARIABLES-------------------------------*/

extern Replica       trial_replica,*replica[MAXREPLICA];
extern Slice         *slice,*trial,*fwtrial,*bwtrial;
extern Sys           sys;
extern Path          path;
extern Langevin      langevin;
extern Stats         nulstat;
extern vector        nulvec;
extern Scr           scr,nulscr;
extern Bivar         bivar              ;
extern State         state[MAXSTATES];


/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/






extern void  set_up             ( void );
extern void  confinput          ( void );
extern void  confoutput         ( void );
extern void  distoutput         ( void );
extern void  targetoutput        ( void );
extern double  energy           ( Slice * );
extern void check_energy      ( Slice *) ;
extern int   langevin_inertia    (Slice *);
extern double total_energy            ( Slice * );
extern void  calculate_forces   ( Slice * );
extern void propagate(Slice *);
extern void propagate_NHL(Slice *);
extern void find_minimum(Slice *);


extern double ran3              ( void );
extern double gssran            ( void );
extern void   gssbivar          (double *, double *);
extern void quicksort(double *, int , int );
extern void  plot_conf          ( void );
extern void  init_graphics      ( void );
extern void  drawline( vector,vector);
extern int trajectory(Slice *, int);
extern int trajectory_state_i(Slice *, Replica *, int, int);
extern int analyse_state_i(Slice *, Replica *, int, int);
extern int in_state(Slice *);
extern int in_upper_window(Slice *, double, int );

extern int shoot(Replica *);
extern int swap_replica(int,int);
extern int swap_replica_0(int,int);
extern int swap_replica_n(int);
extern int swap_states(int ,int );
extern int reverse_replica(int);
extern int add_newstate(Slice *);

extern double create_all_rc(Slice *);
extern double get_rc(Slice *,int);
extern double print_rc(Slice *,int);
extern int analyse ( Slice *, Replica *, int , int );
extern void analyse_full ( Slice *, Replica *, int , int );

extern int initial_path();
extern void reset_center(Slice *);

extern void  choose_new_momenta( Slice *);
extern void fix_energy(Slice *);


extern void Init_Graphics(int , char *[],Slice *);
extern void mainloop_for_graphics();
