MODULE detour_dqmc_cfg_defs
USE dqmc_util
IMPLICIT NONE

  ! string parameters
  character, parameter :: COMMENT = "#"
  character, parameter :: SEPARAT = "="
  character, parameter :: COMMA   = ","

  ! string length
  integer, parameter :: slen = 60
  integer, parameter :: llen = 256
  integer, parameter :: alen = 10    ! array limit

  ! status param
  integer, parameter :: STAT_EOF     = -1
  integer, parameter :: STAT_COMMENT = 0
  integer, parameter :: STAT_NORMAL  = 1

  ! type def
  character(8), parameter :: TYPE_STR(3) =(/"Real*8  ","Integer ","Char(60)"/)

  integer, parameter :: TYPE_REAL    = 1 
  integer, parameter :: TYPE_INTEGER = 2 
  integer, parameter :: TYPE_STRING  = 3 

  ! default parameters
  integer, parameter :: N_Param = 40

  ! name of parameters
  character(len=*), parameter :: PARAM_NAME(N_Param) =  &
       &(/"HSF    ", &    ! methods of how HSF are generated                     
       &  "HSFin  ", &    ! File name of HSF input, if HSF = HSF_FROM_FILE       
       &  "HSFout ", &    ! File name of HSF output                              
       &  "HSFtype", &    ! HSF type, 0 = discrete, 1 = continuous
       &  "L      ", &    ! Number of time slides                                
       &  "U      ", &    ! Parameter for potential energy                       
       &  "accept ", &    ! accept counts        
       &  "bcond  ", &    ! boundary conditions
       &  "debug  ", &    ! flag for debug information output                    
       &  "delta1 ", &    ! parameter for continous HSF
       &  "delta2 ", &    ! parameter for continous HSF
       &  "difflim", &    ! limit of tolerable difference                        
       &  "dmu    ", &    ! perturbation of mu                                   
       &  "dtau   ", &    ! discritize parameter                                 
       &  "errrate", &    ! tolerable error                                      
       &  "fixwrap", &    ! fix nwrap to input value
       &  "gamma  ", &    ! correction of Metropolis ratio                       
       &  "gfile  ", &    ! geometry definition                                  
       &  "mu_dn  ", &    ! parameter for chemical potential                     
       &  "mu_up  ", &    ! parameter for chemical potential                     
       &  "n      ", &    ! number of particles                                  
       &  "nbin   ", &    ! number of statisitical bins                          
       &  "nhist  ", &    ! print history or not
       &  "nitvl  ", &    ! number of interval, for continuous FT integration    
       &  "north  ", &    ! frequence of orthogonalization                       
       &  "npass  ", &    ! number of measurement sweeps                         
       &  "ntry   ", &    ! number of global moves per sweeps                    
       &  "nwarm  ", &    ! number of warmup sweeps                              
       &  "nwrap  ", &    ! frequence of recomputing H                           
       &  "nx     ", &    ! number of sites in x direction                       
       &  "ny     ", &    ! number of sites in y direction                       
       &  "nz     ", &    ! number of sites in y direction                       
       &  "ofile  ", &    ! prefix of output files                               
       &  "reject ", &    ! rejection counts                                     
       &  "seed   ", &    ! random seed                                          
       &  "ssxx   ", &    ! use iterative refinement during sweep
       &  "t_dn   ", &    ! parameter for kinetic energy                         
       &  "t_up   ", &    ! parameter for kinetic energy                         
       &  "tausk  ", &    ! frequence of unequal time measurement                
       &  "tdm    "/)     ! compute time dependent measurement

  ! default values
  character(len=*), parameter :: PARAM_DVAL(N_Param) =  &
       &(/"-1      ", &    ! HSF    
       &  "HSF.in  ", &    ! HSFin  
       &  "HSF.out ", &    ! HSFout 
       &  "0       ", &    ! HSFtype
       &  "12      ", &    ! L      
       &  "0.0     ", &    ! U      
       &  "0       ", &    ! accept 
       &  "0,0,0   ", &    ! bcond
       &  "0       ", &    ! debug  
       &  "1.0     ", &    ! delta1
       &  "1.0     ", &    ! delta2
       &  "0.001   ", &    ! difflim
       &  "0.0     ", &    ! dmu    
       &  "0.125   ", &    ! dtau   
       &  "0.001   ", &    ! errrate
       &  "0       ", &    ! fixwrap
       &  "0.0     ", &    ! gamma  
       &  "geom.def", &    ! gfile  
       &  "0.0     ", &    ! mu_up  
       &  "0.0     ", &    ! mu_dn  
       &  "16      ", &    ! n      
       &  "10      ", &    ! nbin   
       &  "0       ", &    ! nhist  
       &  "4       ", &    ! nitvl  
       &  "12      ", &    ! north  
       &  "5000    ", &    ! npass  
       &  "0       ", &    ! ntry   
       &  "1000    ", &    ! nwarm  
       &  "12      ", &    ! nwrap  
       &  "4       ", &    ! nx     
       &  "4       ", &    ! ny     
       &  "2       ", &    ! nz     
       &  "quest   ", &    ! ofile  
       &  "0       ", &    ! reject 
       &  "0       ", &    ! seed 
       &  "0       ", &    ! ssxx
       &  "1.0     ", &    ! t_up      
       &  "1.0     ", &    ! t_dn      
       &  "10      ", &    ! tausk  
       &  "0       "/)     ! tdm
  
  ! parameter type
  integer, parameter :: PARAM_TYPE(N_Param) = &
       &(/TYPE_INTEGER, &    ! HSF    
       &  TYPE_STRING,  &    ! HSFin  
       &  TYPE_STRING,  &    ! HSFout 
       &  TYPE_INTEGER, &    ! HSFtype
       &  TYPE_INTEGER, &    ! L      
       &  TYPE_REAL,    &    ! U      
       &  TYPE_INTEGER, &    ! accept 
       &  TYPE_REAL   , &    ! bcond
       &  TYPE_INTEGER, &    ! debug  
       &  TYPE_REAL,    &    ! delta1
       &  TYPE_REAL,    &    ! delta2
       &  TYPE_REAL,    &    ! difflim
       &  TYPE_REAL,    &    ! dmu    
       &  TYPE_REAL,    &    ! dtau   
       &  TYPE_REAL,    &    ! errrate
       &  TYPE_INTEGER, &    ! fixwrap
       &  TYPE_REAL,    &    ! gamma  
       &  TYPE_STRING,  &    ! gfile  
       &  TYPE_REAL,    &    ! mu_up     
       &  TYPE_REAL,    &    ! mu_dn     
       &  TYPE_INTEGER, &    ! n      
       &  TYPE_INTEGER, &    ! nbin   
       &  TYPE_INTEGER, &    ! nhist  
       &  TYPE_INTEGER, &    ! nitvl  
       &  TYPE_INTEGER, &    ! north  
       &  TYPE_INTEGER, &    ! npass  
       &  TYPE_INTEGER, &    ! ntry   
       &  TYPE_INTEGER, &    ! nwarm  
       &  TYPE_INTEGER, &    ! nwrap  
       &  TYPE_INTEGER, &    ! nx     
       &  TYPE_INTEGER, &    ! ny     
       &  TYPE_INTEGER, &    ! nz     
       &  TYPE_STRING,  &    ! ofile  
       &  TYPE_INTEGER, &    ! reject 
       &  TYPE_INTEGER, &    ! seed   
       &  TYPE_INTEGER, &    ! ssxx
       &  TYPE_REAL,    &    ! t_up      
       &  TYPE_REAL,    &    ! t_dn      
       &  TYPE_INTEGER, &    ! tausk  
       &  TYPE_INTEGER/)     ! tdm

  ! is array parameter
  logical, parameter :: PARAM_ARRAY(N_Param) = &
       &(/.false.,&           ! HSF    
       &  .false.,&           ! HSFin  
       &  .false.,&           ! HSFout 
       &  .false.,&           ! HSFtype
       &  .false.,&           ! L      
       &  .true. ,&           ! U      
       &  .false.,&           ! accept 
       &  .true. ,&           ! bcond
       &  .false.,&           ! debug  
       &  .false.,&           ! delta1
       &  .false.,&           ! delta2
       &  .false.,&           ! difflim
       &  .false.,&           ! dmu    
       &  .false.,&           ! dtau   
       &  .false.,&           ! errrate
       &  .false.,&           ! fixwrap
       &  .false.,&           ! gamma  
       &  .false.,&           ! gfile  
       &  .true. ,&           ! mu_up     
       &  .true. ,&           ! mu_dn     
       &  .false.,&           ! n      
       &  .false.,&           ! nbin   
       &  .false.,&           ! nhist  
       &  .false.,&           ! nitvl  
       &  .false.,&           ! north  
       &  .false.,&           ! npass  
       &  .false.,&           ! ntry   
       &  .false.,&           ! nwarm  
       &  .false.,&           ! nwrap  
       &  .false.,&           ! nx     
       &  .false.,&           ! ny     
       &  .false.,&           ! nz     
       &  .false.,&           ! ofile  
       &  .false.,&           ! reject 
       &  .false.,&           ! seed   
       &  .false.,&           ! ssxx      
       &  .true. ,&           ! t_up
       &  .true. ,&           ! t_dn
       &  .false.,&           ! tausk  
       &  .false./)           ! tdm

  !
  ! Data Type
  ! =========
  !
  type Param
     integer        :: id           ! hashcode from name
     character(len=slen)  :: pname  ! parameter name
     integer        :: ptype        ! type of parameter
     logical        :: isArray      ! is the param an array?
     logical        :: isSet        ! is the parameter been set?
     character(len=slen)  :: defaultval   ! default value
     type(Param), pointer :: next 

     ! values
     integer        :: ival
     real(wp)       :: rval
     integer, pointer     :: iptr(:) 
     real(wp), pointer    :: rptr(:) 
  end type Param

  type config
     type(Param), pointer :: record(:)    ! head of the linked list
     integer        :: nParam
     logical        :: hasDef = .false.
  end type config

END MODULE
