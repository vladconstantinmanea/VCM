C******************************
C Laplace Tidal Equations, with Rayleigh Dissipation (linear), with ice cap
C
C As discussed in "Solving the Laplace Tidal Equations using Freely
C    Available, Easily Extensible Finite Element Software," Sewell and
C    Manea, Computers and Geosciences, v???, pp??.-
C
C If run without changes, this program produces the plots seen in
C    Figure 6 of that article, but users can easily modify it to solve
C    many other tidal dissipation problems (with ice cap).  Parameters 
C    most likely to be modified by users are preceded by C???????????
C
C Runs with the free PDE2D demo downloadable at www.pde2d.com
C 
C The default output consists of the following plots (in PS file format): 
C U,V, ETA and Q at 0%,25%,50%,75%,100% for the last period, and two plots
C with the geographical distribution of tidal dissipation.  In these 
C plots, P1=longitude, P2=latitude.
C
C******************************
C
C     **************************                                              
C     * PDE2D 9.6 MAIN PROGRAM *                                              
C     **************************                                              
C     *** 2D PROBLEM SOLVED (COLLOCATION METHOD) ***                          
C##############################################################################
C     Is double precision mode to be used?  Double precision is recommended.  #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If double precision mode is used, variables and functions assigned   +#
C     + names beginning with a letter in the range A-H or O-Z will be DOUBLE +#
C     + PRECISION, and you should use double precision constants and FORTRAN +#
C     + expressions throughout; otherwise such variables and functions will  +#
C     + be of type REAL.  In either case, variables and functions assigned   +#
C     + names beginning with I,J,K,L,M or N will be of INTEGER type.         +#
C     +                                                                      +#
C     + It is possible to convert a single precision PDE2D program to double +#
C     + precision after it has been created, using an editor.  Just change   +#
C     + all occurrences of "real" to "double precision"                      +#
C     +                    " tdp" to "dtdp"  (note leading blank)            +#
C     + Any user-written code or routines must be converted "by hand", of    +#
C     + course.  To convert from double to single, reverse the changes.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      implicit double precision (a-h,o-z)
      parameter (neqnmx=  99)
C##############################################################################
C     NP1GRID = number of P1-grid lines                                       #
C##############################################################################
      PARAMETER (NP1GRID = 37)        
C##############################################################################
C     NP2GRID = number of P2-grid lines                                       #
C##############################################################################
      PARAMETER (NP2GRID = 19)        
C##############################################################################
C     How many differential equations (NEQN) are there in your problem?       #
C##############################################################################
      PARAMETER (NEQN = 4)            
      parameter (np3grid = 1)
C        DIMENSIONS OF WORK ARRAYS                                            
C        SET TO 1 FOR AUTOMATIC ALLOCATION                                    
      PARAMETER (IRWK8Z=           1)
      PARAMETER (IIWK8Z=           1)
      PARAMETER (NXP8Z=101,NYP8Z=101,KDEG8Z=1,NZP8Z=KDEG8Z+1)
C##############################################################################
C     The solution is normally saved on an NP1+1 by NP2+1 rectangular grid    #
C     of points,                                                              #
C                    P1 = P1A + I*(P1B-P1A)/NP1,    I = 0,...,NP1             #
C                    P2 = P2A + J*(P2B-P2A)/NP2,    J = 0,...,NP2             #
C     Enter values for NP1 and NP2.  Suggested values: NP1=NP2=25.            #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you want to save the solution at an arbitrary user-specified set  +#
C     + of points, set NP2=0 and NP1+1=number of points.  In this case you   +#
C     + can request tabular output, but no plots can be made.                +#
C     +                                                                      +#
C     + If you set NEAR8Z=1 in the main program, the values saved at each    +#
C     + output point will actually be the solution as evaluated at a nearby  +#
C     + collocation point.  For most problems this obviously will produce    +#
C     + less accurate output or plots, but for certain (rare) problems, a    +#
C     + solution component may be much less noisy when plotted only at       +#
C     + collocation points.                                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      PARAMETER (NP1 = 90)            
      PARAMETER (NP2 = 45)            
C##############################################################################
C     The solution will be saved (for possible postprocessing) at the NSAVE+1 #
C     time points                                                             #
C                          T0 + K*(TF-T0)/NSAVE                               #
C     K=0,...,NSAVE.  Enter a value for NSAVE.                                #
C                                                                             #
C     If a user-specified constant time step is used, NSTEPS must be an       #
C     integer multiple of NSAVE.                                              #
C##############################################################################
      PARAMETER (NSAVE = 200)         
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e     
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu   
      dimension p1grid(np1grid),p2grid(np2grid),p3grid(np3grid),p1out8z(       
     &0:np1,0:np2),p2out8z(0:np1,0:np2),p3out8z(0:np1,0:np2),p1cross(100       
     &),p2cross(100),tout8z(0:nsave)                                           
C      dimension xres8z(nxp8z),yres8z(nyp8z),zres8z(nzp8z),                    
C     & ures8z(neqn,nxp8z,nyp8z,nzp8z)                                         
      allocatable iwrk8z(:),rwrk8z(:)                                          
C      dimension iwrk8z(iiwk8z),rwrk8z(irwk8z)                                 
      character*40 title                                                       
      logical linear,crankn,noupdt,nodist,fillin,evcmpx,adapt,plot,lsqfi       
     &t,fdiff,solid,econ8z,ncon8z,restrt,gridid                                
      common/dtdp14/ sint8z(20),bint8z(20),slim8z(20),blim8z(20)               
      common/dtdp15/ evlr8z,ev0r,evli8z,ev0i,evcmpx                            
      common/dtdp16/ p8z,evr8z(50),evi8z(50)                                   
      common/dtdp19/ toler(neqnmx),adapt                                       
      common/dtdp30/ econ8z,ncon8z                                             
      common/dtdp45/ perdc(neqnmx)                                             
      common/dtdp46/ eps8z,cgtl8z,npmx8z,itype,near8z                          
      common/dtdp52/ nxa8z,nya8z,nza8z,kd8z                                    
      common/dtdp53/ work8z(nxp8z*nyp8z*nzp8z+9)                               
      common/dtdp64/ amin8z(4*neqnmx),amax8z(4*neqnmx)                         
      common/dtdp76/ mdim8z,nx18z,ny18z,p1a,p1b,p2a,p2b,uout(0:np1,0:np2       
     &,4,neqn,0:nsave)                                                         
      pi = 4.0*atan(1.d0)                                                      
      zr8z = 0.0                                                               
      nxa8z = nxp8z                                                            
      nya8z = nyp8z                                                            
      nza8z = nzp8z                                                            
      nx18z = np1+1                                                            
      ny18z = np2+1                                                            
      mdim8z = 4                                                               
      kd8z = kdeg8z                                                            
C##############################################################################
C     If you don't want to read the FINE PRINT, default NPROB.                #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you want to solve several similar problems in the same run, set   +#
C     + NPROB equal to the number of problems you want to solve.  Then NPROB +#
C     + loops through the main program will be done, with IPROB=1,...,NPROB, +#
C     + and you can make the problem parameters vary with IPROB.  NPROB      +#
C     + defaults to 1.                                                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C???????????
C   Work through 10*NPROB periods.  The first run, delete restart file 
C      "pde2d.res" before starting.  If you need to rerun (sometimes it takes
C      many periods for convergence) it will start from the end of the last 
C      run, if you do not delete "pde2d.res".
      NPROB = 2                                                               
      do 78755 iprob=1,nprob                                                   
C##############################################################################
C     PDE2D solves the time-dependent system (note: U,F,G,U0 may be vectors,  #
C     C,RHO may be matrices):                                                 #
C                                                                             #
C        C(X,Y,T,U,Ux,Uy)*d(U)/dT = F(X,Y,T,U,Ux,Uy,Uxx,Uyy,Uxy)              #
C                                                                             #
C     or the steady-state system:                                             #
C                                                                             #
C        F(X,Y,U,Ux,Uy,Uxx,Uyy,Uxy) = 0                                       #
C                                                                             #
C     or the linear and homogeneous eigenvalue system:                        #
C                                                                             #
C        F(X,Y,U,Ux,Uy,Uxx,Uyy,Uxy) = lambda*RHO(X,Y)*U                       #
C                                                                             #
C     with boundary conditions:                                               #
C                                                                             #
C                  G(X,Y,[T],U,Ux,Uy) = 0                                     #
C           (periodic boundary conditions are also permitted)                 #
C                                                                             #
C     For time-dependent problems there are also initial conditions:          #
C                                                                             #
C            U = U0(X,Y)   at T=T0                                            #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If your PDEs involve the solution at points other than (P1,P2), the  +#
C     + function                                                             +#
C     +             (D)OLDSOL2(IDER,IEQ,PP1,PP2,KDEG)                        +#
C     + will interpolate (using interpolation of degree KDEG=1,2 or 3) to    +#
C     + (PP1,PP2) the function saved in UOUT(*,*,IDER,IEQ,ISET) on the last  +#
C     + time step or iteration (ISET) for which it has been saved.  Thus,    +#
C     + for example, if IDER=1, this will return the latest value of         +#
C     + component IEQ of the solution at (PP1,PP2), assuming this has not    +#
C     + been modified using UPRINT... If your equations involve integrals of +#
C     + the solution, for example, you can use (D)OLDSOL2 to approximate     +#
C     + these using the solution from the last time step or iteration.       +#
C     +                                                                      +#
C     + CAUTION: For a steady-state or eigenvalue problem, you must reset    +#
C     + NOUT=1 if you want to save the solution each iteration.              +#
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
C     + A system of NEQN complex partial differential equations must be      +#
C     + written as a system of 2*NEQN real equations, by separating the      +#
C     + equations into their real and imaginary parts.  However, note that   +#
C     + the complex arithmetic abilities of FORTRAN can be used to simplify  +#
C     + this separation.  For example, the complex PDE:                      +#
C     +      I*(Uxx+Uyy) - 1/(1+U**10) = 0,   where U = UR + UI*I            +#
C     + would be difficult to split up analytically, but using FORTRAN       +#
C     + expressions it is easy:                                              +#
C     +   F1 = -(UIxx+UIyy) -  REAL(1.0/(1.0+CMPLX(UR,UI)**10))              +#
C     +   F2 =  (URxx+URyy) - AIMAG(1.0/(1.0+CMPLX(UR,UI)**10))              +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C     You may now define global parameters, which may be referenced in any    #
C     of the "FORTRAN expressions" you input throughout the rest of this      #
C     interactive session.  You will be prompted alternately for parameter    #
C     names and their values; enter a blank name when you are finished.       #
C                                                                             #
C     Parameter names are valid FORTRAN variable names, starting in           #
C     column 1.  Thus each name consists of 1 to 6 alphanumeric characters,   #
C     the first of which must be a letter.  If the first letter is in the     #
C     range I-N, the parameter must be an integer.                            #
C                                                                             #
C     Parameter values are either FORTRAN constants or FORTRAN expressions    #
C     involving only constants and global parameters defined on earlier       #
C     lines.  They may also be functions of the problem number IPROB, if      #
C     you are solving several similar problems in one run (NPROB > 1).  Note  #
C     that you are defining global CONSTANTS, not functions; e.g., parameter  #
C     values may not reference any of the independent or dependent variables  #
C     of your problem.                                                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you define other parameters here later, using an editor, you must +#
C     + add them to COMMON block /PARM8Z/ everywhere this block appears, if  +#
C     + they are to be "global" parameters.                                  +#
C     +                                                                      +#
C     + The variable PI is already included as a global parameter, with an   +#
C     + accurate value 3.14159...                                            +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C???????????
C  *** Enceladus data  ***
C  revolution frequency (radians/sec)
      Omega = 5.31 D-5                                                   
      PERIOD = 2.*pi/Omega
C  ocean density (kg/m^3)
      sigma = 1000.0 d0
C  satellite radius (m)
      RT = 0.252 D6                                                       
C  ocean depth (m)
      h0 = 10000.d0 
C  eccentricity of orbit
      e = 0.0047D0
C  obliquity 
      theta0 = 0.0           
C  surface graviational acceleration (m/sec^2)
      g = 0.11 d0
C  mode = 1 for eccentricity tidal heating
C  mode = 2 for obliquity tidal heating
      mode = 1 
C  Rayleigh dissipation coefficient (1/sec)
      alpha = 1.d-7
C  Youngs modulus for ice (Pa)
      EM     =                                                                 
     & 8.778 D9                                                                 
C  Ice cap thickness (m)
      hi     =                                                                 
     & 1000. 
C  Poisson ratio for ice
      Rnu    =                                                                 
     & 0.33                                                                    
C  bottom drag coefficient 
C  CAUTION: If Cd > 0, you must change LINEAR and NOUPDT to .FALSE.
C           and FDIFF to .TRUE.   You may optionally change ADAPT.
      Cd = 0. d0 
C
      betinv =                                                                 
     & EM*hi/RT**2                                                          
C##############################################################################
C     A collocation finite element method is used, with bi-cubic Hermite      #
C     basis functions on the elements (small rectangles) defined by the grid  #
C     points:                                                                 #
C               P1GRID(1),...,P1GRID(NP1GRID)                                 #
C               P2GRID(1),...,P2GRID(NP2GRID)                                 #
C     You will first be prompted for NP1GRID, the number of P1-grid points,   #
C     then for P1GRID(1),...,P1GRID(NP1GRID).  Any points defaulted will be   #
C     uniformly spaced between the points you define; the first and last      #
C     points cannot be defaulted.  Then you will be prompted similarly        #
C     for the number and values of the P2-grid points.  The limits on the     #
C     parameters are then:                                                    #
C               P1GRID(1) < P1 < P1GRID(NP1GRID)                              #
C               P2GRID(1) < P2 < P2GRID(NP2GRID)                              #
C                                                                             #
C##############################################################################
      call dtdpwx(p1grid,np1grid,0)                                            
      call dtdpwx(p2grid,np2grid,0)                                            
C        P1GRID DEFINED                                                        
      P1GRID(1) =                                                              
     & 0                                                                       
      P1GRID(NP1GRID) =                                                        
     & 2*pi                                                                    
C        P2GRID DEFINED                                                        
      P2GRID(1) =                                                              
     & -pi/2.                                                                  
      P2GRID(NP2GRID) =                                                        
     & pi/2.                                                                   
C                                                                              
      p3grid(1) = 0                                                            
      call dtdpwx(p1grid,np1grid,1)                                            
      call dtdpwx(p2grid,np2grid,1)                                            
C##############################################################################
C     If you don't want to read the FINE PRINT, enter ISOLVE = 1.             #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The following linear system solvers are available:                   +#
C     +                                                                      +#
C     + 1. Sparse direct method                                              +#
C     +               Harwell Library routine MA27 (used by permission) is   +#
C     +               used to solve the (positive definite) "normal"         +#
C     +               equations A**T*A*x = A**T*b.  The normal equations,    +#
C     +               which are essentially the equations which would result +#
C     +               if a least squares finite element method were used     +#
C     +               instead of a collocation method, are substantially     +#
C     +               more ill-conditioned than the original system Ax = b,  +#
C     +               so it may be important to use high precision if this   +#
C     +               option is chosen.                                      +#
C     + 2. Frontal method                                                    +#
C     +               This is an out-of-core band solver.  If you want to    +#
C     +               override the default number of rows in the buffer (11),+#
C     +               set a new value for NPMX8Z in the main program.        +#
C     + 3. Jacobi conjugate gradient iterative method                        +#
C     +               A preconditioned conjugate gradient iterative method   +#
C     +               is used to solve the (positive definite) normal        +#
C     +               equations.  High precision is also important if this   +#
C     +               option is chosen.  (This solver is MPI-enhanced, if    +#
C     +               MPI is available.)  If you want to override the        +#
C     +               default convergence tolerance, set a new relative      +#
C     +               tolerance CGTL8Z in the main program.                  +#
C     + 4. Local solver (normal equations)                                   +#
C     + 5. Local solver (original equations)                                 +#
C     +               Choose these options ONLY if alterative linear system  +#
C     +               solvers have been installed locally.  See subroutines  +#
C     +               (D)TD3M, (D)TD3N in file (d)subs.f for instructions    +#
C     +               on how to add local solvers.                           +#
C     + 6. MPI-based parallel band solver                                    +#
C     +               This is a parallel solver which runs efficiently on    +#
C     +               multiple processor machines, under MPI.  It is a       +#
C     +               band solver, with the matrix distributed over the      +#
C     +               available processors.  Choose this option ONLY if the  +#
C     +               solver has been activated locally.  See subroutine     +#
C     +               (D)TD3O in file (d)subs.f for instructions on how to   +#
C     +               activate this solver and the MPI-enhancements to the   +#
C     +               conjugate gradient solver.                             +#
C     +                                                                      +#
C     + Enter ISOLVE = 1,2,3,4,5 or 6 to select a linear system solver.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     If you don't want to read the FINE PRINT, enter ISOLVE = 1.             #
C##############################################################################
      ISOLVE =          1                                                      
C        *******TIME-DEPENDENT PROBLEM                                         
      itype = 2                                                                
C##############################################################################
C     Enter the initial time value (T0) and the final time value (TF), for    #
C     this time-dependent problem.  T0 defaults to 0.                         #
C                                                                             #
C     TF is not required to be greater than T0.                               #
C##############################################################################
      T0 = 0.0                                                                 
      T0 =                                                                     
     & 0                                                                       
      TF =                                                                     
     & 10*PERIOD                                                               
C##############################################################################
C     Is this a linear problem? ("linear" means all differential equations    #
C     and all boundary conditions are linear).  If you aren't sure, it is     #
C     safer to answer "no".                                                   #
C##############################################################################
      LINEAR = .TRUE.                                                          
C##############################################################################
C     Do you want the time step to be chosen adaptively?  If you answer       #
C     'yes', you will then be prompted to enter a value for TOLER(1), the     #
C     local relative time discretization error tolerance.  The default is     #
C     TOLER(1)=0.01.  If you answer 'no', a user-specified constant time step #
C     will be used.  We suggest that you answer 'yes' and default TOLER(1)    #
C     (although for certain linear problems, a constant time step may be much #
C     more efficient).                                                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If a negative value is specified for TOLER(1), then ABS(TOLER(1)) is +#
C     + taken to be the "absolute" error tolerance.  If a system of PDEs is  +#
C     + solved, by default the error tolerance specified in TOLER(1) applies +#
C     + to all variables, but the error tolerance for the J-th variable can  +#
C     + be set individually by specifying a value for TOLER(J) using an      +#
C     + editor, after the end of the interactive session.                    +#
C     +                                                                      +#
C     + Each time step, two steps of size dt/2 are taken, and that solution  +#
C     + is compared with the result when one step of size dt is taken.  If   +#
C     + the maximum difference between the two answers is less than the      +#
C     + tolerance (for each variable), the time step dt is accepted (and the +#
C     + next step dt is doubled, if the agreement is "too" good); otherwise  +#
C     + dt is halved and the process is repeated.  Note that forcing the     +#
C     + local (one-step) error to be less than the tolerance does not        +#
C     + guarantee that the global (cumulative) error is less than that value.+#
C     + However, as the tolerance is decreased, the global error should      +#
C     + decrease correspondingly.                                            +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ADAPT = .FALSE.                                                          
      IF (IPROB.eq.1) then
         TOLER(1) = 1.D20
      ELSE
         TOLER(1) = 0.05                                                      
      ENDIF
C##############################################################################
C     If you don't want to read the FINE PRINT, it is safe (though possibly   #
C     very inefficient) to enter 'no'.                                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If your time-dependent problem is linear with all PDE and boundary   +#
C     + condition coefficients independent of time except inhomogeneous      +#
C     + terms, then a large savings in execution time may be possible if     +#
C     + this is recognized (the LU decomposition computed on the first step  +#
C     + can be used on subsequent steps).  Is this the case for your         +#
C     + problem?  (Caution: if you answer 'yes' when you should not, you     +#
C     + will get incorrect results with no warning.)                         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NOUPDT = .TRUE.                                                          
C##############################################################################
C     The time stepsize will be constant, DT = (TF-T0)/NSTEPS.  Enter a       #
C     value for NSTEPS, the number of time steps.                             #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you later turn on adaptive step control, the time stepsize will be+#
C     + chosen adaptively, between an upper limit of DTMAX = (TF-T0)/NSTEPS  +#
C     + and a lower limit of 0.0001*DTMAX.                                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C???????????
C   Number of time steps per 10-period run (must be multiple of 200) 
      NSTEPS =                                                                 
     & 1600                                                                     
      dt = (tf-t0)/max(nsteps,1)                                               
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Is the Crank-Nicolson scheme to be used to discretize time?  If you  +#
C     + answer 'no', a backward Euler scheme will be used.                   +#
C     +                                                                      +#
C     + If a user-specified constant time step is chosen, the second order   +#
C     + Crank Nicolson method is recommended only for problems with very     +#
C     + well-behaved solutions, and the first order backward Euler scheme    +#
C     + should be used for more difficult problems.  In particular, do not   +#
C     + use the Crank Nicolson method if the left hand side of any PDE is    +#
C     + zero, for example, if a mixed elliptic/parabolic problem is solved.  +#
C     +                                                                      +#
C     + If adaptive time step control is chosen, however, an extrapolation   +#
C     + is done between the 1-step and 2-step answers which makes the Euler  +#
C     + method second order, and the Crank-Nicolson method strongly stable.  +#
C     + Thus in this case, both methods have second order accuracy, and both +#
C     + are strongly stable.                                                 +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      CRANKN = .FALSE.                                                         
      FDIFF = .FALSE.                                                          
C##############################################################################
C     You may calculate one or more integrals (over the entire region) of     #
C     some functions of the solution and its derivatives.  How many integrals #
C     (NINT), if any, do you want to calculate?                               #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + In the FORTRAN program created by the preprocessor, the computed     +#
C     + values of the integrals will be returned in the vector SINT8Z.  If   +#
C     + several iterations or time steps are done, only the last computed    +#
C     + values are saved in SINT8Z (all values are printed).                 +#
C     +                                                                      +#
C     + A limiting value, SLIM8Z(I), for the I-th integral can be set        +#
C     + below in the main program.  The computations will then stop          +#
C     + gracefully whenever SINT8Z(I) > SLIM8Z(I), for any I=1...NINT.       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NINT =          1                                                        
C##############################################################################
C     You may calculate one or more boundary integrals (over the entire       #
C     boundary) of some functions of the solution and its derivatives.  How   #
C     many boundary integrals (NBINT), if any, do you want to calculate?      #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + In the FORTRAN program created by the preprocessor, the computed     +#
C     + values of the integrals will be returned in the vector BINT8Z.  If   +#
C     + several iterations or time steps are done, only the last computed    +#
C     + values are saved in BINT8Z (all values are printed).                 +#
C     +                                                                      +#
C     + A limiting value, BLIM8Z(I), for the I-th boundary integral can be   +#
C     + set below in the main program.  The computations will then stop      +#
C     + gracefully whenever BINT8Z(I) > BLIM8Z(I), for any I=1...NBINT.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NBINT =          0                                                       
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Normally, interpolation is done to approximate the initial values    +#
C     + using cubic Hermites.  Since some derivatives must be interpolated,  +#
C     + if the initial values are not smooth (ie, have large or infinite     +#
C     + derivatives), the resulting cubic interpolants may have undesired    +#
C     + noise or large spikes.  Do you want to compute a least squares       +#
C     + approximation to the initial values, rather than an interpolant?     +#
C     + The least squares fit is generally much smoother, but requires one   +#
C     + extra linear system solution.                                        +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      LSQFIT = .FALSE.                                                         
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to read the initial conditions from the restart file,    +#
C     + if it exists (and use the conditions supplied above if it does not   +#
C     + exist)?                                                              +#
C     +                                                                      +#
C     + If so, PDE2D will dump the final solution at the end of each run     +#
C     + into a restart file "pde2d.res".  Thus the usual procedure for       +#
C     + using this dump/restart option is to make sure there is no restart   +#
C     + file in your directory left over from a previous job, then the       +#
C     + first time you run this job, the initial conditions supplied above   +#
C     + will be used, but on the second and subsequent runs the restart file +#
C     + from the previous run will be used to define the initial conditions. +#
C     +                                                                      +#
C     + You can do all the "runs" in one program, by setting NPROB > 1.      +#
C     + Each pass through the DO loop, T0,TF,NSTEPS and possibly other       +#
C     + parameters may be varied, by making them functions of IPROB.         +#
C     +                                                                      +#
C     + If the 2D or 3D collocation method is used, the coordinate           +#
C     + transformation should not change between dump and restart.           +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      RESTRT = .TRUE.                                                         
C     GRIDID = .FALSE. IF FINITE ELEMENT GRID CHANGES BETWEEN DUMP, RESTART    
      GRIDID = .TRUE.                                                          
C##############################################################################
C     If you do not have any periodic boundary conditions, enter IPERDC=0.    #
C                                                                             #
C     Enter IPERDC=1 for periodic conditions at P1 = P1GRID(1),P1GRID(NP1GRID)#
C           IPERDC=2 for periodic conditions at P2 = P2GRID(1),P2GRID(NP2GRID)#
C           IPERDC=4 for periodic conditions on both P1 and P2                #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + When periodic boundary conditions are selected, they apply to all    +#
C     + variables by default.  To turn off periodic boundary conditions on   +#
C     + the I-th variable, set PERDC(I) to 0 (or another appropriate value   +#
C     + of IPERDC) below in the main program and set the desired boundary    +#
C     + conditions in subroutine GB8Z, "by hand".                            +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      IPERDC =          1                                                      
C##############################################################################
C     The solution is saved on an NP1+1 by NP2+1 rectangular grid covering    #
C     the rectangle (P1A,P1B) x (P2A,P2B).  Enter values for P1A,P1B,P2A,P2B. #
C     These variables are usually defaulted.                                  #
C                                                                             #
C     The defaults are P1A = P1GRID(1), P1B = P1GRID(NP1GRID)                 #
C                      P2A = P2GRID(1), P2B = P2GRID(NP2GRID)                 #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Solution derivatives may occasionally be calculated incorrectly at   +#
C     + the boundary, if the transformation X(P1,P2),Y(P1,P2) has a singular +#
C     + Jacobian there, for example, at R=0 when polar coordinates are used. +#
C     + In this case it is suggested you move P1A, P1B, P2A or P2B in very   +#
C     + slightly from the boundary.                                          +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C        defaults for p1a,p1b,p2a,p2b                                          
      p1a = p1grid(1)                                                          
      p1b = p1grid(np1grid)                                                    
      p2a = p2grid(1)                                                          
      p2b = p2grid(np2grid)                                                    
C        DEFINE P1A,P1B,P2A,P2B IMMEDIATELY BELOW:                             
      call dtdpx3(np1,np2,0,p1a,p1b,p2a,p2b,zr8z,zr8z,hp18z,hp28z,hp38z,       
     &p1out8z,p2out8z,p3out8z,npts8z)                                          
      call dtdpqx(np1grid,np2grid,np3grid,isolve,neqn,ii8z,ir8z,iperdc)        
      if (iiwk8z.gt.1) ii8z = iiwk8z                                           
      if (irwk8z.gt.1) ir8z = irwk8z                                           
C        *******allocate workspace                                             
      allocate (iwrk8z(ii8z),rwrk8z(ir8z))                                     
C        *******DRAW GRID LINES?                                               
      PLOT = .FALSE.                                                            
C        *******call pde solver                                                
      call dtdp3x(p1grid, p2grid, p3grid, np1grid,np2grid, -1, neqn, p1o       
     &ut8z, p2out8z, p3out8z, uout, tout8z, npts8z, t0, dt, nsteps, nout       
     &, nsave, crankn, noupdt, itype, linear, isolve, rwrk8z, ir8z, iwrk       
     &8z, ii8z, iperdc, plot, lsqfit, fdiff, nint, nbint, restrt, gridid       
     &)                                                                        
      deallocate (iwrk8z,rwrk8z)                                               
C        *******read from restart file to array ures8z                         
C      call dtdpr3(1,xres8z,nxp8z,yres8z,nyp8z,zres8z,nzp8z,ures8z,neqn)       
C        *******write array ures8z back to restart file                        
C      call dtdpr3(2,xres8z,nxp8z,yres8z,nyp8z,zres8z,nzp8z,ures8z,neqn)       
C        *******call user-written postprocessor                                
      call postpr(tout8z,nsave,p1out8z,p2out8z,np1,np2,uout,neqn,
     & iprob,nprob)      
      if (iprob.lt.nprob) go to 78755
C        *******SURFACE PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       V                                                    #
C                5       Vx                                                   #
C                6       Vy                                                   #
C                7       ETA                                                  #
C                8       ETAx                                                 #
C                9       ETAy                                                 #
C               10       Q                                                    #
C               11       Qx                                                   #
C               12       Qy                                                   #
C                .        .                                                   #
C                .        .                                                   #
C##############################################################################
      IVAR =          1                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
      ISET1 =                                                                  
     & 180                                                                     
      ISET2 =                                                                  
     & 200                                                                     
      ISINC =                                                                  
     & 5                                                                       
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0.0                                                               
      UMAX = 0.0                                                               
      UMIN =                                                                   
     & alow                                                                    
      UMAX =                                                                   
     & ahigh                                                                   
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'U (m/sec)                               '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78756 is8z=iset1,iset2,isinc                                          
      call dtdplo(p1out8z,p2out8z,p3out8z,uout(0,0,ivara8z,ivarb8z,is8z)       
     &,np1,np2,0,3,ix8z,jy8z,0,title,vlon,vlat,umin,umax,tout8z(is8z))         
78756 continue                                                                 
C        *******SURFACE PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       V                                                    #
C                5       Vx                                                   #
C                6       Vy                                                   #
C                7       ETA                                                  #
C                8       ETAx                                                 #
C                9       ETAy                                                 #
C               10       Q                                                    #
C               11       Qx                                                   #
C               12       Qy                                                   #
C                .        .                                                   #
C                .        .                                                   #
C##############################################################################
      IVAR =          4                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
      ISET1 =                                                                  
     & 180                                                                     
      ISET2 =                                                                  
     & 200                                                                     
      ISINC =                                                                  
     & 5                                                                       
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0.0                                                               
      UMAX = 0.0                                                               
      UMIN =                                                                   
     & alow                                                                    
      UMAX =                                                                   
     & ahigh                                                                   
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'V (m/sec)                               '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78757 is8z=iset1,iset2,isinc                                          
      call dtdplo(p1out8z,p2out8z,p3out8z,uout(0,0,ivara8z,ivarb8z,is8z)       
     &,np1,np2,0,3,ix8z,jy8z,0,title,vlon,vlat,umin,umax,tout8z(is8z))         
78757 continue                                                                 
C        *******SURFACE PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       V                                                    #
C                5       Vx                                                   #
C                6       Vy                                                   #
C                7       ETA                                                  #
C                8       ETAx                                                 #
C                9       ETAy                                                 #
C               10       Q                                                    #
C               11       Qx                                                   #
C               12       Qy                                                   #
C                .        .                                                   #
C                .        .                                                   #
C##############################################################################
      IVAR =          7                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
      ISET1 =                                                                  
     & 180                                                                     
      ISET2 =                                                                  
     & 200                                                                     
      ISINC =                                                                  
     & 5                                                                       
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0.0                                                               
      UMAX = 0.0                                                               
      UMIN =                                                                   
     & alow                                                                    
      UMAX =                                                                   
     & ahigh                                                                   
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'ETA (m)                                 '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78758 is8z=iset1,iset2,isinc                                          
      call dtdplo(p1out8z,p2out8z,p3out8z,uout(0,0,ivara8z,ivarb8z,is8z)       
     &,np1,np2,0,3,ix8z,jy8z,0,title,vlon,vlat,umin,umax,tout8z(is8z))         
78758 continue                                                                 
C        *******SURFACE PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       V                                                    #
C                5       Vx                                                   #
C                6       Vy                                                   #
C                7       ETA                                                  #
C                8       ETAx                                                 #
C                9       ETAy                                                 #
C               10       Q                                                    #
C               11       Qx                                                   #
C               12       Qy                                                   #
C                .        .                                                   #
C                .        .                                                   #
C##############################################################################
      IVAR =         10                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
      ISET1 =                                                                  
     & 180                                                                     
      ISET2 =                                                                  
     & 200                                                                     
      ISINC =                                                                  
     & 5                                                                       
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0.0                                                               
      UMAX = 0.0                                                               
      UMIN =                                                                   
     & alow                                                                    
      UMAX =                                                                   
     & ahigh                                                                   
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Q (Pa)                                  '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78759 is8z=iset1,iset2,isinc                                          
      call dtdplo(p1out8z,p2out8z,p3out8z,uout(0,0,ivara8z,ivarb8z,is8z)       
     &,np1,np2,0,3,ix8z,jy8z,0,title,vlon,vlat,umin,umax,tout8z(is8z))         
78759 continue                                                                 
C        *******SURFACE PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       V                                                    #
C                5       Vx                                                   #
C                6       Vy                                                   #
C                7       ETA                                                  #
C                8       ETAx                                                 #
C                9       ETAy                                                 #
C##############################################################################
      IVAR =          6                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = NSAVE                                                            
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0                                                            
      UMAX = 0                                                            
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Energy dissipation rate (W/m^2)'                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78779 is8z=iset1,iset2,isinc                                          
      call dtdplo(p1out8z,p2out8z,p3out8z,uout(0,0,ivara8z,ivarb8z,is8z)       
     &,np1,np2,0,3,ix8z,jy8z,0,title,vlon,vlat,umin,umax,tout8z(is8z))         
78779 continue                                                                 
C        *******CONTOUR PLOTS                                                  
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means U  (possibly as modified by UPRINT,..)               #
C                2       Ux                                                   #
C                3       Uy                                                   #
C                4       W                                                    #
C                5       Wx                                                   #
C                6       Wy                                                   #
C                7       H                                                    #
C                8       Hx                                                   #
C                9       Hy                                                   #
C##############################################################################
      IVAR =          6                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
      ISET1 = NSAVE                                                           
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to scale the axes on the plot so that the region is      +#
C     + undistorted?  Otherwise the axes will be scaled so that the figure   +#
C     + approximately fills the plot space.                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NODIST = .FALSE.                                                          
C                                                                              
      ivar8z = 4*(ivarb8z-1)+ivara8z                                           
      alow = amin8z(ivar8z)                                                    
      ahigh = amax8z(ivar8z)                                                   
C##############################################################################
C     Enter lower (UMIN) and upper (UMAX) bounds for the contour values. UMIN #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     Labeled contours will be drawn corresponding to the values              #
C                                                                             #
C                  UMIN + S*(UMAX-UMIN),    for S=0.05,0.15,...0.95.          #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, UMIN and UMAX are set to the minimum and maximum values  +#
C     + of the variable to be plotted.  For a common scaling, you may want   +#
C     + to set UMIN=ALOW, UMAX=AHIGH.  ALOW and AHIGH are the minimum and    +#
C     + maximum values over all output points and over all saved time steps  +#
C     + or iterations.                                                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = 0.0                                                               
      UMAX = 0.0                                                               
C##############################################################################
C     Do you want two additional unlabeled contours to be drawn between each  #
C     pair of labeled contours?                                               #
C##############################################################################
      FILLIN = .FALSE.                                                         
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Energy dissipation rate (W/m^2)'                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78780 is8z=iset1,iset2,isinc                                          
      call dtdpln(uout(0,0,ivara8z,ivarb8z,is8z),np1,np2,0,p1a,p1b,p2a,p       
     &2b,zr8z,zr8z,3,ix8z,jy8z,0,title,umin,umax,nodist,fillin,tout8z(is       
     &8z),zr8z,zr8z,zr8z,zr8z,2,ical8z)                                        
78780 continue                                                                 
78755 continue                                                                 
      call endgks                                                              
      stop                                                                     
      end                                                                      
                                                                               
                                                                               
      subroutine tran8z(itrans,p1,p2,p38z)                                     
      implicit double precision (a-h,o-z)                                      
      common /dtdp41/x,y,z8z,x1,x2,x38z,y1,y2,y38z,z18z,z28z,z38z,x11,x2       
     &1,x31,x12,x22,x32,x13,x23,x33,y11,y21,y31,y12,y22,y32,y13,y23,y33,       
     &z11,z21,z31,z12,z22,z32,z13,z23,z33                                      
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                     
C##############################################################################
C     You can solve problems in your region only if you can describe it by    #
C                         X = X(P1,P2)                                        #
C                         Y = Y(P1,P2)                                        #
C     with constant limits on the parameters P1,P2.  If your region is        #
C     rectangular, enter ITRANS=0 and the trivial parameterization            #
C                         X = P1                                              #
C                         Y = P2                                              #
C     will be used.  Otherwise, you need to read the FINE PRINT below.        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If P1,P2 represent polar or other non-Cartesian coordinates, you can +#
C     + reference the Cartesian coordinates X,Y and derivatives of your      +#
C     + unknowns with respect to these coordinates, when you define your     +#
C     + PDE coefficients, boundary conditions, and volume and boundary       +#
C     + integrals, if you enter ITRANS .NE. 0.  Enter:                       +#
C     +   ITRANS = 1, if P1,P2 are polar coordinates, that is, if            +#
C     +               P1=R, P2=Theta, where    X = R*cos(Theta)              +#
C     +                                        Y = R*sin(Theta)              +#
C     +   ITRANS = -1, same as ITRANS=1, but P1=Theta, P2=R                  +#
C     +   ITRANS = 3, to define your own coordinate transformation.  In      +#
C     +               this case, you will be prompted to define X,Y and      +#
C     +               their first and second derivatives in terms of P1,P2.  +#
C     +               Because of symmetry, you will not be prompted for all  +#
C     +               of the second derivatives.  If you make a mistake in   +#
C     +               computing any of these derivatives, PDE2D will usually +#
C     +               be able to issue a warning message. (X1 = dX/dP1, etc) +#
C     +   ITRANS = -3, same as ITRANS=3, but you will only be prompted to    +#
C     +               define X,Y; their first and second derivatives will    +#
C     +               be approximated using finite differences.              +#
C     +   When ITRANS = -3 or 3, the first derivatives of X,Y must all be    +#
C     +   continuous.                                                        +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ITRANS =          0                                                      
C                                                                              
      z8z = p38z                                                               
      z38z = 1                                                                 
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine pdes8z(yd8z,i8z,j8z,kint8z,p1,p2,p38z,t,uu8z)                 
      implicit double precision (a-h,o-z)                                      
      parameter (neqnmx=  99)                                                  
C        un8z(1,I),un8z(2,I),... hold the (rarely used) values                 
C        of UI,UI1,... from the previous iteration or time step                
      common /dtdp5x/un8z(10,neqnmx)                                           
      common /dtdp18/norm1,norm2,n38z                                          
      double precision norm1,norm2,n38z,normx,normy,nz8z                       
      dimension uu8z(10,neqnmx)                                                
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                         
      zr8z = 0.0                                                               
      U  = uu8z(1, 1)                                                          
      U1 = uu8z(2, 1)                                                          
      U2 = uu8z(3, 1)                                                          
      U11= uu8z(5, 1)                                                          
      U22= uu8z(6, 1)                                                          
      U12= uu8z(8, 1)                                                          
      U21= uu8z(8, 1)                                                          
      V  = uu8z(1, 2)                                                          
      V1 = uu8z(2, 2)                                                          
      V2 = uu8z(3, 2)                                                          
      V11= uu8z(5, 2)                                                          
      V22= uu8z(6, 2)                                                          
      V12= uu8z(8, 2)                                                          
      V21= uu8z(8, 2)                                                          
      ETA  = uu8z(1, 3)                                                        
      ETA1 = uu8z(2, 3)                                                        
      ETA2 = uu8z(3, 3)                                                        
      ETA11= uu8z(5, 3)                                                        
      ETA22= uu8z(6, 3)                                                        
      ETA12= uu8z(8, 3)                                                        
      ETA21= uu8z(8, 3)                                                        
      Q  = uu8z(1, 4)                                                          
      Q1 = uu8z(2, 4)                                                          
      Q2 = uu8z(3, 4)                                                          
      Q11= uu8z(5, 4)                                                          
      Q22= uu8z(6, 4)                                                          
      Q12= uu8z(8, 4)                                                          
      Q21= uu8z(8, 4)                                                          
      call dtdpcd(p1,p2,p38z)                                                  
      call dtdpcb(p1,p2,p38z,norm1,norm2,n38z,x,y,z8z,normx,normy,nz8z,3       
     &)                                                                        
      call dtdpcc(p1,p2,p38z,                                                  
     &          U1,U2,zr8z,U11,U22,zr8z,U12,zr8z,zr8z,                         
     & x,y,z8z,Ux,Uy,uz8z,Uxx,Uyy,uzz8z,Uxy,uxz8z,uyz8z,                       
     & Uyx,uzx8z,uzy8z,dvol,darea)                                             
      Unorm = Ux*normx + Uy*normy                                              
      call dtdpcc(p1,p2,p38z,                                                  
     &          V1,V2,zr8z,V11,V22,zr8z,V12,zr8z,zr8z,                         
     & x,y,z8z,Vx,Vy,uz8z,Vxx,Vyy,uzz8z,Vxy,uxz8z,uyz8z,                       
     & Vyx,uzx8z,uzy8z,dvol,darea)                                             
      Vnorm = Vx*normx + Vy*normy                                              
      call dtdpcc(p1,p2,p38z,                                                  
     &          ETA1,ETA2,zr8z,ETA11,ETA22,zr8z,ETA12,zr8z,zr8z,               
     & x,y,z8z,ETAx,ETAy,uz8z,ETAxx,ETAyy,uzz8z,ETAxy,uxz8z,uyz8z,             
     & ETAyx,uzx8z,uzy8z,dvol,darea)                                           
      ETAnorm = ETAx*normx + ETAy*normy                                        
      call dtdpcc(p1,p2,p38z,                                                  
     &          Q1,Q2,zr8z,Q11,Q22,zr8z,Q12,zr8z,zr8z,                         
     & x,y,z8z,Qx,Qy,uz8z,Qxx,Qyy,uzz8z,Qxy,uxz8z,uyz8z,                       
     & Qyx,uzx8z,uzy8z,dvol,darea)                                             
      Qnorm = Qx*normx + Qy*normy                                              
                          if (i8z.eq.0) then                                   
      yd8z = 0.0                                                               
C##############################################################################
C     Enter FORTRAN expressions for the functions whose integrals are to be   #
C     calculated and printed.  They may be functions of                       #
C                                                                             #
C        X,Y,U,Ux,Uy,Uxx,Uyy,Uxy                                              #
C            V,Vx,Vy,Vxx,Vyy,Vxy                                              #
C            ETA,ETAx,ETAy,ETAxx,ETAyy,ETAxy                                  #
C            Q,Qx,Qy,Qxx,Qyy,Qxy and (if applicable) T                        #
C             .  .   .   .   .    .    .                                      #
C                                                                             #
C     The parameters P1,P2 and derivatives with respect to these may also     #
C     be referenced (U1 = dU/dP1, etc):                                       #
C           U1,U2,U11,U22,U12                                                 #
C           V1,V2,V11,V22,V12                                                 #
C           ETA1,ETA2,ETA11,ETA22,ETA12                                       #
C           Q1,Q2,Q11,Q22,Q12                                                 #
C             .    .    .    .     .                                          #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you only want to integrate a function over part of the region,    +#
C     + define that function to be zero in the rest of the region.           +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C                                                  INTEGRAL DEFINED            
      if (kint8z.eq.    1) yd8z =
     & Cd*sigma*(u**2+v**2)**1.5*RT**2*cos(y)
     & + alpha*sigma*h0*(u**2+v**2)*RT**2*cos(y)
C##############################################################################
C     Enter FORTRAN expressions for the functions whose integrals are to be   #
C     calculated and printed.  They may be functions of                       #
C                                                                             #
C        X,Y,U,Ux,Uy,Uxx,Uyy,Uxy                                              #
C            V,Vx,Vy,Vxx,Vyy,Vxy                                              #
C            ETA,ETAx,ETAy,ETAxx,ETAyy,ETAxy                                  #
C            Q,Qx,Qy,Qxx,Qyy,Qxy and (if applicable) T                        #
C              .   .    .    .    .     .                                     #
C                                                                             #
C     The components (NORMx,NORMy) of the unit outward normal vector          #
C     may also be referenced.                                                 #
C                                                                             #
C     The parameters P1,P2 and derivatives with respect to these may also     #
C     be referenced:                                                          #
C           U1,U2,U11,U22,U12                                                 #
C           V1,V2,V11,V22,V12                                                 #
C           ETA1,ETA2,ETA11,ETA22,ETA12                                       #
C           Q1,Q2,Q11,Q22,Q12                                                 #
C             .    .    .    .     .                                          #
C     You can also reference the normal derivatives Unorm,Vnorm,ETAnorm,      #
C     Qnorm...                                                                #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you only want to integrate a function over part of the boundary,  +#
C     + define that function to be zero on the rest of the boundary.         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C                                                  BND. INTEGRAL1 DEFINED      
C     if (kint8z.eq.-1) yd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
      if (kint8z.gt.0) yd8z = yd8z*dvol                                        
      if (kint8z.lt.0) yd8z = yd8z*darea                                       
                          else                                                 
C##############################################################################
C     Now enter FORTRAN expressions to define the PDE coefficients, which     #
C     may be functions of                                                     #
C                                                                             #
C        X,Y,T,U,Ux,Uy,Uxx,Uyy,Uxy                                            #
C              V,Vx,Vy,Vxx,Vyy,Vxy                                            #
C              ETA,ETAx,ETAy,ETAxx,ETAyy,ETAxy                                #
C              Q,Qx,Qy,Qxx,Qyy,Qxy                                            #
C                .   .    .    .    .     .                                   #
C                                                                             #
C     Recall that the PDEs have the form                                      #
C                                                                             #
C     C11*d(U)/dT + C12*d(V)/dT + C13*d(ETA)/dT + C14*d(Q)/dT +...= F1        #
C     C21*d(U)/dT + C22*d(V)/dT + C23*d(ETA)/dT + C24*d(Q)/dT +...= F2        #
C     C31*d(U)/dT + C32*d(V)/dT + C33*d(ETA)/dT + C34*d(Q)/dT +...= F3        #
C     C41*d(U)/dT + C42*d(V)/dT + C43*d(ETA)/dT + C44*d(Q)/dT +...= F4        #
C                   .               .               .                     .   #
C                                                                             #
C     The parameters P1,P2 and derivatives with respect to these may also     #
C     be referenced (U1 = dU/dP1, etc):                                       #
C           U1,U2,U11,U22,U12                                                 #
C           V1,V2,V11,V22,V12                                                 #
C           ETA1,ETA2,ETA11,ETA22,ETA12                                       #
C           Q1,Q2,Q11,Q22,Q12                                                 #
C             .    .    .    .     .     .     .     .     .                  #
C##############################################################################
            if (mode.eq.1) then
c
c      POT  =  0.75*Omega**2*RT**2*e*( (1-3*sin(y)**2)*cos(omega*t)
c     & + cos(y)**2*(3*cos(2*x)*cos(omega*t)+4*sin(2*x)*sin(omega*t)))
      POTx =  0.75*Omega**2*RT**2*e*cos(y)**2*
     & (-6*sin(2*x)*cos(omega*t) + 8*cos(2*x)*sin(omega*t))
      POTy = -0.75*Omega**2*RT**2*e*( 6*cos(y)*sin(y)*cos(omega*t) +
     &2*sin(y)*cos(y)*(3*cos(2*x)*cos(omega*t)+4*sin(2*x)*sin(omega*t)))  
            else
c
c      POT  = -1.5*Omega**2*RT**2*theta0*sin(y)*cos(y)*
c     & (cos(x-Omega*t) + cos(x+Omega*t)) 
      POTx =  1.5*Omega**2*RT**2*theta0*sin(y)*cos(y)*
     & (sin(x-Omega*t) + sin(x+Omega*t))
      POTy = -1.5*Omega**2*RT**2*theta0*(cos(y)**2-sin(y)**2)*
     & (cos(x-Omega*t) + cos(x+Omega*t))
            endif
                if (j8z.eq.0) then                                             
      yd8z = 0.0                                                               
C                                                  C(1,1) DEFINED              
      if (i8z.eq. -101) yd8z =                                                 
     & 1                                                                       
C                                                  C(1,2) DEFINED              
      if (i8z.eq. -102) yd8z =                                                 
     & 0                                                                       
C                                                  C(1,3) DEFINED              
      if (i8z.eq. -103) yd8z =                                                 
     & 0                                                                       
C                                                  C(1,4) DEFINED              
      if (i8z.eq. -104) yd8z =                                                 
     & 0                                                                       
C                                                  F1 DEFINED
      if (i8z.eq.    1) yd8z =                                                 
     &  2*Omega*V*sin(y) - Cd/h0*sqrt(u**2+v**2)*u              
     & -alpha*U - g/(RT*cos(y))*ETAx + POTx/(RT*cos(y))  
     & -1.0/sigma/(RT*cos(y))*Qx
C                                                  C(2,1) DEFINED              
      if (i8z.eq. -201) yd8z =                                                 
     & 0                                                                       
C                                                  C(2,2) DEFINED              
      if (i8z.eq. -202) yd8z =                                                 
     & 1                                                                       
C                                                  C(2,3) DEFINED              
      if (i8z.eq. -203) yd8z =                                                 
     & 0                                                                       
C                                                  C(2,4) DEFINED              
      if (i8z.eq. -204) yd8z =                                                 
     & 0                                                                       
C                                                  F2 DEFINED                  
      if (i8z.eq.    2) yd8z =                                                 
     & -2*Omega*U*sin(y) - Cd/h0*sqrt(u**2+v**2)*v                  
     & -alpha*V - g/RT*ETAy + POTy/RT
     & - 1.0/sigma/RT*Qy
C                                                  C(3,1) DEFINED              
      if (i8z.eq. -301) yd8z =                                                 
     & 0                                                                       
C                                                  C(3,2) DEFINED              
      if (i8z.eq. -302) yd8z =                                                 
     & 0                                                                       
C                                                  C(3,3) DEFINED              
      if (i8z.eq. -303) yd8z =                                                 
     & 1                                                                       
C                                                  C(3,4) DEFINED              
      if (i8z.eq. -304) yd8z =                                                 
     & 0                                                                       
C                                                  F3 DEFINED                  
      if (i8z.eq.    3) yd8z =                                                 
     & -h0/(RT*cos(y))*(Ux + cos(y)*Vy - sin(y)*V)                             
C                                                  C(4,1) DEFINED              
      if (i8z.eq. -401) yd8z =                                                 
     & 0                                                                       
C                                                  C(4,2) DEFINED              
      if (i8z.eq. -402) yd8z =                                                 
     & 0                                                                       
C                                                  C(4,3) DEFINED              
      if (i8z.eq. -403) yd8z =                                                 
     & 0                                                                       
C                                                  C(4,4) DEFINED              
      if (i8z.eq. -404) yd8z =                                                 
     & 0                                                                       
C                                                  F4 DEFINED                  
      DelETA = ETAxx/cos(y)**2 + ETAyy - ETAy*sin(y)/cos(y)
      DelQ   = Qxx/cos(y)**2 + Qyy - Qy*sin(y)/cos(y)
      if (i8z.eq.    4) yd8z =                                                 
     & betinv*(DelETA + 2*ETA) - (DelQ+2*Q) + (1+Rnu)*Q      
                else                                                           
                endif                                                          
                          endif                                                
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      function u8z(i8z,p1,p2,p38z,t0)                                          
      implicit double precision (a-h,o-z)                                      
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                       
      call dtdpcd(p1,p2,p38z)                                                  
      call dtdpcb(p1,p2,p38z,z18z,z28z,z38z,x,y,z8z,d18z,d28z,d38z,1)          
      u8z = 0.0                                                                
C##############################################################################
C     Now the initial values must be defined using FORTRAN expressions.       #
C     They may be functions of X and Y (and the parameters P1,P2), and may    #
C     also reference the initial time T0.                                     #
C##############################################################################
C                                                  U0 DEFINED                  
      if (i8z.eq.    1) u8z =                                                  
     & 0                                                                       
C                                                  V0 DEFINED                  
      if (i8z.eq.    2) u8z =                                                  
     & 0                                                                       
C                                                  ETA0 DEFINED                
      if (i8z.eq.    3) u8z =                                                  
     & 0                                                                       
C                                                  Q0 DEFINED                  
      if (i8z.eq.    4) u8z =                                                  
     & 0                                                                       
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine gb8z(gd8z,ifac8z,i8z,j8z,p1,p2,p38z,t,uu8z)                   
      implicit double precision (a-h,o-z)                                      
      parameter (neqnmx=  99)                                                  
      dimension uu8z(10,neqnmx)                                                
C        un8z(1,I),un8z(2,I),... hold the (rarely used) values                 
C        of UI,UI1,... from the previous iteration or time step                
      common /dtdp5x/ un8z(10,neqnmx)                                          
      common /dtdp18/norm1,norm2,n38z                                          
      double precision none,norm1,norm2,n38z,normx,normy,nz8z                  
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                   
      none = dtdplx(2)                                                         
      zr8z = 0.0                                                               
      U  = uu8z(1, 1)                                                          
      U1 = uu8z(2, 1)                                                          
      U2 = uu8z(3, 1)                                                          
      V  = uu8z(1, 2)                                                          
      V1 = uu8z(2, 2)                                                          
      V2 = uu8z(3, 2)                                                          
      ETA  = uu8z(1, 3)                                                        
      ETA1 = uu8z(2, 3)                                                        
      ETA2 = uu8z(3, 3)                                                        
      Q  = uu8z(1, 4)                                                          
      Q1 = uu8z(2, 4)                                                          
      Q2 = uu8z(3, 4)                                                          
      call dtdpcd(p1,p2,p38z)                                                  
      call dtdpcb(p1,p2,p38z,norm1,norm2,n38z,x,y,z8z,normx,normy,nz8z,3       
     &)                                                                        
      call dtdpcb(                                                             
     & p1,p2,p38z,U1,U2,zr8z,x,y,z8z,Ux,Uy,uz8z,2)                             
      Unorm = Ux*normx + Uy*normy                                              
      call dtdpcb(                                                             
     & p1,p2,p38z,V1,V2,zr8z,x,y,z8z,Vx,Vy,uz8z,2)                             
      Vnorm = Vx*normx + Vy*normy                                              
      call dtdpcb(                                                             
     & p1,p2,p38z,ETA1,ETA2,zr8z,x,y,z8z,ETAx,ETAy,uz8z,2)                     
      ETAnorm = ETAx*normx + ETAy*normy                                        
      call dtdpcb(                                                             
     & p1,p2,p38z,Q1,Q2,zr8z,x,y,z8z,Qx,Qy,uz8z,2)                             
      Qnorm = Qx*normx + Qy*normy                                              
      if (j8z.eq.0) gd8z = 0.0                                                 
C##############################################################################
C     Enter FORTRAN expressions to define the boundary condition functions,   #
C     which may be functions of                                               #
C                                                                             #
C              X,Y,U,Ux,Uy,                                                   #
C                  V,Vx,Vy,                                                   #
C                  ETA,ETAx,ETAy,                                             #
C                  Q,Qx,Qy and (if applicable) T                              #
C                    .   .    .                                               #
C                                                                             #
C     Recall that the boundary conditions have the form                       #
C                                                                             #
C                           G1 = 0                                            #
C                           G2 = 0                                            #
C                           G3 = 0                                            #
C                           G4 = 0                                            #
C                            .   .                                            #
C     Enter NONE to indicate "no" boundary condition.                         #
C                                                                             #
C     The parameters P1,P2 and derivatives with respect to these may also     #
C     be referenced (U1 = dU/dP1, etc):                                       #
C                              U1,U2                                          #
C                              V1,V2                                          #
C                              ETA1,ETA2                                      #
C                              Q1,Q2                                          #
C                                .    .                                       #
C     The components (NORMx,NORMy) of the unit outward normal vector          #
C     may also be referenced, as well as the normal derivatives Unorm,        #
C     Vnorm,ETAnorm,Qnorm...                                                  #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If "no" boundary condition is specified, the corresponding PDE is    +#
C     + enforced at points just inside the boundary (exactly on the          +#
C     + boundary, if EPS8Z is set to 0 in the main program).                 +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
            if (ifac8z.eq. 1) then                                             
C##############################################################################
C                                                                             #
C     First define the boundary conditions on the face P1 = P1GRID(1).        #
C##############################################################################
                          if (j8z.eq.0) then                                   
C        PERIODIC BOUNDARY CONDITIONS SET (SEE IPERDC)                         
C                                                  G1 DEFINED                  
C     if (i8z.eq.    1) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G2 DEFINED                  
C     if (i8z.eq.    2) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G3 DEFINED                  
C     if (i8z.eq.    3) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G4 DEFINED                  
C     if (i8z.eq.    4) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
                          else                                                 
                          endif                                                
            endif                                                              
            if (ifac8z.eq. 2) then                                             
C##############################################################################
C                                                                             #
C     Now define the boundary conditions on the face P1 = P1GRID(NP1GRID).    #
C##############################################################################
                          if (j8z.eq.0) then                                   
C        PERIODIC BOUNDARY CONDITIONS SET (SEE IPERDC)                         
C                                                  G1 DEFINED                  
C     if (i8z.eq.    1) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G2 DEFINED                  
C     if (i8z.eq.    2) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G3 DEFINED                  
C     if (i8z.eq.    3) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
C                                                  G4 DEFINED                  
C     if (i8z.eq.    4) gd8z =                                                 
C    & [DEFAULT SELECTED, DEFINITION COMMENTED OUT]                            
                          else                                                 
                          endif                                                
            endif                                                              
            if (ifac8z.eq. 3) then                                             
C##############################################################################
C                                                                             #
C     Now define the boundary conditions on the face P2 = P2GRID(1).          #
C##############################################################################
                          if (j8z.eq.0) then                                   
C                                                  G1 DEFINED                  
      if (i8z.eq.    1) gd8z =                                                 
     & none                                                                    
C                                                  G2 DEFINED                  
      if (i8z.eq.    2) gd8z =                                                 
     & none                                                                    
C                                                  G3 DEFINED                  
      if (i8z.eq.    3) gd8z =                                                 
     & none                                                                    
C                                                  G4 DEFINED                  
      if (i8z.eq.    4) gd8z =                                                 
     & none                                                                    
                          else                                                 
                          endif                                                
            endif                                                              
            if (ifac8z.eq. 4) then                                             
C##############################################################################
C                                                                             #
C     Now define the boundary conditions on the face P2 = P2GRID(NP2GRID).    #
C##############################################################################
                          if (j8z.eq.0) then                                   
C                                                  G1 DEFINED                  
      if (i8z.eq.    1) gd8z =                                                 
     & none                                                                    
C                                                  G2 DEFINED                  
      if (i8z.eq.    2) gd8z =                                                 
     & none                                                                    
C                                                  G3 DEFINED                  
      if (i8z.eq.    3) gd8z =                                                 
     & none                                                                    
C                                                  G4 DEFINED                  
      if (i8z.eq.    4) gd8z =                                                 
     & none                                                                    
                          else                                                 
                          endif                                                
            endif                                                              
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine pmod8z(p1,p2,p38z,t,uu8z,uprint,uxprint,uyprint,uzp8z)        
      implicit double precision (a-h,o-z)                                      
      dimension uu8z(10,*),uprint(*),uxprint(*),uyprint(*),uzp8z(*)            
      common/dtdp14/sint(20),bint(20),slim8z(20),blim8z(20)                    
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                          
      zr8z = 0.0                                                               
      U  = uu8z(1, 1)                                                          
      U1 = uu8z(2, 1)                                                          
      U2 = uu8z(3, 1)                                                          
      U11= uu8z(5, 1)                                                          
      U22= uu8z(6, 1)                                                          
      U12= uu8z(8, 1)                                                          
      U21= uu8z(8, 1)                                                          
      V  = uu8z(1, 2)                                                          
      V1 = uu8z(2, 2)                                                          
      V2 = uu8z(3, 2)                                                          
      V11= uu8z(5, 2)                                                          
      V22= uu8z(6, 2)                                                          
      V12= uu8z(8, 2)                                                          
      V21= uu8z(8, 2)                                                          
      ETA  = uu8z(1, 3)                                                        
      ETA1 = uu8z(2, 3)                                                        
      ETA2 = uu8z(3, 3)                                                        
      ETA11= uu8z(5, 3)                                                        
      ETA22= uu8z(6, 3)                                                        
      ETA12= uu8z(8, 3)                                                        
      ETA21= uu8z(8, 3)                                                        
      Q  = uu8z(1, 4)                                                          
      Q1 = uu8z(2, 4)                                                          
      Q2 = uu8z(3, 4)                                                          
      Q11= uu8z(5, 4)                                                          
      Q22= uu8z(6, 4)                                                          
      Q12= uu8z(8, 4)                                                          
      Q21= uu8z(8, 4)                                                          
      call dtdpcd(p1,p2,p38z)                                                  
      call dtdpcc(p1,p2,p38z,                                                  
     &          U1,U2,zr8z,U11,U22,zr8z,U12,zr8z,zr8z,                         
     & x,y,z8z,Ux,Uy,uz8z,Uxx,Uyy,uzz8z,Uxy,uxz8z,uyz8z,                       
     & Uyx,uzx8z,uzy8z,dvol8z,dare8z)                                          
      uxprint( 1) = Ux                                                         
      uyprint( 1) = Uy                                                         
      call dtdpcc(p1,p2,p38z,                                                  
     &          V1,V2,zr8z,V11,V22,zr8z,V12,zr8z,zr8z,                         
     & x,y,z8z,Vx,Vy,uz8z,Vxx,Vyy,uzz8z,Vxy,uxz8z,uyz8z,                       
     & Vyx,uzx8z,uzy8z,dvol8z,dare8z)                                          
      uxprint( 2) = Vx                                                         
      uyprint( 2) = Vy                                                         
      call dtdpcc(p1,p2,p38z,                                                  
     &          ETA1,ETA2,zr8z,ETA11,ETA22,zr8z,ETA12,zr8z,zr8z,               
     & x,y,z8z,ETAx,ETAy,uz8z,ETAxx,ETAyy,uzz8z,ETAxy,uxz8z,uyz8z,             
     & ETAyx,uzx8z,uzy8z,dvol8z,dare8z)                                        
      uxprint( 3) = ETAx                                                       
      uyprint( 3) = ETAy                                                       
      call dtdpcc(p1,p2,p38z,                                                  
     &          Q1,Q2,zr8z,Q11,Q22,zr8z,Q12,zr8z,zr8z,                         
     & x,y,z8z,Qx,Qy,uz8z,Qxx,Qyy,uzz8z,Qxy,uxz8z,uyz8z,                       
     & Qyx,uzx8z,uzy8z,dvol8z,dare8z)                                          
      uxprint( 4) = Qx                                                         
      uyprint( 4) = Qy                                                         
C##############################################################################
C     If you don't want to read the FINE PRINT, default all of the following  #
C     variables.                                                              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Normally, PDE2D saves the values of U,Ux,Uy,V,Vx,Vy,                 +#
C     + ETA,ETAx,ETAy,Q,Qx,Qy...at the output points.  If different          +#
C     + variables are to be saved (for later printing or plotting) the       +#
C     + following functions can be used to re-define the output variables:   +#
C     +    define UPRINT(1) to replace  U                                    +#
C     +           UXPRINT(1)            Ux                                   +#
C     +           UYPRINT(1)            Uy                                   +#
C     +           UPRINT(2)             V                                    +#
C     +           UXPRINT(2)            Vx                                   +#
C     +           UYPRINT(2)            Vy                                   +#
C     +           UPRINT(3)             ETA                                  +#
C     +           UXPRINT(3)            ETAx                                 +#
C     +           UYPRINT(3)            ETAy                                 +#
C     +           UPRINT(4)             Q                                    +#
C     +           UXPRINT(4)            Qx                                   +#
C     +           UYPRINT(4)            Qy                                   +#
C     +                   .              .                                   +#
C     +                   .              .                                   +#
C     + Each function may be a function of                                   +#
C     +                                                                      +#
C     +    X,Y,U,Ux,Uy,Uxx,Uyy,Uxy                                           +#
C     +        V,Vx,Vy,Vxx,Vyy,Vxy                                           +#
C     +        ETA,ETAx,ETAy,ETAxx,ETAyy,ETAxy                               +#
C     +        Q,Qx,Qy,Qxx,Qyy,Qxy and (if applicable) T                     +#
C     +            .   .    .    .    .     .                                +#
C     +                                                                      +#
C     + Each may also be a function of the integral estimates SINT(1),...,   +#
C     + BINT(1),...                                                          +#
C     +                                                                      +#
C     + The parameters P1,P2 and derivatives with respect to these may also  +#
C     + be referenced (U1 = dU/dP1, etc):                                    +#
C     +       U1,U2,U11,U22,U12                                              +#
C     +       V1,V2,V11,V22,V12                                              +#
C     +       ETA1,ETA2,ETA11,ETA22,ETA12                                    +#
C     +       Q1,Q2,Q11,Q22,Q12                                              +#
C     +         .    .    .    .     .                                       +#
C     +                                                                      +#
C     + The default for each variable is no change, for example, UPRINT(1)   +#
C     + defaults to U.  Enter FORTRAN expressions for each of the            +#
C     + following functions (or default).                                    +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C        DEFINE UPRINT(*),UXPRINT(*),UYPRINT(*) HERE:                          
      UXPRINT(1) = SINT(1)
      UYPRINT(1) =
     & Cd*sigma*(u**2+v**2)**1.5
     & + alpha*sigma*h0*(u**2+v**2)
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      function axis8z(i8z,p1,p2,p38z,ical8z)                                   
      implicit double precision (a-h,o-z)                                      
      call dtdpcd(p1,p2,p38z)                                                  
      call dtdpcb(p1,p2,p38z,z18z,z28z,z38z,x,y,z8z,d18z,d28z,d38z,1)          
      if (i8z.eq.1) axis8z = x                                                 
      if (i8z.eq.2) axis8z = y                                                 
      return                                                                   
      end                                                                      
C        dummy routines                                                        
      subroutine xy8z(i8z,iarc8z,s,x,y,s0,sf)                                  
      implicit double precision (a-h,o-z)                                      
      return                                                                   
      end                                                                      
      subroutine dis8z(x,y,ktri,triden,shape)                                  
      implicit double precision (a-h,o-z)                                      
      return                                                                   
      end                                                                      
      function fb8z(i8z,iarc8z,ktri,s,x,y,t)                                   
      implicit double precision (a-h,o-z)                                      
      fb8z = 0                                                                 
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine postpr(tout,nsave,p1out,p2out,np1,np2,uout,neqn,
     & iprob,nprob)
      implicit double precision (a-h,o-z)                                      
      dimension p1out(0:np1,0:np2),p2out(0:np1,0:np2),tout(0:nsave)            
      dimension uout(0:np1,0:np2,4,neqn,0:nsave)                               
      common/parm8z/ pi,Omega ,Cd    ,sigma ,RT    ,h0    ,e                   
     &,theta0,g     ,mode  ,alpha ,hi    ,betinv,Rnu                     
      common /dtdp27/ itask,npes,icomm                                         
      common /dtdp46/ eps8z,cgtl8z,npmx8z,itype,near8z                         
      data lun,lud/0,47/                                                       
      if (itask.gt.0) return                                                   
C     UOUT(I,J,IDER,IEQ,L) = U_IEQ,  if IDER=1                                 
C                            Ux_IEQ, if IDER=2                                 
C                            Uy_IEQ, if IDER=3                                 
C       (possibly as modified by UPRINT,..)                                    
C       at the point (P1OUT(I,J) , P2OUT(I,J))                                 
C       at time/iteration TOUT(L).                                             
C       ******* ADD POSTPROCESSING CODE HERE:                                  
C       IN THE EXAMPLE BELOW, MATLAB PLOTFILES pde2d.m,                        
C      pde2d.rdm CREATED (REMOVE  COMMENTS TO ACTIVATE)                      
c   compute energy dissipation
      L = 0
      do i=1,10
         sum = 0
         do j=1,20
            L = L+1
            sum = sum + uout(1,1,2,1,L) 
         enddo
         sum = sum/20.
         print *,' average dissipation rate this period (W) = ',sum
      enddo
      do i=0,np1
      do j=0,np2
         uout(i,j,3,2,nsave) = 0
      enddo
      enddo
      do L=181,200
         do i=0,np1
         do j=0,np2
            uout(i,j,3,2,nsave) = uout(i,j,3,2,nsave) + uout(i,j,3,1,L)
         enddo
         enddo
      enddo
      do i=0,np1
      do j=0,np2
         uout(i,j,3,2,nsave) = uout(i,j,3,2,nsave)/20.d0
      enddo
      enddo
      if (iprob.lt.nprob) return
C!      if (lun.eq.0) then                                                     
C!         lun = 46                                                            
C!         open (lun,file='pde2d.m')                                           
C!         open (lud,file='pde2d.rdm')                                         
C!         write (lun,*) 'fid = fopen(''pde2d.rdm'');'                         
C!      endif                                                                  
C!      do 78753 l=0,nsave                                                     
C!         if (tout(l).ne.dtdplx(2)) nsave0 = l                                
C!78753 continue                                                               
C!      write (lud,78754) nsave0                                               
C!      write (lud,78754) neqn                                                 
C!      write (lud,78754) np1                                                  
C!      write (lud,78754) np2                                                  
C!78754 format (i8)                                                            
C!      do 78756 i=0,np1                                                       
C!      do 78755 j=0,np2                                                       
C!         p1 = p1out(i,j)                                                     
C!         p2 = p2out(i,j)                                                     
C!         p38z = 0.0                                                          
C!         call dtdpcd(p1,p2,p38z)                                             
C!         call dtdpcb(p1,p2,p38z,z18z,z28z,z38z,x,y,z8z,                      
C!     &   d18z,d28z,d38z,1)                                                   
C!         write (lud,78762) p1,p2,x,y                                         
C!78755 continue                                                               
C!78756 continue                                                               
C!      do 78761 l=0,nsave0                                                    
C!         write (lud,78762) tout(l)                                           
C!         do 78760 ieq=1,neqn                                                 
C!         do 78759 ider=1,3                                                   
C!         do 78758 i=0,np1                                                    
C!         do 78757 j=0,np2                                                    
C!            write (lud,78762) uout(i,j,ider,ieq,l)                           
C!78757    continue                                                            
C!78758    continue                                                            
C!78759    continue                                                            
C!78760    continue                                                            
C!78761 continue                                                               
C!78762 format (e16.8)                                                         
C       ******* WRITE pde2d.m                                                  
C!      call mtdp2dc(itype,lun)                                                
      return                                                                   
      end                                                                      
