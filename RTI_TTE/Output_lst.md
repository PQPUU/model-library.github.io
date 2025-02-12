# Time-to-first respiratory tract infection

## Model output

<a href="./Output_TTE_RTI_BCG_only.lst">Download here</a>

```
Tue Apr 26 11:31:48 CEST 2022
$PROBLEM    Time-to-event analysis of respiratory tract infection
;           data,ITT dataset
; the same model code for PP dataset

; MTIME is used for simulations https://www.page-meeting.org/?abstract=3166

; To create simulation dataset for VPC, use update_inits command from PsN with -flip_comments option
$INPUT      ID DV TIME EVID TYPE FTIME GROUP
; DV - status in the original data, categorical, DV=0 is "no event", DV=1 "event"

; TIME - time in the original data, unit=days

; EVID - Event ID, 3=the first record of the individual, 0=other records

; TYPE - 1 - time zero, 2=event/censoring time, 3=maximum time for simulations

; FTIME - final time for VPCs

; GROUP - vaccination group, 1=placebo, 2=BCG
$DATA      RTI_TTE_ITT_20211203.csv IGNORE=@ IGNORE=(TYPE.EQ.3)

;Sim_start

;IGNORE=(TYPE.EQ.2)
$SUBROUTINE ADVAN=13 TOL=9
;$ABBREVIATED COMRES=7

;$SUBROUTINE ADVAN=13 TOL=9

;Sim_end
$MODEL      COMP=(HAZARD)
$PK
IF(GROUP.EQ.1) RFGROUP = 0
IF(GROUP.EQ.2) RFGROUP = THETA(3)
 LAM = THETA(1)*EXP(ETA(1))/100
 SHP = THETA(2)
;Sim_start
;;;;;;;;; THE SIMULATION PART FOR TTE SIMULATIONS ;;
; IF (ICALL.EQ.4) THEN                         ; The event time sim $problem
;     IF (NEWIND.EQ.0) THEN                    ; Only for the first record
;         COM(6)  =  1                         ; Reset simulation ID counter
;        COM(4)  =  459                        ; Set max time/censoring time
;     ENDIF
;     IF (NEWIND.EQ.1) THEN                    ; For every new ind except first in dataset
;         ICOUNT  =  COM(6) + 1                  ; Update individual counter over simulations
;         COM(6)  =  ICOUNT
;     ENDIF
;     IF (NEWIND.NE.2) THEN                    ; For every new individual
;         CALL RANDOM(2,R)
;         COM(3)  =  -1                          ; Variable for survival at event time
;         COM(2)  =  R                           ; Store the random number
;         COM(1)  =  -1                          ; Variable for the event time
;         COM(7)  =  0                           ; Individual event counter
;     ENDIF
; ENDIF
;;;;;;---------MTIME for increasing precision in $DES --------
; IF (NEWIND.NE.2) THEN
;     TEMP = 0
; ENDIF
; TEMP = TEMP+0.1
; MTIME(1) = TEMP
; MTDIFF = 1
;Sim_end
$DES
 DADT(1) = LAM*EXP(SHP*T)*EXP(RFGROUP)
;Sim_start
; SUR = EXP(-A(1))
; IF(COM(2).GT.SUR.AND.COM(1).EQ.-1) THEN ; If event save event time in COM(1)
;    COM(1) = T
;    COM(3) = SUR
; ENDIF
;Sim_end
$ERROR
;Sim_start
 IF(NEWIND.NE.2) OLDCHZ = 0           ;reset the cumulative hazard
 CHZ = A(1)-OLDCHZ                  ;cumulative hazard from previous time point in data set
 OLDCHZ = A(1)                      ;rename old cumulative hazard
 SUR = EXP(-CHZ)                    ;survival probability
 HAZNOW=LAM*EXP(SHP*TIME)*EXP(RFGROUP)              ; rate of event each time pt NB: update with each new model
 IF(DV.EQ.0) Y = SUR                ;censored event (prob of survival)
 IF(DV.NE.0) Y = SUR*HAZNOW         ;prob density function of event
;"FIRST
;"@CHARACTER(LEN=100)::FMT ! Define FORMAT string for writing dataset
;;;;;;; NORMAL TTE MODEL
;CHZ = A(1)
;SURX = EXP(-CHZ)              ;survival probability
;IF (COM(1).GT.COM(4)) THEN      ;IF T > ENDTIME, T=ENDTIME
;;;;;;; Check survival again at endtime
;IF (COM(2).GT.SURX) THEN
;COM(1) = COM(4)
;ELSE
;COM(1) = -1              ;Integrated too far, reset event
;ENDIF
;ENDIF
;EVT = COM(1)        ; Save Event time
;RNM = COM(2)        ; Save random number, just for debugging
;ENDTIME = COM(4)          ; Endtime of study
;TT = COM(5)            ; Analytic event time
;;;;;;;;; ADD RTTE, DV TO OUTPUT, SET DV=0 IF NO EVENT OR CENSORED, DV=1 IF EVENT, RTTE = 1 IF EVENT OR CENSORED
;IF (ICALL.EQ.4) THEN          ; Initate DV to 0 (No event)
;DV=0
;ENDIF
;TMDV = 0
;IF (EVID.GE.2) THEN         ;Set MDV variable for output
;TMDV=1
;ENDIF
;ICOUNT = COM(6)+(IREP-1)*NINDR
;ITER = IREP
;;;;;;;;;;; Define the format of the output file
;"LAST
;"FMT='(E13.7,9(1XE13.7))' ! The output FORMAT
;" ! Write all events
;" IF (NEWIND.EQ.0) THEN !Open file at first record
;" OPEN (99, FILE = 'simtab.dat', POSITION='APPEND')
;" IF (IREP.EQ.1) THEN !Write header for 1st subproblem
;" WRITE (99,'(A,9(1XA))') 'ID','DV','TIME','RTTE','SURX','ICOUNT','ITER','RAND','ENDTIME','GROUP'
;" ENDIF
;" ENDIF
;" IF (EVT.NE.-1) THEN !If an EVENT
;" DV=1
;" RTTE=1
;" TMDV=0
;" ! Write SIM specific output
;" WRITE (99,FMT) ID,DV,EVT,RTTE,COM(3),ICOUNT,ITER,COM(2),ENDTIME,GROUP
;" COM(1) = -1 !Reset Event time variable
;" COM(2) = 0 !Reset Random variable
;" COM(3) = -1 !Reset survival variable
;" COM(7) = COM(7) + 1 !Update Event counter
;" ELSE IF (LIREC.EQ.NDREC.AND.COM(7).EQ.0) THEN !Right Censoring (if no previous events)
;" DV=0
;" TMDV=0
;" RTTE=1
;" TMP=COM(4)
;" WRITE (99,FMT) ID,DV,TMP,RTTE,SURX,ICOUNT,ITER,COM(2),ENDTIME,GROUP
;" ENDIF
;" IF (NDREC.EQ.LIREC.AND.NIREC.EQ.NINDR) THEN ! Last record for last individual
;" CLOSE(99) ! Close File pointer
;" ENDIF
;Sim_end
$THETA  (0,0.400657) ; LAMBDA/100
 -0.00349062 ; ALPHA
 0.071172 ; RFGROUP
$OMEGA  0  FIX
;Sim_start
$ESTIMATION MAXEVAL=9990 METHOD=0 LIKE PRINT=1 SIGL=9 NSIG=3
$COVARIANCE PRINT=E
$TABLE      ID TIME SUR HAZNOW LAM SHP RFGROUP EVID NOPRINT ONEHEADER
            FILE=mytab
;$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION NSUB = 1

;Sim_end


NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  3) THERE MAY BE AN ERROR IN THE ABBREVIATED CODE. THE FOLLOWING
 ONE OR MORE RANDOM VARIABLES ARE DEFINED WITH "IF" STATEMENTS THAT DO NOT
 PROVIDE DEFINITIONS FOR BOTH THE "THEN" AND "ELSE" CASES. IF ALL
 CONDITIONS FAIL, THE VALUES OF THESE VARIABLES WILL BE ZERO.
  
   Y


License Registered to: Uppsala University
Expiration Date:    14 AUG 2022
Current Date:       26 APR 2022
Days until program expires : 108
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 Time-to-event analysis of respiratory tract infection
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     2000
 NO. OF DATA ITEMS IN DATA SET:   8
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   2
 MDV DATA ITEM IS DATA ITEM NO.:  8
0INDICES PASSED TO SUBROUTINE PRED:
   4   3   0   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID DV TIME EVID TYPE FTIME GROUP MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 RFGROUP LAM SHP SUR HAZNOW
0FORMAT FOR DATA:
 (7E5.0,1F2.0)

 TOT. NO. OF OBS RECS:     1000
 TOT. NO. OF INDIVIDUALS:     1000
0LENGTH OF THETA:   3
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.4007E+00     0.1000E+07
 -0.1000E+07    -0.3491E-02     0.1000E+07
 -0.1000E+07     0.7117E-01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.0000E+00
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:             YES
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:             3U
 MC SAMPLES (ESAMPLE):    300
 WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
0-- TABLE   1 --
0RECORDS ONLY:    ALL
04 COLUMNS APPENDED:    YES
 PRINTED:                NO
 HEADERS:               ONE
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 IDFORMAT:
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME SUR HAZNOW LAM SHP RFGROUP EVID
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (LSODA, ADVAN13)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         HAZARD       ON         YES        YES        YES        YES
    2         OUTPUT       OFF        YES        NO         NO         NO
 INITIAL (BASE) TOLERANCE SETTINGS:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      4
   TIME DATA ITEM IS DATA ITEM NO.:          3

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1


 #TBLN:      1
 #METH: First Order

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 EPS-ETA INTERACTION:                     NO
 PRED F SET TO A LIKELIHOOD:              YES
 NO. OF FUNCT. EVALS. ALLOWED:            9990
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      9
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     9
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): psn.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR TABLE/SCATTER STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   9
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12

 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=NPRED
 RES=NRES
 WRES=NWRES
 IWRS=NIWRES
 IPRD=NIPRED
 IRS=NIRES

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   7854.48733585122        NO. OF FUNC. EVALS.:   4
 CUMULATIVE NO. OF FUNC. EVALS.:        4
 NPARAMETR:  4.0066E-01 -3.4906E-03  7.1172E-02
 PARAMETER:  1.0000E-01 -1.0000E-01  1.0000E-01
 GRADIENT:  -2.9079E-02 -3.9207E-01 -1.2692E-02

0ITERATION NO.:    1    OBJECTIVE VALUE:   7854.48733284285        NO. OF FUNC. EVALS.:  12
 CUMULATIVE NO. OF FUNC. EVALS.:       16
 NPARAMETR:  4.0066E-01 -3.4902E-03  7.1172E-02
 PARAMETER:  1.0000E-01 -9.9987E-02  1.0000E-01
 GRADIENT:   3.8902E-02  4.8045E-02  1.4421E-02

0ITERATION NO.:    2    OBJECTIVE VALUE:   7854.48733284285        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       29
 NPARAMETR:  4.0066E-01 -3.4902E-03  7.1172E-02
 PARAMETER:  1.0000E-01 -9.9987E-02  1.0000E-01
 GRADIENT:   3.4044E-02 -2.5633E-03  1.1856E-02

0ITERATION NO.:    3    OBJECTIVE VALUE:   7854.48733262230        NO. OF FUNC. EVALS.:  10
 CUMULATIVE NO. OF FUNC. EVALS.:       39
 NPARAMETR:  4.0062E-01 -3.4899E-03  7.1151E-02
 PARAMETER:  9.9916E-02 -9.9980E-02  9.9971E-02
 GRADIENT:  -4.4462E-02 -2.5107E-01 -2.1131E-02

0ITERATION NO.:    4    OBJECTIVE VALUE:   7854.48733219206        NO. OF FUNC. EVALS.:  11
 CUMULATIVE NO. OF FUNC. EVALS.:       50
 NPARAMETR:  4.0061E-01 -3.4897E-03  7.1146E-02
 PARAMETER:  9.9894E-02 -9.9972E-02  9.9963E-02
 GRADIENT:  -3.4992E-02 -1.1387E-01 -1.8869E-02

0ITERATION NO.:    5    OBJECTIVE VALUE:   7854.48733212944        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  4.0061E-01 -3.4896E-03  7.1168E-02
 PARAMETER:  9.9880E-02 -9.9972E-02  9.9994E-02
 GRADIENT:  -3.6059E-02 -1.1312E-01 -1.4898E-02

0ITERATION NO.:    6    OBJECTIVE VALUE:   7854.48733212944        NO. OF FUNC. EVALS.:   0
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  4.0061E-01 -3.4896E-03  7.1168E-02
 PARAMETER:  9.9880E-02 -9.9972E-02  9.9994E-02
 GRADIENT:  -3.6059E-02 -1.1312E-01 -1.4898E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:       59
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):            0
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    0.0000000000000000     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    7854.4873321294372     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       7854.4873321294372     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                             0
  
 #TERE:
 Elapsed estimation  time in seconds:     1.21
 Elapsed covariance  time in seconds:     0.50
 Elapsed postprocess time in seconds:     0.07
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     7854.487       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         4.01E-01 -3.49E-03  7.12E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1     
 
 ETA1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1     
 
 ETA1
+        0.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         3.34E-02  4.60E-04  8.42E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1     
 
 ETA1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1     
 
 ETA1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        1.12E-03
 
 TH 2
+       -1.06E-05  2.12E-07
 
 TH 3
+       -1.51E-03  1.13E-06  7.09E-03
 
 OM11
+       ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        3.34E-02
 
 TH 2
+       -6.91E-01  4.60E-04
 
 TH 3
+       -5.37E-01  2.91E-02  8.42E-02
 
 OM11
+       ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        3.52E+03
 
 TH 2
+        1.73E+05  1.32E+07
 
 TH 3
+        7.24E+02  3.48E+04  2.90E+02
 
 OM11
+       ......... ......... ......... .........
 
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                      EIGENVALUES OF COR MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

             1         2         3
 
         1.38E-01  9.72E-01  1.89E+00
 
 Elapsed finaloutput time in seconds:     0.08
 #CPUT: Total CPU Time in Seconds,        2.107
Stop Time:
Tue Apr 26 11:32:04 CEST 2022
```

[Back](../rti_tte_main)

[Home](../../model-library.github.io/)
