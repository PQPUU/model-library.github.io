# Respiratory tract infection over time Markov model

## Model output

<a href="./Output_Markov_RTI_BCGonly.lst">Download here</a>

```
03/03/2022 
18:06
;; 1. Based on: 268
;; 2. Description:
;;    BCG on P03, P23; IIV P21, P20, P10; Logit transformation + $PRED parameterization
;; 3. Label:
;;    Base Markov Chain for Respiratory Tract Infection Health Status Score 
;; 4. Structural model:
;;    4 state Markov Chain for HS0 (healthy), HS1 (mild), HS2 (moderage), and (lumped) HS3-4-5-6-7 (severe and up)        
;; 5. Covariate model: 
;;    NA
;; 6. Interindividual variability:
;;    P21, P20, P10
;; 7. Interoccasion variability:
;;    NA
;; 8. Residual variabi1lity:
;;    NA
;; 9. Estimation:
;;    Laplace
$PROBLEM    Markov Chain for Respiratory Tract Infection Health Status
$INPUT      ID EVENT TIME DV PDV GROUP
;ID trial participant identifier

;EVENT event number in case of multiple events

;TIME timepoint in weeks

;DV health status (HS) score from 0 (healthy), 1 (mild), 2 (moderate), 3 (severe), 4 (hospitalization), 5 (hospitalization with oxygen), 6 (hospitalization with ventilation) to 7 (death)

;PDV previous DV

;GROUP placebo (0) or BCG (1)

;Sim_start
$DATA      Real_Markov_RTI_BCGonly.csv IGNORE=@

;$DATA      Real_Markov_RTI_BCGonly_SIM.csv IGNORE=@

;Sim_end
$PRED

;Typical value of the probability (between 0 and 1) of transitioning between Markov states

TVP01 = THETA(1)
TVP02 = THETA(2)
TVP03 = THETA(3)
IF(GROUP.EQ.2) TVP03 = THETA(14)
TVP10 = THETA(4)
TVP12 = THETA(5)
TVP13 = THETA(6)
TVP20 = THETA(7)
TVP21 = THETA(8)
TVP23 = THETA(9)
IF(GROUP.EQ.2) TVP23 = THETA(13)
TVP30 = THETA(10)
TVP31 = THETA(11)
TVP32 = THETA(12)

;Logit transformation

PL01 = LOG(TVP01 / (1 - TVP01))
PL02 = LOG(TVP02 / (1 - TVP02))
PL03 = LOG(TVP03 / (1 - TVP03))
PL10 = LOG(TVP10 / (1 - TVP10))
PL12 = LOG(TVP12 / (1 - TVP12))
PL13 = LOG(TVP13 / (1 - TVP13))
PL20 = LOG(TVP20 / (1 - TVP20))
PL21 = LOG(TVP21 / (1 - TVP21))
PL23 = LOG(TVP23 / (1 - TVP23))
PL30 = LOG(TVP30 / (1 - TVP30))
PL31 = LOG(TVP31 / (1 - TVP31))
PL32 = LOG(TVP32 / (1 - TVP32))

;IIV 

LGT01 = PL01 
LGT02 = PL02 
LGT03 = PL03 
LGT10 = PL10 + ETA(1)
LGT12 = PL12  
LGT13 = PL13 
LGT20 = PL20 + ETA(2)
LGT21 = PL21 + ETA(3) 
LGT23 = PL23 
LGT30 = PL30 
LGT31 = PL31 
LGT32 = PL32

;Probability backcalculated by logistic transformation

P01 = EXP(LGT01) / (1 + EXP(LGT01))
P02 = EXP(LGT02) / (1 + EXP(LGT02))
P03 = EXP(LGT03) / (1 + EXP(LGT03))
P10 = EXP(LGT10) / (1 + EXP(LGT10))
P12 = EXP(LGT12) / (1 + EXP(LGT12))
P13 = EXP(LGT13) / (1 + EXP(LGT13))
P20 = EXP(LGT20) / (1 + EXP(LGT20))
P21 = EXP(LGT21) / (1 + EXP(LGT21))
P23 = EXP(LGT23) / (1 + EXP(LGT23))
P30 = EXP(LGT30) / (1 + EXP(LGT30))
P31 = EXP(LGT31) / (1 + EXP(LGT31))
P32 = EXP(LGT32) / (1 + EXP(LGT32))

;transition probabilities

;from state 0
T01 = P01
T02 = P02 * (1 - P01) ;fraction from remaining probability
T03 = P03 * (1 - P01 - (P02 * (1 - P01))) ;fraction from remaining probability
T00= 1 - P01 - P02 * (1 - P01) - P03 * (1 - P01 - (P02 * (1 - P01))) ;fraction from remaining probability 

;from state 1
T10 = P10
T12 = P12 * (1 - P10) ;fraction from remaining probability
T13 = P13 * (1 - P10 - (P12 * (1 - P10))) ;fraction from remaining probability
T11 = 1 - P10 - P12 * (1 - P10) - P13 * (1 - P10 - (P12 * (1 - P10))) ;fraction from remaining probability

;from state 2
T20 = P20
T21 = P21 * (1 - P20) ;fraction from remaining probability
T23 = P23 * (1 - P20 - (P21 * (1 - P20))) ;fraction from remaining probability
T22 = 1 - P20 - P21 * (1 - P20) - P23 * (1 - P20 - (P21 * (1 - P20))) ;fraction from remaining probability

;from state 3
T30 = P30
T31 = P31 * (1 - P30) ;fraction from remaining probability
T32 = P32 * (1 - P30 - (P31 * (1 - P30))) ;fraction from remaining probability
T33 = 1 - P30 - P31 * (1 - P30) - P32 * (1 - P30 - (P31 * (1 - P30))) ;fraction from remaining probability

;cumulative probabilities

CP00 = T00
CP01 = CP00 + T01
CP02 = CP01 + T02
CP03 = 1 - CP02

CP10 = T10
CP11 = CP10 + T11
CP12 = CP11 + T12
CP13 = 1 - CP12

CP20 = T20
CP21 = CP20 + T21
CP22 = CP21 + T22
CP23 = 1 - CP22

CP30 = T30
CP31 = CP30 + T31
CP32 = CP31 + T32
CP33 = 1 - CP32

;Sim_start

IF(PDV.EQ.0.AND.DV.EQ.0) Y = T00
IF(PDV.EQ.0.AND.DV.EQ.1) Y = T01
IF(PDV.EQ.0.AND.DV.EQ.2) Y = T02
IF(PDV.EQ.0.AND.DV.EQ.3) Y = T03

IF(PDV.EQ.1.AND.DV.EQ.0) Y = T10
IF(PDV.EQ.1.AND.DV.EQ.1) Y = T11
IF(PDV.EQ.1.AND.DV.EQ.2) Y = T12 
IF(PDV.EQ.1.AND.DV.EQ.3) Y = T13 

IF(PDV.EQ.2.AND.DV.EQ.0) Y = T20
IF(PDV.EQ.2.AND.DV.EQ.1) Y = T21
IF(PDV.EQ.2.AND.DV.EQ.2) Y = T22
IF(PDV.EQ.2.AND.DV.EQ.3) Y = T23

IF(PDV.EQ.3.AND.DV.EQ.0) Y = T30
IF(PDV.EQ.3.AND.DV.EQ.1) Y = T31 
IF(PDV.EQ.3.AND.DV.EQ.2) Y = T32 
IF(PDV.EQ.3.AND.DV.EQ.3) Y = T33

;;-------- simulation block ---------------------------------------------
;
;   IF (ICALL.EQ.4) THEN
;      CALL RANDOM (2,R)  ; generate a random number between 0 and 1 (uniform distribution)
;
;      IF (NEWIND.NE.2) THEN   ; reset for new individual
;         PREV   = 0  ; everyone starts with HS0
;         DV     = 0  ; observe HS0 at t=0
;		  PDVOUT = 0 ;output PDV
;      ENDIF
;
;     ;if the preceding DV is 0, allocate 0, 1, 2, 3 or higher 
;      IF (PREV.EQ.0) THEN    
;         IF(R.LE.CP00) DV=0
;         IF(R.GT.CP00.AND.R.LE.CP01) DV=1
;         IF(R.GT.CP01.AND.R.LE.CP02) DV=2
;         IF(R.GT.CP02) DV=3
;      ENDIF
;
;     ;if the preceding DV is 1, allocate 0, 1, 2, 3 or higher 
;      IF (PREV.EQ.1) THEN    
;         IF(R.LE.CP10) DV=0
;         IF(R.GT.CP10.AND.R.LE.CP11) DV=1
;         IF(R.GT.CP11.AND.R.LE.CP12) DV=2
;         IF(R.GT.CP12) DV=3
;      ENDIF
;	  
;     ;if the preceding DV is 2, allocate 0, 1, 2, 3 or higher 
;      IF (PREV.EQ.2) THEN    
;         IF(R.LE.CP20) DV=0
;         IF(R.GT.CP20.AND.R.LE.CP21) DV=1
;         IF(R.GT.CP21.AND.R.LE.CP22) DV=2
;         IF(R.GT.CP22) DV=3
;      ENDIF
;
;     ;if the preceding DV is 3, allocate 0, 1, 2, 3 or higher 
;      IF (PREV.EQ.3) THEN    
;         IF(R.LE.CP30) DV=0
;         IF(R.GT.CP30.AND.R.LE.CP31) DV=1
;         IF(R.GT.CP31.AND.R.LE.CP32) DV=2
;         IF(R.GT.CP32) DV=3
;      ENDIF	  
;
;   ENDIF
;
;;-------- previous score ----------------------------------------------
;
;PDVOUT = PREV
;PREV = DV
;
;Sim_end


$THETA  (0,0.015,1) ; P01
 (0,0.007,1) ; P02
 (0,0.0009,1) ; P03
 (0,0.45,1) ; P10
 (0,0.02,1) ; P12
 (0,0.005,1) ; P13
 (0,0.3,1) ; P20
 (0,0.6,1) ; P21
 (0,0.1,1) ; P23
 (0,0.05,1) ; P30
 (0,0.5,1) ; P31
 (0,0.6,1) ; P32
 (0,0.1,1) ; P23_BCG
 (0,0.0009,1) ; P03_BCG
$OMEGA  0.09  ;    IIV_P10
 0.09  ;    IIV_P20
 0.09  ;    IIV_P21
;Sim_start
;
;$SIGMA 0 FIX ; DUMMY for simulation
;
$ESTIMATION METH=COND LAPLACE LIKE MAXEVAL=9999
$COVARIANCE
$TABLE      ID EVENT TIME DV P01 P02 P03 P10 P12 P13 P20 P21 P23 P30
            P31 P32 CP00 CP01 CP02 CP03 CP10 CP11 CP12 CP13 CP20 CP21
            CP22 CP23 CP30 CP31 CP32 CP33 T00 T01 T02 T03 T10 T11 T12
            T13 T20 T21 T22 T23 T30 T31 T32 T33 ETAS(1:LAST) NOPRINT
            ONEHEADER FILE=mytab
;$SIMULATION (12345) (54321 UNIFORM) SUB=50 ONLYSIM

;$TABLE      ID EVENT TIME P01 P02 P03 P10 P12 P13 P20

;            P21 P23 P30 P31 P32 CP00 CP01 CP02 CP03 CP10 CP11 CP12 CP13

;			 CP20 CP21 CP22 CP23 CP30 CP31 CP32 CP33 T00 T01 T02 T03

;            T10 T11 T12 T13 T20 T21 T22 T23 T30 T31 T32 T33 ETAS(1:LAST) 

;            NOPRINT ONEHEADER FILE=mysimtab

;Sim_end

  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  13) WITH USER-WRITTEN PRED OR $PRED, NM-TRAN CANNOT APPEND THE
 MDV DATA ITEM.
  
License Registered to: Uppsala University
Expiration Date:    14 AUG 2022
Current Date:        3 MAR 2022
Days until program expires : 161
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 Markov Chain for Respiratory Tract Infection Health Status
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:    51368
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   4
0LABELS FOR DATA ITEMS:
 ID EVENT TIME DV PDV GROUP
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 P01 P02 P03 P10 P12 P13 P20 P21 P23 P30 P31 P32 T01 T02 T03 T00 T10 T12 T13 T11 T20 T21 T23 T22 T30 T31 T32 T33 CP00 CP01
 CP02 CP03 CP10 CP11 CP12 CP13 CP20 CP21 CP22 CP23 CP30 CP31 CP32 CP33
0FORMAT FOR DATA:
 (6E5.0)

 TOT. NO. OF OBS RECS:    51368
 TOT. NO. OF INDIVIDUALS:     1000
0LENGTH OF THETA:  14
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1500E-01     0.1000E+01
  0.0000E+00     0.7000E-02     0.1000E+01
  0.0000E+00     0.9000E-03     0.1000E+01
  0.0000E+00     0.4500E+00     0.1000E+01
  0.0000E+00     0.2000E-01     0.1000E+01
  0.0000E+00     0.5000E-02     0.1000E+01
  0.0000E+00     0.3000E+00     0.1000E+01
  0.0000E+00     0.6000E+00     0.1000E+01
  0.0000E+00     0.1000E+00     0.1000E+01
  0.0000E+00     0.5000E-01     0.1000E+01
  0.0000E+00     0.5000E+00     0.1000E+01
  0.0000E+00     0.6000E+00     0.1000E+01
  0.0000E+00     0.1000E+00     0.1000E+01
  0.0000E+00     0.9000E-03     0.1000E+01
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
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
 ID EVENT TIME DV P01 P02 P03 P10 P12 P13 P20 P21 P23 P30 P31 P32 CP00 CP01 CP02 CP03 CP10 CP11 CP12 CP13 CP20 CP21 CP22
 CP23 CP30 CP31 CP32 CP33 T00 T01 T02 T03 T10 T11 T12 T13 T20 T21 T22 T23 T30 T31 T32 T33 ETA1 ETA2 ETA3
1


 #TBLN:      1
 #METH: Laplacian Conditional Estimation

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     NO
 LAPLACIAN OBJ. FUNC.:                    YES
 NUMERICAL 2ND DERIVATIVES:               NO
 PRED F SET TO A LIKELIHOOD:              YES
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
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


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=NPRED
 RES=NRES
 WRES=NWRES
 IWRS=NIWRES
 IPRD=NIPRED
 IRS=NIRES

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   14285.4790244904        NO. OF FUNC. EVALS.:  19
 CUMULATIVE NO. OF FUNC. EVALS.:       19
 NPARAMETR:  1.5000E-02  7.0000E-03  9.0000E-04  4.5000E-01  2.0000E-02  5.0000E-03  3.0000E-01  6.0000E-01  1.0000E-01  5.0000E-02
             5.0000E-01  6.0000E-01  1.0000E-01  9.0000E-04  9.0000E-02  9.0000E-02  9.0000E-02
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   4.1994E+02 -9.8022E+01  2.1915E+01  2.0577E+02  8.8171E-01  5.5043E-01  3.5795E+00  8.4688E+01  1.1001E+01  3.0011E-01
             2.0005E+00  2.4003E+00 -7.7991E+00 -4.4234E+00 -1.1083E+02 -1.0231E+01 -3.1777E+01

0ITERATION NO.:   27    OBJECTIVE VALUE:   14021.5557421650        NO. OF FUNC. EVALS.: 652
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0761E-02  8.0007E-03  4.5089E-04  4.7050E-01  1.9495E-02  4.6783E-03  3.0344E-01  5.2846E-01  4.7616E-02  4.3477E-02
             4.5454E-01  4.9999E-01  1.3861E-01  9.9141E-04  1.2633E+00  1.0268E+00  9.0011E-01
 PARAMETER: -2.3638E-01  2.3463E-01 -5.9163E-01  1.8255E-01  7.3927E-02  3.3171E-02  1.1633E-01 -1.9151E-01 -6.9858E-01 -4.6637E-02
            -8.2350E-02 -3.0552E-01  4.7035E-01  1.9682E-01  1.4208E+00  1.3172E+00  1.2514E+00
 GRADIENT:   7.9916E-04  3.9502E-03 -3.9543E-04  3.3212E-03 -1.6418E-04 -1.2539E-04 -7.8617E-04 -2.6496E-04 -6.4421E-04 -6.3710E-05
            -3.1356E-04 -3.4625E-04 -6.1907E-04 -9.6878E-05 -3.8096E-03  1.3815E-03  1.5120E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      671
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.6507E-03  9.7255E-03 -5.2021E-04
 SE:             1.3337E-02  8.2214E-03  6.3934E-03
 N:                    1000        1000        1000

 P VAL.:         7.8430E-01  2.3683E-01  9.3515E-01

 ETASHRINKSD(%)  6.2458E+01  7.4330E+01  7.8680E+01
 ETASHRINKVR(%)  8.5906E+01  9.3411E+01  9.5454E+01
 EBVSHRINKSD(%)  6.0389E+01  7.3410E+01  7.7859E+01
 EBVSHRINKVR(%)  8.4310E+01  9.2930E+01  9.5098E+01
 RELATIVEINF(%)  1.5690E+01  7.0701E+00  4.9021E+00
 EPSSHRINKSD(%)  1.0000E+02
 EPSSHRINKVR(%)  1.0000E+02

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):            0
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    0.0000000000000000     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    14021.555742165043     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       14021.555742165043     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1623
  
 #TERE:
 Elapsed estimation  time in seconds:   176.11
 Elapsed covariance  time in seconds:   318.16
 Elapsed postprocess time in seconds:     0.84
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    14021.556       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12     
          TH13      TH14     
 
         1.08E-02  8.00E-03  4.51E-04  4.71E-01  1.95E-02  4.68E-03  3.03E-01  5.28E-01  4.76E-02  4.35E-02  4.55E-01  5.00E-01
          1.39E-01  9.91E-04
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        1.26E+00
 
 ETA2
+        0.00E+00  1.03E+00
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        1.12E+00
 
 ETA2
+        0.00E+00  1.01E+00
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12     
          TH13      TH14     
 
         5.18E-04  4.47E-04  2.35E-04  2.78E-02  4.60E-03  2.56E-03  2.78E-02  3.69E-02  3.61E-02  4.40E-02  1.17E-01  1.18E-01
          5.30E-02  2.72E-04
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        2.21E-01
 
 ETA2
+       .........  4.24E-01
 
 ETA3
+       ......... .........  4.12E-01
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.84E-02
 
 ETA2
+       .........  2.09E-01
 
 ETA3
+       ......... .........  2.17E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 TH 1
+        2.69E-07
 
 TH 2
+        2.43E-08  2.00E-07
 
 TH 3
+       -2.55E-10 -2.02E-10  5.53E-08
 
 TH 4
+       -1.36E-06 -2.21E-06 -1.66E-07  7.74E-04
 
 TH 5
+       -3.22E-08 -7.35E-08 -7.27E-09 -2.76E-07  2.12E-05
 
 TH 6
+        3.30E-08 -6.59E-08 -7.07E-09 -1.69E-06  2.59E-06  6.57E-06
 
 TH 7
+       -2.01E-06 -1.62E-07  7.80E-08  1.02E-04  1.21E-06  2.58E-06  7.73E-04
 
 TH 8
+        1.77E-06  2.69E-07 -1.18E-07  1.22E-05 -1.08E-05 -7.41E-06  2.61E-05  1.36E-03
 
 TH 9
+       -4.97E-07 -1.44E-07 -2.40E-08  5.07E-06 -1.92E-05  9.30E-09  4.93E-05 -1.11E-04  1.31E-03
 
 TH10
+        4.89E-07  7.28E-07 -3.74E-07  1.31E-05 -8.38E-06 -1.19E-05 -3.55E-05  1.21E-04  2.48E-06  1.93E-03
 
 TH11
+       -1.63E-06  2.01E-06  1.33E-06  1.16E-04  5.10E-06  5.55E-07 -7.79E-05 -1.95E-04  1.74E-07  3.09E-04  1.37E-02
 
 TH12
+        1.79E-06 -7.76E-07 -1.09E-10 -1.69E-04 -1.01E-05 -2.13E-06 -7.88E-05  4.01E-04 -6.41E-07  2.46E-04  6.90E-04  1.39E-02
 
 TH13
+       -1.19E-06 -2.89E-08  7.38E-09  1.08E-04  2.41E-06 -2.60E-05 -3.84E-05 -1.43E-04 -9.60E-07 -4.56E-05  7.86E-04 -3.60E-04
          2.81E-03
 
 TH14
+        1.15E-09 -4.28E-09 -5.59E-13 -9.22E-07 -4.99E-08 -3.85E-09 -7.40E-07 -6.56E-07  1.65E-09  6.85E-07 -3.56E-06  1.74E-06
         -1.86E-06  7.39E-08
 
 OM11
+       -3.72E-06  4.90E-06  3.16E-07  2.75E-03 -2.19E-04 -9.70E-05 -4.54E-05 -2.20E-04  4.04E-04  2.22E-04 -1.21E-03 -5.81E-04
          9.06E-04 -2.31E-06  4.90E-02
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM22
+        4.88E-06 -1.40E-06 -1.78E-06 -8.85E-04  2.65E-04  1.45E-05 -2.98E-03 -2.37E-04 -4.11E-03 -7.73E-05  1.47E-03 -4.87E-04
          2.83E-04  1.12E-05 -4.79E-03 ......... .........  1.79E-01
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 OM33
+        7.46E-06 -9.47E-06 -4.77E-07 -4.70E-04  2.17E-04  2.68E-04 -2.08E-03  2.48E-03 -3.58E-03 -5.14E-05 -7.44E-04  1.21E-03
         -1.90E-03  4.99E-06 -3.88E-03 ......... .........  1.07E-01 .........  1.70E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 TH 1
+        5.18E-04
 
 TH 2
+        1.05E-01  4.47E-04
 
 TH 3
+       -2.09E-03 -1.92E-03  2.35E-04
 
 TH 4
+       -9.46E-02 -1.77E-01 -2.53E-02  2.78E-02
 
 TH 5
+       -1.35E-02 -3.57E-02 -6.71E-03 -2.16E-03  4.60E-03
 
 TH 6
+        2.49E-02 -5.75E-02 -1.17E-02 -2.38E-02  2.20E-01  2.56E-03
 
 TH 7
+       -1.40E-01 -1.30E-02  1.19E-02  1.32E-01  9.48E-03  3.62E-02  2.78E-02
 
 TH 8
+        9.25E-02  1.63E-02 -1.37E-02  1.19E-02 -6.38E-02 -7.85E-02  2.55E-02  3.69E-02
 
 TH 9
+       -2.66E-02 -8.89E-03 -2.82E-03  5.04E-03 -1.15E-01  1.00E-04  4.91E-02 -8.31E-02  3.61E-02
 
 TH10
+        2.15E-02  3.70E-02 -3.62E-02  1.07E-02 -4.14E-02 -1.05E-01 -2.91E-02  7.44E-02  1.56E-03  4.40E-02
 
 TH11
+       -2.68E-02  3.84E-02  4.83E-02  3.55E-02  9.45E-03  1.85E-03 -2.39E-02 -4.51E-02  4.11E-05  6.00E-02  1.17E-01
 
 TH12
+        2.92E-02 -1.47E-02 -3.92E-06 -5.16E-02 -1.85E-02 -7.06E-03 -2.40E-02  9.22E-02 -1.50E-04  4.75E-02  5.00E-02  1.18E-01
 
 TH13
+       -4.32E-02 -1.22E-03  5.92E-04  7.30E-02  9.87E-03 -1.92E-01 -2.60E-02 -7.32E-02 -5.01E-04 -1.96E-02  1.27E-01 -5.76E-02
          5.30E-02
 
 TH14
+        8.17E-03 -3.52E-02 -8.75E-06 -1.22E-01 -3.99E-02 -5.53E-03 -9.80E-02 -6.55E-02  1.68E-04  5.73E-02 -1.12E-01  5.42E-02
         -1.29E-01  2.72E-04
 
 OM11
+       -3.25E-02  4.96E-02  6.06E-03  4.47E-01 -2.15E-01 -1.71E-01 -7.39E-03 -2.69E-02  5.05E-02  2.28E-02 -4.67E-02 -2.23E-02
          7.72E-02 -3.83E-02  2.21E-01
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM22
+        2.22E-02 -7.41E-03 -1.79E-02 -7.51E-02  1.36E-01  1.34E-02 -2.53E-01 -1.52E-02 -2.68E-01 -4.15E-03  2.96E-02 -9.77E-03
          1.26E-02  9.73E-02 -5.11E-02 ......... .........  4.24E-01
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 OM33
+        3.49E-02 -5.13E-02 -4.92E-03 -4.10E-02  1.14E-01  2.53E-01 -1.82E-01  1.63E-01 -2.40E-01 -2.83E-03 -1.54E-02  2.48E-02
         -8.69E-02  4.45E-02 -4.25E-02 ......... .........  6.14E-01 .........  4.12E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 TH 1
+        3.91E+06
 
 TH 2
+       -4.23E+05  5.40E+06
 
 TH 3
+       -1.39E+04  1.11E+05  1.82E+07
 
 TH 4
+        4.06E+03  2.23E+04  6.70E+03  1.81E+03
 
 TH 5
+        5.34E+03 -3.41E+03 -1.34E+03 -9.25E+02  5.37E+04
 
 TH 6
+       -3.44E+04  1.91E+04  3.74E+04 -7.60E+02 -1.82E+04  1.94E+05
 
 TH 7
+        1.02E+04 -2.27E+03 -2.37E+03 -1.99E+02 -1.43E+02 -7.02E+02  1.46E+03
 
 TH 8
+       -5.08E+03 -1.44E+03  1.63E+03 -5.35E+01  3.86E+02  1.97E+03 -4.52E+01  8.29E+02
 
 TH 9
+        9.92E+02  1.48E+03  1.13E+03  3.59E+01  5.47E+02 -6.22E+02  1.88E+01  5.30E+01  8.47E+02
 
 TH10
+       -5.92E+02 -1.74E+03  3.95E+03 -2.07E+01  4.00E+01  9.95E+02  1.81E+01 -3.91E+01 -3.99E+00  5.34E+02
 
 TH11
+        4.75E+02 -1.10E+03 -2.09E+03 -2.53E+01  4.69E+01 -1.17E+02  1.11E+01  1.14E+01 -1.67E+00 -1.48E+01  7.69E+01
 
 TH12
+       -2.74E+02  5.83E+02  1.09E+02  1.69E+01  5.61E+00  4.20E+01  4.70E+00 -2.07E+01 -8.45E-01 -6.30E+00 -5.16E+00  7.38E+01
 
 TH13
+        9.47E+02  6.44E+02  6.49E+02 -2.43E+01 -2.12E+02  1.53E+03  2.59E+01  3.78E+01  1.95E+00  1.50E+01 -1.98E+01  8.00E+00
          3.89E+02
 
 TH14
+        7.40E+04  4.66E+05 -9.53E+04  1.46E+04  4.03E+04  2.26E+04  9.27E+03  8.72E+03 -1.35E+03 -5.71E+03  3.59E+03 -1.70E+03
          8.79E+03  1.46E+07
 
 OM11
+        5.19E+01 -1.81E+03 -5.12E+02 -1.09E+02  2.54E+02  3.08E+02  1.38E+01  1.11E+01 -5.70E+00 -8.83E-02  4.02E+00 -4.24E-01
         -3.63E+00 -3.79E+01  2.87E+01
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM22
+        3.78E+01 -1.58E+02  3.61E+02  3.08E+00 -7.38E+01  2.82E+02  1.87E+01  1.43E+01  1.38E+01  1.53E+00 -1.32E+00  7.60E-01
         -2.20E+00 -9.53E+02  4.95E-01 ......... .........  1.03E+01
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      OM11      OM12      OM13      OM22      OM23      OM33  
 
 OM33
+        7.31E+01  4.22E+02 -2.47E+02  2.19E+00  1.02E+01 -4.86E+02  7.63E+00 -2.35E+01  8.92E+00 -1.52E+00  8.66E-01 -5.39E-01
          2.87E+00  2.31E+02 -1.00E+00 ......... ......... -6.55E+00 .........  1.14E+01
 
 Elapsed finaloutput time in seconds:    24.27
 #CPUT: Total CPU Time in Seconds,      524.391
Stop Time: 
03/03/2022 
18:15
```

[Back](../rti_markov_main)

[Home](../../model-library.github.io/)
