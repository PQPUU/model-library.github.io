# Respiratory tract infection over time Markov model

## Execution model

<a href="./Execution_code_Markov_RTI_BCGonly.mod">Download here</a>

```

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
$DATA      Simulated_Markov_RTI_BCGonly.csv IGNORE=@
;$DATA      Simulated_Markov_RTI_BCGonly_SIM.csv IGNORE=@
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
$OMEGA  0.09  ; IIV_P10
0.09 ; IIV_P20
0.09 ; IIV_P21
;Sim_start
;
;$SIGMA 0 FIX ; DUMMY for simulation
;

$ESTIMATION METH=COND LAPLACE LIKE MAXEVAL=9999
$COV 
$TABLE      ID EVENT TIME DV P01 P02 P03 P10 P12 P13 P20
            P21 P23 P30 P31 P32 CP00 CP01 CP02 CP03 CP10 CP11 CP12 CP13
			CP20 CP21 CP22 CP23 CP30 CP31 CP32 CP33 T00 T01 T02 T03
            T10 T11 T12 T13 T20 T21 T22 T23 T30 T31 T32 T33 ETAS(1:LAST) 
			NOPRINT ONEHEADER FILE=mytab
;$SIMULATION (12345) (54321 UNIFORM) SUB=50 ONLYSIM
;$TABLE      ID EVENT TIME P01 P02 P03 P10 P12 P13 P20
;            P21 P23 P30 P31 P32 CP00 CP01 CP02 CP03 CP10 CP11 CP12 CP13
;			 CP20 CP21 CP22 CP23 CP30 CP31 CP32 CP33 T00 T01 T02 T03
;            T10 T11 T12 T13 T20 T21 T22 T23 T30 T31 T32 T33 ETAS(1:LAST) 
;            NOPRINT ONEHEADER FILE=mysimtab
;Sim_end
```

[Back](../rti_markov_main)

[Home](../../model-library.github.io/)