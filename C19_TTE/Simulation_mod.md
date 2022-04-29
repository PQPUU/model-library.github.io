# Time-to-first hospitalization

## Simulation model

<a href="./Simulation_code_TTE_hospitalization_BCG_only.mod">Download here</a>

```
$PROBLEM    Time-to-event analysis of hospitalization (all causes)
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
$DATA      TOTALHOSP_ITT_20211203.csv IGNORE=@ IGNORE=(TYPE.EQ.2)

;Sim_start

;IGNORE=(TYPE.EQ.3)

;$SUBROUTINE ADVAN=13 TOL=9
$ABBREVIATED COMRES=7
$SUBROUTINE ADVAN=13 TOL=9
;Sim_end
$MODEL      COMP=(HAZARD)
$PK
 IF(GROUP.EQ.1) RFGROUP = 0
 IF(GROUP.EQ.2) RFGROUP = THETA(2)
 BASE = THETA(1)*EXP(ETA(1))/100 
;Sim_start
;;;;;;;; THE SIMULATION PART FOR TTE SIMULATIONS ;;
 IF (ICALL.EQ.4) THEN                            ; The event time sim $problem
     IF (NEWIND.EQ.0) THEN                       ; Only for the first record
        COM(6)  =  1                            ; Reset simulation ID counter
        COM(4)  =  459                           ; Set max time/censoring time in the dataset
     ENDIF
     IF (NEWIND.EQ.1) THEN                       ; For every new ind except first in dataset
         ICOUNT  =  COM(6) + 1                   ; Update individual counter over simulations
         COM(6)  =  ICOUNT
     ENDIF
     IF (NEWIND.NE.2) THEN                       ; For every new individual
         CALL RANDOM(2,R)
         COM(3)  =  -1                           ; Variable for survival at event time
         COM(2)  =  R                            ; Store the random number
         COM(1)  =  -1                           ; Variable for the event time
         COM(7)  =  0                            ; Individual event counter
     ENDIF
 ENDIF
;;;;;---------MTIME for increasing precision in $DES --------
 IF (NEWIND.NE.2) THEN
 TEMP = 0
 ENDIF
 TEMP = TEMP+0.1
 MTIME(1) = TEMP
 MTDIFF = 1
;Sim_end
$DES
 DADT(1) = BASE*EXP(RFGROUP)
;Sim_start
 SUR = EXP(-A(1))
 IF(COM(2).GT.SUR.AND.COM(1).EQ.-1) THEN ; If event save event time in COM(1)
    COM(1) = T
    COM(3) = SUR
 ENDIF
;Sim_end
$ERROR
;Sim_start
; IF(NEWIND.NE.2) OLDCHZ = 0           ;reset the cumulative hazard
; CHZ = A(1)-OLDCHZ                  ;cumulative hazard from previous time point in data set
; OLDCHZ = A(1)                      ;rename old cumulative hazard
; SUR = EXP(-CHZ)                    ;survival probability
; HAZNOW=BASE*EXP(RFGROUP)             ; rate of event each time pt NB: update with each new model
; IF(DV.EQ.0) Y = SUR                ;censored event (prob of survival)
; IF(DV.NE.0) Y = SUR*HAZNOW         ;prob density function of event
"FIRST
"@CHARACTER(LEN=100)::FMT ! Define FORMAT string for writing dataset
;;;;;; NORMAL TTE MODEL
CHZ = A(1)
SURX = EXP(-CHZ)              ;survival probability
IF (COM(1).GT.COM(4)) THEN      ;IF T > ENDTIME, T=ENDTIME
;;;;;; Check survival again at endtime
IF (COM(2).GT.SURX) THEN
COM(1) = COM(4)
ELSE
COM(1) = -1              ;Integrated too far, reset event
ENDIF
ENDIF
EVT = COM(1)        ; Save Event time
RNM = COM(2)        ; Save random number, just for debugging
ENDTIME = COM(4)          ; Endtime of study
TT = COM(5)            ; Analytic event time
;;;;;;;; ADD RTTE, DV TO OUTPUT, SET DV=0 IF NO EVENT OR CENSORED, DV=1 IF EVENT, RTTE = 1 IF EVENT OR CENSORED
IF (ICALL.EQ.4) THEN          ; Initate DV to 0 (No event)
DV=0
ENDIF
TMDV = 0
IF (EVID.GE.2) THEN         ;Set MDV variable for output
TMDV=1
ENDIF
ICOUNT = COM(6)+(IREP-1)*NINDR
ITER = IREP
;;;;;;;;;; Define the format of the output file
;Sim_end
"LAST
"FMT='(E13.7,9(1XE13.7))' ! The output FORMAT
" ! Write all events
" IF (NEWIND.EQ.0) THEN !Open file at first record
" OPEN (99, FILE = 'simtab.dat', POSITION='APPEND')
" IF (IREP.EQ.1) THEN !Write header for 1st subproblem
" WRITE (99,'(A,9(1XA))') 'ID','DV','TIME','RTTE','SURX','ICOUNT','ITER','RAND','ENDTIME','GROUP'
" ENDIF
" ENDIF
" IF (EVT.NE.-1) THEN !If an EVENT
" DV=1
" RTTE=1
" TMDV=0
" ! Write SIM specific output
" WRITE (99,FMT) ID,DV,EVT,RTTE,COM(3),ICOUNT,ITER,COM(2),ENDTIME,GROUP
" COM(1) = -1 !Reset Event time variable
" COM(2) = 0 !Reset Random variable
" COM(3) = -1 !Reset survival variable
" COM(7) = COM(7) + 1 !Update Event counter
" ELSE IF (LIREC.EQ.NDREC.AND.COM(7).EQ.0) THEN !Right Censoring (if no previous events)
" DV=0
" TMDV=0
" RTTE=1
" TMP=COM(4)
" WRITE (99,FMT) ID,DV,TMP,RTTE,SURX,ICOUNT,ITER,COM(2),ENDTIME,GROUP
" ENDIF
" IF (NDREC.EQ.LIREC.AND.NIREC.EQ.NINDR) THEN ! Last record for last individual
" CLOSE(99) ! Close File pointer
" ENDIF
$THETA  (0,0.0114921) ; BASE/100
$THETA  0.309921 ; RFGROUP
$OMEGA  0  FIX
;Sim_start
;$ESTIMATION MAXEVAL=9990 METHOD=0 LIKE PRINT=1 SIGL=9 NSIG=3
;$COVARIANCE PRINT=E
;$TABLE      ID TIME SUR HAZNOW BASE RFGROUP EVID NOPRINT ONEHEADER
;            FILE=mytab
$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION NSUB=1
;Sim_end
```

[Back](../hospitalization_tte_main)
