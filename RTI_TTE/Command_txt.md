# Time-to-first PCR-confirmed COVID-19

## Commands used

<a download="Command_TTE_RTI_BCG_only" href="./Command_TTE_RTI_BCG_only.txt">Download here</a>

```
; Execution
execute Execution_code_TTE_RTI_BCG_only.mod

; Creating simulation model
update_inits Execution_code_TTE_RTI_BCG_only.mod -output_model=Simulation_code_TTE_RTI_BCG_only.mod -flip_comments

; Simulate data
execute Simulation_code_TTE_RTI_BCG_only.mod -clean=0
```

[Back](../rti_tte_main)

[Home](../../model-library.github.io/)