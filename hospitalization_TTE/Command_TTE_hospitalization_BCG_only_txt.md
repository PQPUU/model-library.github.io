# Time-to-first hospitalization
## Commands used
```![image]
; Execution
execute Executable_TTE_hospitalization_BCG_only.mod

; Creating simulation model
update_inits Executable_TTE_hospitalization_BCG_only.mod -output_model=Simulations_TTE_hospitalization_BCG_only.mod -flip_comments

; Simulate data
execute Simulations_TTE_hospitalization_BCG_only.mod -clean=0
```

[Back](../hosp_main)
