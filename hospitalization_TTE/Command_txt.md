# Time-to-first hospitalization

## Commands used

<a download="Command_TTE_hospitalization_BCG_only" href="./Command_TTE_hospitalization_BCG_only.txt">Download here</a>

```![image]
; Execution
execute Executable_TTE_hospitalization_BCG_only.mod

; Creating simulation model
update_inits Executable_TTE_hospitalization_BCG_only.mod -output_model=Simulations_TTE_hospitalization_BCG_only.mod -flip_comments

; Simulate data
execute Simulations_TTE_hospitalization_BCG_only.mod -clean=0
```

[Back](../hospitalization_tte_main)
