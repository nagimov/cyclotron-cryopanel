# heat exchanger configuration

```
a_in  ------------>------------- a_out
b_out ------------<------------- b_in
```

# ode45 output matrix configuration 

The solver produces matrix T with solution data. 
Each row contains T at each stage of the heat exchanger, for a fixed t.
For N=3, columns 1-4 contain stream A, columns 5-8 contain stream B.

**********************************************************************************************
```
[N2_t0_input, N2_t0_slice1, ..., N2_t0_slice3, He_t0_input, ..., He_t0_slice3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,     ...     ]
[N2_ti_input, N2_ti_slice1, ..., N2_ti_slice3, He_ti_input, ..., He_ti_slice3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,      ...    ]
[N2_tf_input, N2_tf_slice1, ..., N2_tf_slice3, He_tf_input, ..., He_tf_slice3]
```
**********************************************************************************************

Where t0 is the intial time and tf is the final time 

