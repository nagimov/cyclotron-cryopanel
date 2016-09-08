## Matrix structure 


heateq.m generates algebraic equations `heateq(T)` that are set-up as per the following matrix. 
Then solver.m looks to see when `F=heateq(T)=0`

**********************************************************************************************
```
[N2_t0_input, N2_t0_stage1, ..., N2_t0_stage3, He_t0_input, ..., He_t0_stage3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,     ...     ]
[N2_ti_input, N2_ti_stage1, ..., N2_ti_stage3, He_ti_input, ..., He_ti_stage3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,      ...    ]
[N2_tf_input, N2_tf_stage1, ..., N2_tf_stage3, He_tf_input, ..., He_tf_stage3]
```
**********************************************************************************************
