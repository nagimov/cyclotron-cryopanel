## Governing Equations 

`d(MU)/dt = m_in * H_in - m_out * H_out + \sum Q`


## Matrix structure 


heateq.m generates algebraic equations `heateq(T)` that are set-up as per the following matrix. 
Then solver.m looks to see when `F=heateq(T)=0`

**********************************************************************************************
```
[N2_t0_input, N2_t0_slice1, ..., N2_t0_slice3, He_t0_input, ..., He_t0_slice3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,     ...     ]
[N2_ti_input, N2_ti_slice1, ..., N2_ti_slice3, He_ti_input, ..., He_ti_slice3]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,      ...    ]
[N2_tf_input, N2_tf_slice1, ..., N2_tf_slice3, He_tf_input, ..., He_tf_slice3]
```
**********************************************************************************************
