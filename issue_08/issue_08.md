## Governing Equations 

`d(MU)/dt = m_in * H_in - m_out * H_out + \sum Q`


## Matrix structure 


heateq.m generates algebraic equations `heateq(T)` that are set-up as per the following matrix. 
Then solver.m looks to see when `F=heateq(T)=0`

**********************************************************************************************
```
[N2_t0_input, N2_t0_slice1, ..., N2_t0_sliceN, He_t0_input, ..., He_t0_sliceN]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,     ...     ]
[N2_ti_input, N2_ti_slice1, ..., N2_ti_sliceN, He_ti_input, ..., He_ti_sliceN]
[   ...     ,     ...     , ...,    ...      ,    ...     , ...,      ...    ]
[N2_tf_input, N2_tf_slice1, ..., N2_tf_sliceN, He_tf_input, ..., He_tf_sliceN]
```
**********************************************************************************************


## List of MATLAB solvers

# Solution of non-linear equations 

- 'fsolve': multi variable system
- 'fzero': single variable system
 
# Non-linear Optimization 

- 'fminunc' : find min of a multi variable function
- 'fminsearch' : find min of a multi variable function, derrivative-free method
- 'fminbnd' : constrained optimization, single variable function
- 'fmincon' : constrained optimization, multi variable function
- 'fgoalattain' : multi-goal optimization 
