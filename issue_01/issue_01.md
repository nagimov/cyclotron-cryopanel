# heat exchanger configuration

a_in  ------------>------------- a_out

b_out ------------<------------- b_in


# ode45 output matrix configuration 

The solver produces matrix T with solution data. 

Each row contains T at each stage of the heat exchanger, for a fixed t.

For N=3, columns 1-4 contain stream A, columns 5-8 contain stream B.



[N2_t0_input  N2_t0_stage1 ... N2_t0_stage3 He_t0_input ... He_t0_stage3]

[   .................................................................   ]

[N2_ti_input  N2_ti_stage1 ... N2_ti_stage3 He_ti_input ... He_ti_stage3]

[   .................................................................   ]

[N2_tf_input  N2_tf_stage1 ... N2_tf_stage3 He_tf_input ... He_tf_stage3]

Where t0 is the intial time and tf is the final time 

