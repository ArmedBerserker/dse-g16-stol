# TODO: link these to yaml files
AR = 11
MAC = 1.5
S = 23
l_f = 8.9

X_aft_cg = 1 # TBD

X_v = 0.9*l_f
V_v = 0.06275 # from roskam
AR_v = 1.5 # avg from adsee
Sweep_v = 25 # degrees
lambda_v = 0.4 # from adsee

T_over_C_v = 0.15 # from roskam

S_v = (0.36*MAC)/(X_v-1)*S 
b_v = (AR_v*S_v)**0.5
Cr_v = (2*S_v)/(b_v*(1+lambda_v))
Ct_v = lambda_v*Cr_v