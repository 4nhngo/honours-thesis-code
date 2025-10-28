display "LINE BUILT ============================================== ";
Parameter built_line(i,j);
built_line(i,j)$und(i,j) = round(x.l(i,j));
display built_line;

display "STORAGE BUILT ============================================== ";
Parameter built_store(i);
built_store(i) = round(E.l(i));
display built_store;

*----- Cost breakdown ------------------------------------------------------
display "COST BREAKDOWN ============================================== ";
Scalar
   cost_lines   "Line CAPEX"
   cost_store   "Storage energy CAPEX"
   cost_voll    "Shed penalty"
   cost_curt    "Curtailment penalty"
   cost_gen     "Generation operating cost"
   cost_total   "Objective value (OF)";
cost_lines = sum((i,j)$(ord(i)<ord(j) and und(i,j)), CL*dist(i,j)*x.l(i,j));
cost_store = sum(i, CS*E.l(i));
cost_voll  = sum((i,t), VOLL * L.l(i,t));
cost_curt  = sum((i,t), PCURT * W.l(i,t));
cost_gen   = sum((i,t), cg(i) * g.l(i,t));
cost_total = OF.l;
display cost_lines, cost_store, cost_voll, cost_curt, cost_gen, cost_total;

*----- Node operations -----------------------------------------------------
display "OPERATIONS SUMMARY - NODE LEVEL ===================================== ";
Parameter
   flows(i,j,t)   "flow [MW]"
   gen(i,t)       "generation [MW]"
   dem(i,t)       "demand [MW]"
   shed(i,t)      "load shedding [MW]"
   curt(i,t)      "generation rejection [MW]"
   charge(i,t)    "charging [MW]"
   discharge(i,t) "discharging [MW]"
   soc(i)       "initial SOC[MWh]";
flows(i,j,t)   = f.l(i,j,t);
gen(i,t)       = g.l(i,t);
dem(i,t)       = d(i,t);
shed(i,t)      = L.l(i,t);
curt(i,t)      = W.l(i,t);
charge(i,t)    = C.l(i,t);
discharge(i,t) = R.l(i,t);
soc(i)       = h0.l(i);
display dem, gen, shed, curt, charge, discharge, soc, flows;

*----- System operations ---------------------------------------------------
display "OPERATIONS SUMMARY - SYSTEM LEVEL ===================================== ";
Parameter
   gen_t(t)    "system generation by time [MW]"
   dem_t(t)    "system demand by time [MW]"
   shed_t(t)   "system shed by time [MW]"
   curt_t(t)   "system curtailment by time [MW]"
   inflow(i,t) "sum flows into i [MW]"
   outflow(i,t) "sum flows out of i [MW]"
   flow_node(i,t) "sum flows at node i";
gen_t(t)  = sum(i, gen(i,t));
dem_t(t)  = sum(i, dem(i,t));
shed_t(t) = sum(i, shed(i,t));
curt_t(t) = sum(i, curt(i,t));
inflow(i,t)  = sum(j$a(j,i), f.l(j,i,t));
outflow(i,t) = sum(j$a(i,j), f.l(i,j,t));
flow_node(i,t) = sum(j$a(i,j), f.l(i,j,t));
display dem_t, gen_t, shed_t, curt_t, flow_node;

*----- Storage summary -----------------------------------------------------
display "STORAGE ============================================== ";
Parameter
   h_start(i)       'SOC at START of t1 [MWh] = alpha0*E'
   h_t1_end(i)      'SOC at END of t1 [MWh]'
   h_end(i)         'SOC at END of last period [MWh]'
   R_sum(i)         'sum discharge over horizon [MWh]'
   C_sum(i)         'sum charge over horizon [MWh]'
   Storage_check(i) 'h_start + ΣC - ΣR - h_end (should = 0)';

h_start(i)  = h0.l(i);        !! start of horizon
h_t1_end(i) = sum(t$(ord(t)=1), h.l(i,t)); !! end of t1
h_end(i)    = sum(t$(ord(t)=card(t)), h.l(i,t));     !! end of last period
R_sum(i)    = sum(t, R.l(i,t));
C_sum(i)    = sum(t, C.l(i,t));
Storage_check(i) = h_start(i) + C_sum(i) - R_sum(i) - h_end(i);

display h_start, h_t1_end, h_end, R_sum, C_sum, Storage_check;

*----- Balance check -------------------------------------------------------
display "BALANCE CHECK ============================================== ";
Parameter balanceResidual(i,t);
balanceResidual(i,t) = gen(i,t) - dem(i,t) + inflow(i,t)
                     - curt(i,t) + shed(i,t) - charge(i,t) + discharge(i,t);
display balanceResidual;