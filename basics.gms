$title SETUP

$onText
• Objective = CAPEX_lines + CAPEX_storage + VOLL + Generation Rejection Penalty + Generation Cost
• Power balance: g_it - d_it + Σ_j f_jit  - W_it + L_it - C_it + R_it = 0
• Box constraints for gen, shedding, curtailment (W ≤ g), storage power, SOC
• Antisymmetric flows and line capacities: -Fmax*x_ij ≤ f_ij,t ≤ Fmax*x_ij
$offText
$onEolCom
$onInline

*=====================================================
* 1) SETS
*=====================================================
Set n "nodes" / n1*n3 /;
Alias (n,i,j);
Set t "time periods" / t1*t2 /;


Set und(i,j) "undirected candidate lines (i<j)"
/ n1.n2, n1.n3, n2.n3 /;

Set a(i,j) "directed arcs";
a(i,j) = yes$(und(i,j) or und(j,i));

*=====================================================
* 2) PARAMETERS
*=====================================================
Scalar
    CL     "Line cost [$ per km]"                / 150000 /
    CS     "Storage energy cost [$/MWh]"         / 120 /
    VOLL   "Value of lost load [$/MWh]"          / 1e17 /
    PCURT  "Penalty for generation rejection [$/MWh]" / 500000 /;
Parameter cg(n) "gen marginal cost [$/MWh]";
cg('n1') = 20;  cg('n2') = 60;  cg('n3') = 45;

Table dist(i,j) "distance [km]"
         n1  n2  n3
   n1     0  60  90
   n2     0   0  50
   n3     0   0   0 ;
dist(j,i)$(dist(i,j) and not dist(j,i)) = dist(i,j);

Parameter Fmax(i,j) "line capacity [MW]";
Fmax(i,j) = 100$und(i,j);
Fmax(j,i)$(und(i,j)) = Fmax(i,j);

Parameter
    d(n,t)    "demand [MW]"
    gmax(n,t) "dispatchable generation max [MW]" ;

* ---- Node data ----
d('n1','t1') = 20; d('n1','t2') = 15;
d('n2','t1') = 60; d('n2','t2') = 70;
d('n3','t1') = 40; d('n3','t2') = 30;

gmax('n1','t1') = 120; gmax('n1','t2') = 120;
gmax('n2','t1') =  10; gmax('n2','t2') =  15;
gmax('n3','t1') =  20; gmax('n3','t2') =  20;

* ---- SOC parameters ----
Scalar alpha0 "initial SOC fraction of E(i)" / 0.1 /;
Scalar h_end_frac "final SOC fraction of E(i)" / 0 /;
*=====================================================
* 3) VARIABLES
*=====================================================
Binary Variable x(i,j) "line built (1=yes)";
Variable f(i,j,t) "flow [MW]";

Positive Variables
    g(n,t) "generation [MW]"
    L(n,t) "load shedding [MW]"
    W(n,t) "generation rejection [MW]"
    C(n,t) "charging [MW]"
    R(n,t) "discharging [MW]"
    h(n,t) "state of charge [MWh]"
    E(n)   "storage energy [MWh]";

Variable OF "objective function";

*=====================================================
* 4) SCENARIO RESTRICTIONS
*=====================================================
* --- Scenario toggles (0 = off, 1 = on) ---
Scalar ONLY_LINES   / 0 /,
       ONLY_STORAGE / 1 /;
       
* Only-line analysis
E.up(i)$(ONLY_LINES=1)     = 0;
C.up(i,t)$(ONLY_LINES=1)   = 0;
R.up(i,t)$(ONLY_LINES=1)   = 0;
h.up(i,t)$(ONLY_LINES=1)   = 0;

* Only-storage analysis
x.fx(i,j)$(ONLY_STORAGE=1 and a(i,j))   = 0;
f.fx(i,j,t)$(ONLY_STORAGE=1 and a(i,j)) = 0;

*=====================================================
* 5) EQUATIONS
*=====================================================

*----- Objective ------------------------------------------------------------
Equation obj;
obj..
    OF =e=
        sum((i,j)$(ord(i)<ord(j) and und(i,j)), CL*dist(i,j)*x(i,j))
      + sum(i, CS*E(i))
      + sum((i,t), VOLL*L(i,t))
      + sum((i,t), PCURT*W(i,t))
      + sum((i,t), cg(i)*g(i,t));

*----- Symmetry and antisymmetry -------------------------------------------
Equation symm(i,j), antisym(i,j,t);
symm(i,j)$(a(i,j) and ord(i)<ord(j)).. x(i,j) =e= x(j,i);
antisym(i,j,t)$(a(i,j) and ord(i)<ord(j)).. f(i,j,t) + f(j,i,t) =e= 0;

*----- Line capacities -----------------------------------------------------
Equation capPos(i,j,t), capNeg(i,j,t);
capPos(i,j,t)$(a(i,j))..  f(i,j,t) =l= Fmax(i,j)*x(i,j);
capNeg(i,j,t)$(a(i,j)).. -f(i,j,t) =l= Fmax(i,j)*x(i,j);

*----- Power balance -------------------------------------------------------
Equation balance(i,t);
balance(i,t)..
    g(i,t) - d(i,t) + sum(j$a(j,i), f(j,i,t))
  - W(i,t) + L(i,t) - C(i,t) + R(i,t) =e= 0;

*----- SOC dynamics --------------------------------------------------------
Positive Variable h0(i) "SOC at start of t1 [MWh]";
Equation def_h0(i);
def_h0(i).. h0(i) =e= alpha0*E(i);

Equation socEnd1(i), socEndF(i,t), terminalEnd(i);
socEnd1(i)..              h(i,'t1') =e= h0(i) + C(i,'t1') - R(i,'t1');
socEndF(i,t)$(ord(t)>1).. h(i,t)    =e= h(i,t-1) + C(i,t) - R(i,t);
h.lo(i,t) = 0; !! No negative state of charge

* Terminal SOC condition (connot empty storage?)
Equation terminalEnd(i);
terminalEnd(i)..  sum(t$(ord(t)=card(t)), h(i,t)) =e= h_end_frac * E(i);

*----- Storage power feasibility ------------------------------------------
Equation disEnergyCap0(i), disEnergyCapF(i,t), chgEnergyCap0(i), chgEnergyCapF(i,t);
disEnergyCap0(i)..              R(i,'t1') =l= alpha0*E(i);
disEnergyCapF(i,t)$(ord(t)>1).. R(i,t)    =l= h(i,t-1);
chgEnergyCap0(i)..              C(i,'t1') =l= E(i) - h0(i);
chgEnergyCapF(i,t)$(ord(t)>1).. C(i,t)    =l= E(i) - h(i,t-1);

*----- SOC limits ----------------------------------------------------------
Equation socBox(i,t);
socBox(i,t).. h(i,t) =l= E(i);

*----- Operational box constraints ----------------------------------------
Equation gBox(i,t);
gBox(i,t)..  g(i,t) =l= gmax(i,t); !! 0 ≤ g ≤ gmax = Generators can’t exceed their max.

Equation LBox(i,t);
LBox(i,t)..  L(i,t) =l= d(i,t); !! 0 ≤ L ≤ demand = Can’t shed more load than exists.

Equation WleqG(i,t);
WleqG(i,t).. W(i,t) =l= g(i,t); !! 0 ≤ W ≤ g = Can’t reject more generation than you produce

*=====================================================
* 6) DOMAINING & BOUNDS
*=====================================================
* For any (i,j) that is not a candidate arc, fix the build variable and flow to 0.
x.fx(i,j)$(not a(i,j)) = 0;
f.fx(i,j,t)$(not a(i,j)) = 0;
* Forbid i=i
x.fx(i,i) = 0;  f.fx(i,i,t) = 0;


*=====================================================
* 7) SOLVE
*=====================================================
Model plan / all /;
option optCR = 0;
solve plan using mip minimizing OF;

*=====================================================
* 8) REPORT
*=====================================================
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
display dem, gen, shed, curt, charge, discharge, soc;

*----- System operations ---------------------------------------------------
display "OPERATIONS SUMMARY - SYSTEM LEVEL ===================================== ";
Parameter
   gen_t(t)    "system generation by time [MW]"
   dem_t(t)    "system demand by time [MW]"
   shed_t(t)   "system shed by time [MW]"
   curt_t(t)   "system curtailment by time [MW]"
   inflow(i,t) "sum flows into i [MW]"
   outflow(i,t) "sum flows out of i [MW]";
gen_t(t)  = sum(i, gen(i,t));
dem_t(t)  = sum(i, dem(i,t));
shed_t(t) = sum(i, shed(i,t));
curt_t(t) = sum(i, curt(i,t));
inflow(i,t)  = sum(j$a(j,i), f.l(j,i,t));
outflow(i,t) = sum(j$a(i,j), f.l(i,j,t));
display dem_t, gen_t, shed_t, curt_t, inflow, outflow;

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

$offInline
$offEolCom
