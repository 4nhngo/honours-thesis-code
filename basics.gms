$title SETUP

$onText
• Objective = CAPEX_lines + CAPEX_storage + VOLL + Generation Rejection Penalty + Generation Cost
• Power balance: g_it - d_it + Σ_j f_jit  - W_it + L_it - C_it + R_it = 0
• Box constraints for gen, shedding, curtailment (W ≤ g), storage power, SOC
• Antisymmetric flows and line capacities: -Fmax*x_ij ≤ f_ij,t ≤ Fmax*x_ij
$offText
$onEolCom
$onInline

*======================
* 1) SETS
*======================
Set n "nodes" / n1*n3 /;
Alias (n,i,j);
Set t "time periods" / t1*t2 /;

Set und(i,j) "undirected candidate lines (i<j)"
/ n1.n2, n1.n3, n2.n3 / ;

Set a(i,j) "directed arcs";
a(i,j) = yes$(und(i,j) or und(j,i));

*======================
* 2) PARAMETERS
*======================
Scalar
    CL  "Line cost [$ per km]"        / 150000 /
    CS  "Storage energy cost [$/MWh]" / 120    /
    VOLL "Value of lost load [$/MWh]"  / 100000000000000000  /
    PCURT "penalty for generation rejection [$/MWh]" / 500000 /

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
    
Parameter cg(n) "gen marginal cost [$/MWh]";
cg('n1') = 20;  cg('n2') = 60;  cg('n3') = 45; !! node 1 = generator so cheap

* --- Scenario toggles (0 = off, 1 = on) ---
Scalar ONLY_LINES   / 1 /
       ONLY_STORAGE / 0 /;
       

Scalar alpha0 "initial SOC fraction of E(i)" / 0 /;

* ---- Random data ----
* Node 1 = generator-heavy
* Node 2 = load-heavy
* Node 3 = storage-prone mixed load

d('n1','t1') = 20; d('n1','t2') = 15; 
d('n2','t1') = 60; d('n2','t2') = 70; 
d('n3','t1') = 40; d('n3','t2') = 30; 

gmax('n1','t1') = 120; gmax('n1','t2') = 120; 
gmax('n2','t1') =  10; gmax('n2','t2') =  15; 
gmax('n3','t1') =  20; gmax('n3','t2') =  20; 

*======================
* 3) VARIABLES
*======================
Binary Variable x(i,j) "line built (1=yes)";

* flows must be free (+30 MW i->j, −30 MW j->i)
Variable
    f(i,j,t)  "flow [MW]";
    
Positive Variables
    g(n,t)    "generation [MW]"
    L(n,t)    "load shedding [MW]"
    W(n,t)    "generation rejection [MW]"
    C(n,t)    "charging [MW]"
    R(n,t)    "discharging [MW]"
    h(n,t)    "state of charge [MWh]"
    E(n)      "storage energy [MWh]" ;

Variable OF;

*======================
* ONLY LINES or STORAGE ANALYSIS
*======================
* If ONLY_LINES = 1  → forbid storage everywhere
E.up(i)$(ONLY_LINES=1)     = 0;
C.up(i,t)$(ONLY_LINES=1)   = 0;
R.up(i,t)$(ONLY_LINES=1)   = 0;
h.up(i,t)$(ONLY_LINES=1)   = 0;

* If ONLY_STORAGE = 1 → no lines and no flows
x.fx(i,j)$(ONLY_STORAGE=1 and a(i,j)) = 0;
f.fx(i,j,t)$(ONLY_STORAGE=1 and a(i,j)) = 0;

*======================
* 4) EQUATIONS
*======================
Equation
    obj                 !! Objective = CAPEX_lines + CAPEX_storage + VOLL + Generation Rejection Penalty + Generation Cost
    symm(i,j)           !! x(i,j)=x(j,i)
    antisym(i,j,t)      !! f(i,j,t) + f(j,i,t) = 0
    capPos(i,j,t)       !! f ≤ Fmax·x
    capNeg(i,j,t)       !! −Fmax·x ≤ f
    balance(n,t)        !! nodal balance
    socDyn1(n)          !! t=1: h(i,'t1') = alpha0*E(i)
    socDynF(n,t)        !! t>1: h(i,t+1) = h(i,t) + C - R    
    socBox(n,t)         !! h ≤ E
    chgEnergyCap(n,t)   !! C ≤ E - h                         
    disEnergyCap(n,t)   !! R ≤ h                            
    gBox(n,t)
    LBox(n,t)
    WleqG(n,t);

* min Σ x_ij·CL·d_ij + Σ CS·E_i
* = min line CAPEX + storage CAPEX 
* + VOLL [$ / MWh] * L [MW] = $/Δt & Δt=1h currently
* + cost for generation rejection [$/MWh]
* + penalty on generation cost (varies by node [$/MWh])
obj..
    OF =e=
        sum((i,j)$(ord(i)<ord(j) and und(i,j)), CL*dist(i,j)*x(i,j))
      + sum(i, CS*E(i))
      + sum((i,t), VOLL*L(i,t))
      + sum((i,t), PCURT*W(i,t))
      + sum((i,t), cg(i)*g(i,t));
      
*symm: building i→j equals building j→i (one physical line).
symm(i,j)$(a(i,j) and ord(i)<ord(j))..  x(i,j) =e= x(j,i);
*antisym: flows are opposite: if +30 MW from i→j, then −30 MW from j→i.
antisym(i,j,t)$(a(i,j) and ord(i)<ord(j))..  f(i,j,t) + f(j,i,t) =e= 0;

* Line capacity
* −Fmax·x ≤ f ≤ Fmax·x
capPos(i,j,t)$(a(i,j))..  f(i,j,t) =l= Fmax(i,j)*x(i,j);
capNeg(i,j,t)$(a(i,j)).. -f(i,j,t) =l= Fmax(i,j)*x(i,j);

* Power balance git − dit + Σ flow_in − Wit + Lit − Cit + Rit = 0
* @ node i, time t:
*    g (local gen) - d (local demand)
*   + sum of f(j→i) (net flow into i)
*   - W (gen rejection: throw away gen)
*   + L (shed: we reduce load)
*   - C (charge storage, consumes power)
*   + R (discharge storage, supplies power)
*    = 0.
balance(i,t)..
    g(i,t) - d(i,t) + sum(j$a(j,i), f(j,i,t))
  - W(i,t) + L(i,t) - C(i,t) + R(i,t) =e= 0;

* GENERATION REJECTION: h_it = h_i,t-1 + C_it - R_it
* h(i,t) is START-of-period SOC
socDyn1(i).. h(i,'t1') =e= alpha0*E(i);

Equation socDynF(i,t);
socDynF(i,t)$(ord(t)<card(t))..
    h(i,t+1) =e= h(i,t) + C(i,t) - R(i,t);

* Physical feasibility per period
Equation disEnergyCap(i,t), chgEnergyCap(i,t);
disEnergyCap(i,t)..  R(i,t) =l= h(i,t); !! cannot discharge more than soc
chgEnergyCap(i,t)..  C(i,t) =l= E(i) - h(i,t); !! cannot charge more than availabe space left after e_max-soc

* SOC within installed energy
socBox(i,t)..  h(i,t) =l= E(i);

* Box constraints
gBox(i,t)..    g(i,t) =l= gmax(i,t); !! 0 ≤ g ≤ gmax = Generators can’t exceed their max.
LBox(i,t)..    L(i,t) =l= d(i,t); !! 0 ≤ L ≤ demand = Can’t shed more load than exists.
WleqG(i,t)..   W(i,t) =l= g(i,t); !! 0 ≤ W ≤ g = Can’t reject more generation than you produce.

*======================
* 5) DOMAINING & BOUNDS
*======================
* For any (i,j) that is not a candidate arc, fix the build variable and flow to 0.
x.fx(i,j)$(not a(i,j)) = 0;
f.fx(i,j,t)$(not a(i,j)) = 0;
* Forbid i=i
x.fx(i,i) = 0;  f.fx(i,i,t) = 0;
* No negative state of charge
h.lo(i,t) = 0;

*======================
* 6) SOLVE
*======================
Model plan / all /;
option optCR = 0;
solve plan using mip minimizing OF; !! min Σ x_ij·CL·d_ij + Σ CS·E_i = min line CAPEX + storage CAPEX

*======================
* 7) REPORT
*======================

display '======================= LINE BUILT =======================';
Parameter built_line(i,j)  "1 if line built";
built_line(i,j)  = 0;
built_line(i,j)$und(i,j) = round(x.l(i,j));
display built_line;

display '======================= STORAGE BUILT =======================';
Parameter built_store(i)   "1 if storage built";
built_store(i)   = 0;
built_store(i)            = round(E.l(i));
display built_store;

display '======================= COST BREAKDOWN =======================';

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

display '=================== OPERATIONS SUMMARY - NODE LEVEL ==================';
Parameter
   flows(i,j,t)   "flow [MW]"
   gen(i,t)       "generation [MW]"
   dem(i,t)       "demand [MW]"
   shed(i,t)      "load shedding [MW]"
   curt(i,t)      "generation rejection [MW]"
   charge(i,t)    "charging [MW]"
   discharge(i,t) "discharging [MW]"
   soc(i,t)       "state of charge [MWh]";

flows(i,j,t)   = f.l(i,j,t);
gen(i,t)       = g.l(i,t);
dem(i,t)       = d(i,t);
shed(i,t)      = L.l(i,t);
curt(i,t)      = W.l(i,t);
charge(i,t)    = C.l(i,t);
discharge(i,t) = R.l(i,t);
soc(i,t)       = h.l(i,t);

display  dem,gen, shed, curt, charge, discharge, soc;

display '=================== OPERATIONS SUMMARY - SYSTEM LEVEL ==================';
Parameter
   dem_t(t)    "system demand by time [MW]"
   gen_t(t)    "system generation by time [MW]"
   shed_t(t)   "system shed by time [MW]"
   curt_t(t)   "system curtailment by time [MW]"
   net_injection "net injection at node (g - d - W + L - C + R) [MW]";
   
gen_t(t)  = sum(i, gen(i,t));
dem_t(t)  = sum(i, dem(i,t));
shed_t(t) = sum(i, shed(i,t));
curt_t(t) = sum(i, curt(i,t));

Parameter inflow(i,t) "sum flows into i [MW]", outflow(i,t) "sum flows out of i [MW]";
inflow(i,t)  = sum(j$a(j,i), f.l(j,i,t));
outflow(i,t) = sum(j$a(i,j), f.l(i,j,t));
net_injection(i,t) = g.l(i,t) - d(i,t) - W.l(i,t) + L.l(i,t) - C.l(i,t) + R.l(i,t);

display dem_t, gen_t, shed_t, curt_t, net_injection;

display '======================= STORAGE =======================';
Alias (t,tt);

Parameter
   h_first(i)    "SOC at first period [MWh]"
   h_last(i)     "SOC at last period [MWh]"
   R_sum(i)      "sum discharge over horizon [MWh]"
   C_sum(i)      "sum charge over horizon [MWh]"
   Storage_check(i) "h(t1) + Σ_t C - Σ_t R - h(t_last) (should = 0)";

h_first(i) = sum(t$(ord(t)=1), h.l(i,t));
h_last(i)  = sum(t$(ord(t)=card(t)+1), h.l(i,t));
R_sum(i)   = sum(t, R.l(i,t));
C_sum(i)   = sum(t, C.l(i,t));
Storage_check(i) = h_first(i) - h_last(i) + C_sum(i) - R_sum(i);


display h_first, h_last, R_sum, C_sum, Storage_check;


display '======================= DEMAND SUPPLY DECOMPOSITION =======================';
Parameter
   local_gen_used(i,t)  "g - W - C [MW]  (net local gen used at node)"
   storage_supply(i,t)  "R [MW]"
   demand_served(i,t)   "d - L [MW]"
   contribResidual(i,t) "check: local_gen_used + storage_supply + importMW - demand_served";

local_gen_used(i,t) = g.l(i,t) - W.l(i,t) - C.l(i,t);
storage_supply(i,t) = R.l(i,t);
demand_served(i,t)  = d(i,t) - L.l(i,t);

display local_gen_used, storage_supply, demand_served;

display '======================= BALANCE CHECK =======================';
Parameter balanceResidual(i,t) "KCL residual (should be 0)";
balanceResidual(i,t) = gen(i,t) - dem(i,t) + inflow(i,t)
                     - curt(i,t) + shed(i,t) - charge(i,t) + discharge(i,t);
display balanceResidual;



$offInline
$offEolCom
