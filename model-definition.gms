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
    VOLL   "Value of lost load [$/MWh]"          / 1e5 /
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
Scalar alpha0 "initial SOC fraction of E(i)" / 0 /;
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
       ONLY_STORAGE / 0 /;
       
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

*$batinclude "includes/basic-report.gms" plan
$batinclude "includes/voll-tests.gms" plan

$offInline
$offEolCom
