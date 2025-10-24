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
Set t "time periods" / t1*t4 /;

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
    PCURT "penalty for generation rejection [$/MWh]" / 5000000000000 /

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

* ---- Random but structured test data ----
* Node 1 = generator-heavy
* Node 2 = load-heavy
* Node 3 = storage-prone mixed load

d('n1','t1') = 20; d('n1','t2') = 15; d('n1','t3') = 20; d('n1','t4') = 25;
d('n2','t1') = 60; d('n2','t2') = 70; d('n2','t3') = 80; d('n2','t4') = 60;
d('n3','t1') = 40; d('n3','t2') = 30; d('n3','t3') = 20; d('n3','t4') = 50;

gmax('n1','t1') = 120; gmax('n1','t2') = 120; gmax('n1','t3') = 120; gmax('n1','t4') = 120;
gmax('n2','t1') =  10; gmax('n2','t2') =  15; gmax('n2','t3') =  10; gmax('n2','t4') =  15;
gmax('n3','t1') =  20; gmax('n3','t2') =  20; gmax('n3','t3') =  25; gmax('n3','t4') =  25;

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
* 4) EQUATIONS
*======================
Equation
    obj !! Objective = CAPEX_lines + CAPEX_storage + VOLL
    symm(i,j) !! x(i,j)=x(j,i) (one physical line per corridor)
    antisym(i,j,t) !! +30 MW i->j, −30 MW j->i
    capPos(i,j,t) !! f ≤ Fmax·x
    capNeg(i,j,t) !! −Fmax·x ≤ f 
    balance(n,t) !! g − d + Σ inflow − W + L − C + R = 0
    socDyn(n,t) !! h(i,t) = h(i,t−1) + C − R (t=1)
    socBox(n,t) !! h ≤ E
    cRate(n,t) !! C ≤ E
    rRate(n,t) !! R ≤ E
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
*    g (local gen) minus d (local demand)
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
socDyn(i,t)..  h(i,t) =e= h(i,t-1)$(ord(t)>1) + C(i,t) - R(i,t);
* For t>1: today’s energy = yesterday’s energy + charge − discharge.
* For t=1, h(i,t-1) doesn’t exist, so it’s skipped; effectively h(i,1) = C(i,1) − R(i,1)

* SOC & power limits
socBox(i,t)..  h(i,t) =l= E(i);
cRate(i,t)..   C(i,t) =l= E(i); !! Charging capped by E(i)/Δ & Δ=1h 
rRate(i,t)..   R(i,t) =l= E(i); !!Discharging capped by E(i)/Δ & Δ=1h

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

*   f(i,j,t)  "flow [MW]"
*   g(n,t)    "generation [MW]"
*   L(n,t)    "load shedding [MW]"
*   W(n,t)    "generation rejection [MW]"
*   C(n,t)    "charging [MW]"
*   R(n,t)    "discharging [MW]"
*   h(n,t)    "state of charge [MWh]"
*   E(n)      "storage energy [MWh]" ;

* FLOW CHECK
Parameter inflow(i,t) "sum of flows into i", outflow(i,t) "sum of flows out of i";
inflow(i,t)  = sum(j$a(j,i), f.l(j,i,t));
outflow(i,t) = sum(j$a(i,j), f.l(i,j,t));
display inflow, outflow;

Parameter netflow(i,t) "positive = net interchange";
netflow(i,t) = inflow(i,t);

Parameter balanceResidual(i,t) "should be 0";
balanceResidual(i,t) = g.l(i,t) - d(i,t) + netflow(i,t)
                     - W.l(i,t) + L.l(i,t) - C.l(i,t) + R.l(i,t);
display netflow, balanceResidual;

display OF.l, x.l, f.l, g.l, L.l, W.l, C.l, R.l, h.l, E.l;



$offInline
$offEolCom
