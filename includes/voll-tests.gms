ONLY_LINES   = 0;
ONLY_STORAGE = 0;

x.lo(i,j)$(a(i,j)) = 0;  
x.up(i,j)$(a(i,j)) = 1;
f.lo(i,j,t)$(a(i,j)) = -Fmax(i,j);  
f.up(i,j,t)$(a(i,j)) = Fmax(i,j);

Set s "VOLL scenarios" / v1*v8 /;

Parameter VOLL_s(s);
VOLL_s("v1") = 1e2;
VOLL_s("v2") = 1e3;
VOLL_s("v3") = 1e4;
VOLL_s("v4") = 1e5;
VOLL_s("v5") = 1e6;
VOLL_s("v6") = 1e7;
VOLL_s("v7") = 1e8;
VOLL_s("v8") = 1e9;

Parameter
  shed_total(s)     "Σ L(i,t)"
  curt_total(s)     "Σ W(i,t)"
  lines_cost(s)     "Line CAPEX"
  store_cost(s)     "Storage energy CAPEX"
  gen_cost(s)       "Generation operating cost"
  voll_cost(s)      "Shed penalty"
  curt_cost(s)      "Curtailment penalty"
  of_val(s)         "Objective value (OF)"
  Etot(s)           "Total storage energy built [MWh]"
  nb_lines(s)       "# of lines built"
  built_line_s(s,i,j) "x(i,j)"
  E_s(s,i)            "E(i)"
  shed_t_s(s,t)       "Σ_i L(i,t)";

Parameter sweep_summary(s,*)  "Params in function of VOLL";

!! LOOP
loop(s,

  !! set scenario VOLL
  VOLL = VOLL_s(s);

  !! solve model
  solve plan using mip minimizing OF;

  !! aggregate metrics
  shed_total(s) = sum((i,t), L.l(i,t));
  curt_total(s) = sum((i,t), W.l(i,t));
  lines_cost(s) = sum((i,j)$(ord(i)<ord(j) and und(i,j)), CL*dist(i,j)*x.l(i,j));
  store_cost(s) = sum(i, CS*E.l(i));
  gen_cost(s)   = sum((i,t), cg(i)*g.l(i,t));
  voll_cost(s)  = sum((i,t), VOLL*L.l(i,t));
  curt_cost(s)  = sum((i,t), PCURT*W.l(i,t));
  of_val(s)     = OF.l;

  !! investment indicators
  built_line_s(s,i,j)$und(i,j) = round(x.l(i,j));
  nb_lines(s) = sum((i,j)$(ord(i)<ord(j) and und(i,j)), round(x.l(i,j)));
  E_s(s,i)    = E.l(i);
  Etot(s)     = sum(i, E.l(i));

  !! time profile of shedding (optional)
  shed_t_s(s,t) = sum(i, L.l(i,t));

);
!! END LOOP

!! summary table fill
*sweep_summary(s,"VOLL [$ /MWh]")             = VOLL_s(s);
*sweep_summary(s,"Σ L(i,t)")                  = shed_total(s);
*sweep_summary(s,"Σ W(i,t)")                  = curt_total(s);
*sweep_summary(s,"# Lines built")             = nb_lines(s);
*sweep_summary(s,"Σ E(i) [MWh]")              = Etot(s);
*sweep_summary(s,"Line CAPEX [$]")            = lines_cost(s);
*sweep_summary(s,"Storage CAPEX [$]")         = store_cost(s);
*sweep_summary(s,"Generation cost [$]")       = gen_cost(s);
*sweep_summary(s,"Shedding penalty [$]")      = voll_cost(s);
*sweep_summary(s,"Curtailment penalty [$]")   = curt_cost(s);
*sweep_summary(s,"Objective function (OF) [$]")    = of_val(s);

sweep_summary(s,"VOLL")        = VOLL_s(s);
sweep_summary(s,"ΣL")          = shed_total(s);
sweep_summary(s,"ΣW")          = curt_total(s);
sweep_summary(s,"#Ln")         = nb_lines(s);
sweep_summary(s,"ΣE[MWh]")     = Etot(s);
sweep_summary(s,"Line$")       = lines_cost(s);
sweep_summary(s,"Stor$")       = store_cost(s);
sweep_summary(s,"Gen$")        = gen_cost(s);
sweep_summary(s,"VOLL$")       = voll_cost(s);
sweep_summary(s,"Curt$")       = curt_cost(s);
sweep_summary(s,"OF$")         = of_val(s);


display VOLL_s, built_line_s, E_s, shed_t_s;
display sweep_summary;