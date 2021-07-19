include "lex_lesseq.mzn";
include "cumulative.mzn";
include "diffn.mzn";
include "arg_sort.mzn";

%% Inputs
int: w;  % Width of plate
int: n;  % Number of circuits
set of int: CIRCUITS = 1..n; 
array[CIRCUITS, 1..2] of int: dims;  % Dimensions xy for each circuit

% Utilities
int: max_width = min(w, sum(c in CIRCUITS)(dims[c, 1]));
int: max_height = 2*max(max(dim_ver), min_height);
array[CIRCUITS] of float: areas = [1/(dims[c, 1]*dims[c, 2]) | c in CIRCUITS];  % area for each circuit
array[CIRCUITS] of CIRCUITS: per = arg_sort(areas); % Sorting areas from biggest to smallest
array[CIRCUITS] of int: dim_hor = [dims[per[c], 1] | c in CIRCUITS];  % x values dims
array[CIRCUITS] of int: dim_ver = [dims[per[c], 2] | c in CIRCUITS];  % y values dims
int: min_height = floor(sum(c in CIRCUITS)(dims[c, 1]*dims[c, 2]) / w); % Minimum height (objective) computed from minimum area

%% Variables
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol_tmp;  % xy positions of left-bottom vertex for each circuit

% Utilities
array[CIRCUITS] of var int: max_y;  % Max y coordinates for each circuit
array[CIRCUITS] of var 0..max_width-min([dims[c, 1] | c in CIRCUITS]): sol_hor = [sol_tmp[c, 1] | c in CIRCUITS];  % x values sol
array[CIRCUITS] of var 0..max_height-min([dims[c, 2] | c in CIRCUITS]): sol_ver = [sol_tmp[c, 2] | c in CIRCUITS];  % y values sol
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol;  % Ordered solution

% Objective value
var min_height..max_height: obj = max(max_y);  % Minimizing area giving precedence to height

%% Constraints
constraint diffn(sol_hor, sol_ver, dim_hor, dim_ver);  % Avoiding intersections

% Better propagation
constraint cumulative(sol_ver, dim_ver, dim_hor, max_width);  % Cumulative on y
constraint cumulative(sol_hor, dim_hor, dim_ver, obj);  % Cumulative on x

% Avoid symmetries
constraint sol_tmp[1, 1] <= 1+(w-dims[per[1], 1])/2 /\ sol_tmp[1, 2] <= 1+(max_height-dims[per[1], 2])/2;

% Not covered by cumulatives
constraint forall(c in CIRCUITS)(sol_tmp[c, 1] <= max_width - dims[per[c], 1]);  % Constraint on max width

% Better propagation
constraint forall(c in CIRCUITS)(sol_tmp[c, 2] <= obj - dims[per[c], 2]);  % Constraint on max height
           
% Equal shape horizontal & vertical
constraint forall(c1 in CIRCUITS)(symmetry_breaking_constraint(lex_lesseq([sol_hor[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]], [sol_hor[c2] | c2 in reverse(c1..n) where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]])));  % horizontal

constraint forall(c1 in CIRCUITS)(symmetry_breaking_constraint(lex_lesseq([sol_ver[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]], [sol_ver[c2] | c2 in reverse(c1..n) where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]])));  % vertical

% Equal shape only one dimension
constraint forall(c1, c2 in CIRCUITS where c1 < c2 /\ dims[per[c1], 1] == dims[per[c2], 1])(symmetry_breaking_constraint(sol_hor[c1]==sol_hor[c2] -> sol_ver[c1]<=sol_ver[c2]));  % horizontal

constraint forall(c1, c2 in CIRCUITS where c1 < c2 /\ dims[per[c1], 2] == dims[per[c2], 2])(symmetry_breaking_constraint(sol_ver[c1]==sol_ver[c2] -> sol_hor[c1]<=sol_hor[c2]));  % vertical
           
% Utilities
constraint forall(c in CIRCUITS)(max_y[c] == sol[c, 2]+dims[c, 2]);  % Consistence for max_y
constraint forall(c in CIRCUITS)(sol[per[c],1] == sol_tmp[c, 1] /\ sol[per[c],2] == sol_tmp[c, 2]);  % Consistence for ordered solution

%% Search
ann: search_ann;
ann: search_ann1;
search_ann = int_search(sol_hor, input_order, indomain_random);
search_ann1 = int_search(sol_ver, input_order, indomain_random);

%% Solve & output
solve
  :: search_ann
  :: search_ann1
  :: restart_luby(250) % L = 250
  minimize obj;
  
output [
  "Height: ", show(obj), "\n",
  "Solution: ", 
  show(sol_hor),
  show(sol_ver), "\n",
  "Permutation: ", show(per)
];
