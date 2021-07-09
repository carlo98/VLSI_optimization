include "lex_lesseq.mzn";
include "cumulative.mzn";
include "diffn.mzn";
include "arg_sort.mzn";
%include "gecode.mzn";

%% Inputs
int: w;  % Width of plate
int: n;  % Number of circuits
set of int: CIRCUITS = 1..n;
set of int: CIRCUITS_R = 1..2*n; 
array[CIRCUITS, 1..2] of int: dims;  % Dimensions xy for each circuit

% Utilities
int: max_width = min(w, sum(c in CIRCUITS)(dims[c, 1]));
int: max_height = sum(c in CIRCUITS)(dims[c, 2]);
array[CIRCUITS] of float: areas = [1/(dims[c, 1]*dims[c, 2]) | c in CIRCUITS];  % area for each circuit
array[CIRCUITS] of CIRCUITS: per = arg_sort(areas); % Sorting areas from biggest to smallest
int: min_height = floor(sum(c in CIRCUITS)(dims[c, 1]*dims[c, 2]) / w); % Minimum height (objective) computed from minimum area

%% Variables
array[CIRCUITS_R, 1..2] of var 0..max(max_width, max_height): dims_r;  % Dimensions xy for each circuit rotated or not
array[CIRCUITS] of var 0..max_width: dim_hor = [dims_r[per[taken[c]], 1] | c in CIRCUITS];  % x values dims solution
array[CIRCUITS] of var 0..max_height: dim_ver = [dims_r[per[taken[c]], 2] | c in CIRCUITS];  % y values dims solution
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol_tmp;  % xy positions of left-bottom vertex for each circuit
array[CIRCUITS_R, 1..2] of var 0..max(max_width, max_height): sol_r;  % xy positions of left-bottom vertex for each circuit for rotated and original

% Utilities
array[CIRCUITS] of var int: max_y;  % Max y coordinates for each circuit
array[CIRCUITS] of var 0..max_width-min([dims[c, 1] | c in CIRCUITS]): sol_hor = [sol_tmp[c, 1] | c in CIRCUITS];  % x values sol
array[CIRCUITS] of var 0..max_height-min([dims[c, 2] | c in CIRCUITS]): sol_ver = [sol_tmp[c, 2] | c in CIRCUITS];  % y values sol
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol;  % Ordered solution
array[CIRCUITS] of var CIRCUITS_R: taken;  % taken[i] gives circuit i rotated or original

% Objective value
var min_height..max_height: obj = max(max_y);  % Minimizing area giving precedence to height

%% Constraints
constraint diffn(sol_hor, sol_ver, dim_hor, dim_ver);  % Avoiding intersections
constraint cumulative(sol_ver, dim_ver, dim_hor, max_width);  % Cumulative on y
constraint cumulative(sol_hor, dim_hor, dim_ver, obj);  % Cumulative on x

% Implicit
constraint forall(c in CIRCUITS)(sol_tmp[c, 1] <= max_width - dims_r[per[taken[c]], 1]);  % Implicit constraint on max width
constraint forall(c in CIRCUITS)(sol_tmp[c, 2] <= max_height - dims_r[per[taken[c]], 2]);  % Implicit constraint on max height
           
% Equal shape horizontal & vertical
constraint forall(c1 in CIRCUITS)(lex_lesseq([sol_hor[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]], reverse([sol_hor[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]])));  % horizontal

constraint forall(c1 in CIRCUITS)(lex_lesseq([sol_ver[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]], reverse([sol_ver[c2] | c2 in c1..n where [dims[per[c1], 1], dims[per[c1], 2]] == [dims[per[c2], 1], dims[per[c2], 2]]])));  % vertical

constraint forall(c1, c2 in CIRCUITS where c1 <= c2 /\ dims[per[c1], 1] == dims[per[c2], 1])(symmetry_breaking_constraint(sol_hor[c1]==sol_hor[c2] -> sol_ver[c1,]<=sol_ver[c2]));  % horizontal

constraint forall(c1, c2 in CIRCUITS where c1 <= c2 /\ dims[per[c1], 2] == dims[per[c2], 2])(symmetry_breaking_constraint(sol_ver[c1]==sol_ver[c2] -> sol_hor[c1]<=sol_hor[c2]));  % vertical
           
% Utilities
constraint forall(c in CIRCUITS)(sol[per[c],1] == sol_tmp[c, 1] /\ sol[per[c],2] == sol_tmp[c, 2]);  % Consistence for ordered solution
constraint forall(c in CIRCUITS where dims[per[c], 1] == dims[per[c], 2])(taken[c] = c); % Do no rotate squares
constraint forall(c in CIRCUITS)(max_y[c] = sol_tmp[c, 2]+dims_r[per[taken[c]], 2]);  % Consistence for max_y
constraint forall(c in CIRCUITS)(taken[c] = c \/ taken[c] = c+n); % Choosing between rotated or original
constraint forall(c in CIRCUITS)(dims_r[c, 1] = dims[c, 1] /\ 
                                 dims_r[c, 2] = dims[c, 2] /\
                                 dims_r[c+n, 1] = dims[c, 2] /\
                                 dims_r[c+n, 2] = dims[c, 1]);  % Consistence for original and rotated dimensions
constraint forall(c in CIRCUITS)(sol_tmp[c, 1] = sol_r[taken[c], 1] /\ sol_tmp[c, 2] = sol_r[taken[c], 2]); % Consistence for solution

%% Search
ann: search_ann;
ann: search_ann1;
%search_ann = int_search(q, input_order, indomain_min);
%search_ann = int_search(q, first_fail, indomain_min);
search_ann = int_search(sol_hor, input_order, indomain_min);
search_ann1 = int_search(sol_ver, input_order, indomain_min);
%search_ann = int_search(q, input_order, indomain_random);
%search_ann = int_search(sol_tmp, dom_w_deg, indomain_random);

%% Solve & output
solve
  :: search_ann
  :: search_ann1
  :: restart_luby(250) % L = 250
  %:: relax_and_reconstruct(array1d(sol_tmp), 45)
  minimize obj;
  
output [
  "Height: ", show(obj), "\n",
  "Solution: ", 
  show(sol_hor),
  show(sol_ver), "\n",
  "Permutation: ", show(per)
];
