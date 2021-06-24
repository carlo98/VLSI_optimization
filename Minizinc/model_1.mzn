include "lex_lesseq.mzn";
include "cumulative.mzn";
include "diffn.mzn";

%% Inputs
int: w;  % Width of plate
int: n;  % Number of circuits
set of int: CIRCUITS = 1..n; 
array[CIRCUITS, 1..2] of int: dims;  % Dimensions xy for each circuit

% Utilities
int: max_width = min(w, sum(c in CIRCUITS)(dims[c, 1]));
int: max_height = sum(c in CIRCUITS)(dims[c, 2]);
array[CIRCUITS] of int: dim_hor = [dims[c, 1] | c in CIRCUITS];  % x values dims
array[CIRCUITS] of int: dim_ver = [dims[c, 2] | c in CIRCUITS];  % y values dims

%% Variables
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol;  % xy positions of left-bottom vertex for each circuit

% Utilities
array[CIRCUITS] of var int: heights;  % Height for each circuit
array[CIRCUITS] of var int: sol_hor = [sol[c, 1] | c in CIRCUITS];  % x values sol
array[CIRCUITS] of var int: sol_ver = [sol[c, 2] | c in CIRCUITS];  % y values sol

% Objective value
var int: obj = max(heights);

%% Constraints
constraint forall(c in CIRCUITS)(heights[c] == sol[c, 2]+dims[c, 2]);  % Consistence for heights
constraint diffn(sol_hor, sol_ver, dim_hor, dim_ver);  % Avoiding intersections
constraint cumulative(sol_ver, dim_ver, dim_hor, max_width);  % Cumulative on y
constraint cumulative(sol_hor, dim_hor, dim_ver, obj);  % Cumulative on x

% Implicit
constraint forall(c in CIRCUITS)(sol[c, 1] <= max_width - dims[c, 1]);  % Implicit constraint on max width
constraint forall(c in CIRCUITS)(sol[c, 2] <= max_height - dims[c, 2]);  % Implicit constraint on max height

%% Search
ann: search_ann;
%search_ann = int_search(q, input_order, indomain_min);
%search_ann = int_search(q, first_fail, indomain_min);
search_ann = int_search(sol, dom_w_deg, indomain_min);
%search_ann = int_search(q, input_order, indomain_random);
%search_ann = int_search(sol, dom_w_deg, indomain_random);

%% Solve & output
solve
  :: search_ann
  :: restart_luby(250) % L = 250
  minimize obj;
  
output [
  "obj: ", show(obj), "\n",
  "Solution: ", 
  show(sol_hor),
  show(sol_ver)
];