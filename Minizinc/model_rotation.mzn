include "lex_lesseq.mzn";
include "cumulative.mzn";
include "diffn.mzn";
include "gecode.mzn";

%% Inputs
int: w;  % Width of plate
int: n;  % Number of circuits
set of int: CIRCUITS = 1..n;
set of int: CIRCUITS_R = 1..2*n; 
array[CIRCUITS, 1..2] of int: dims;  % Dimensions xy for each circuit

% Utilities
int: max_width = min(w, sum(c in CIRCUITS)(dims[c, 1]));
int: max_height = sum(c in CIRCUITS)(dims[c, 2]);

%% Variables
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol;  % xy positions of left-bottom vertex for each circuit
array[CIRCUITS_R, 1..2] of var 0..max(max_width, max_height): sol_r;  % xy positions of left-bottom vertex for each circuit for rotated and original

% Utilities
array[CIRCUITS] of var int: max_y;  % Max y coordinates for each circuit
array[CIRCUITS] of var int: sol_hor = [sol[c, 1] | c in CIRCUITS];  % x values sol
array[CIRCUITS] of var int: sol_ver = [sol[c, 2] | c in CIRCUITS];  % y values sol
array[CIRCUITS] of var CIRCUITS_R: taken;  % taken[i] gives circuit i rotated or original
array[CIRCUITS_R, 1..2] of var 0..max(max_width, max_height): dims_r;  % Dimensions xy for each circuit rotated or not
array[CIRCUITS] of var 0..max_width: dim_hor = [dims_r[taken[c], 1] | c in CIRCUITS];  % x values dims solution
array[CIRCUITS] of var 0..max_height: dim_ver = [dims_r[taken[c], 2] | c in CIRCUITS];  % y values dims solution
%array[CIRCUITS, 1..2] of var int: dims_v;  % Used for plot in python
%var int: w_v;  % Used for plot in python

% Objective value
var min([dims[c, 2] | c in CIRCUITS])..max_height: obj = max(max_y);  % Minimizing area giving precedence to height

%% Constraints
constraint diffn(sol_hor, sol_ver, dim_hor, dim_ver);  % Avoiding intersections
constraint cumulative(sol_ver, dim_ver, dim_hor, max_width);  % Cumulative on y
constraint cumulative(sol_hor, dim_hor, dim_ver, obj);  % Cumulative on x

% Implicit
constraint forall(c in CIRCUITS)(sol[c, 1] <= max_width - dims_r[taken[c], 1]);  % Implicit constraint on max width
constraint forall(c in CIRCUITS)(sol[c, 2] <= max_height - dims_r[taken[c], 2]);  % Implicit constraint on max height

% Symmetry breaking
constraint lex_lesseq(sol_hor, reverse(sol_hor)) /\  % Horizontal
           lex_lesseq(sol_ver, reverse(sol_ver));  % Vertical
           
% Equal shape horizontal & vertical
constraint forall(c1 in CIRCUITS)(lex_lesseq([sol_hor[c2] | c2 in CIRCUITS where c1 <= c2 /\ {dims[c1, 1], dims[c1, 2]} == {dims[c2, 1], dims[c2, 2]}], reverse([sol_hor[c2] | c2 in CIRCUITS where c1 <= c2 /\ {dims[c1, 1], dims[c1, 2]} == {dims[c2, 1], dims[c2, 2]}])));  % horizontal

constraint forall(c1 in CIRCUITS)(lex_lesseq([sol_ver[c2] | c2 in CIRCUITS where c1 <= c2 /\ {dims[c1, 1], dims[c1, 2]} == {dims[c2, 1], dims[c2, 2]}], reverse([sol_ver[c2] | c2 in CIRCUITS where c1 <= c2 /\ {dims[c1, 1], dims[c1, 2]} == {dims[c2, 1], dims[c2, 2]}])));  % vertical

constraint forall(c in CIRCUITS where dims[c, 1] == dims[c, 2])(taken[c] = c); % Do no rotate squares
           
% Utilities
constraint forall(c in CIRCUITS)(max_y[c] = sol[c, 2]+dims_r[taken[c], 2]);  % Consistence for max_y
constraint forall(c in CIRCUITS)(taken[c] = c \/ taken[c] = c+n); % Chossing between rotated or original
constraint forall(c in CIRCUITS)(dims_r[c, 1] = dims[c, 1] /\ 
                                 dims_r[c, 2] = dims[c, 2] /\
                                 dims_r[c+n, 1] = dims[c, 2] /\
                                 dims_r[c+n, 2] = dims[c, 1]);  % Consistence for original and rotated dimensions
constraint forall(c in CIRCUITS)(sol[c, 1] = sol_r[taken[c], 1] /\ sol[c, 2] = sol_r[taken[c], 2]); % Consistence for solution
%constraint forall(c in CIRCUITS)(dims_v[c, 1] == dims[c, 1] /\ dims_v[c, 2] == dims[c, 2]);
%constraint w_v == w;

%% Search
ann: search_ann;
ann: search_ann1;
%search_ann = int_search(q, input_order, indomain_min);
%search_ann = int_search(q, first_fail, indomain_min);
search_ann = int_search([sol_r[c, 1] | c in CIRCUITS], dom_w_deg, indomain_min);
search_ann1 = int_search([sol_r[c, 2] | c in CIRCUITS], dom_w_deg, indomain_min);
%search_ann = int_search(q, input_order, indomain_random);
%search_ann = int_search(sol_r, dom_w_deg, indomain_random);

%% Solve & output
solve
  :: search_ann
  :: search_ann1
  :: restart_luby(250) % L = 250
  %:: relax_and_reconstruct(array1d(sol_r), 75)
  minimize obj;
  
output [
  "Height: ", show(obj), "\n",
  "Solution: ", 
  show(sol_hor),
  show(sol_ver), "\n",
  "Taken: ",show(taken), "\n",
  "Dim Hor: ",show(dim_hor), "\n",
  "Dim Ver: ", show(dim_ver)
];