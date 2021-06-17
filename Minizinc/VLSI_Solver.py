from minizinc import Instance, Model, Solver
import matplotlib.pyplot as plt
import random
import matplotlib.patches as patches
import numpy as np
import datetime


def plot_result(result):
    dims = result["dims_v"]
    w = result["w_v"]
    ys = np.array(result["sol"])[:,1]
    xs = np.array(result["sol"])[:,0]

    lims = (0, max(w, max([ys[i]+dims[i][1] for i in range(len(ys))])))

    fig1 = plt.figure(figsize=(10, 10))
    ax1 = fig1.add_subplot(111, aspect='equal')
    for i in range(len(xs)):
        ax1.add_patch(patches.Rectangle((xs[i], ys[i]), dims[i][0], dims[i][1], edgecolor='black', facecolor=random.choice(['r', 'g', 'y'])))

    plt.xticks(range(w+1))
    plt.yticks(range(max([ys[i]+dims[i][1] for i in range(len(ys))])+1))
    #plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.show()


def solve_instance(instance_id):
    # Load model from file
    model = Model("./model_2.mzn")
    
    # Find the MiniZinc solver configuration for Gecode
    solver = Solver.lookup("gecode")

    # Assign data
    model.add_file("./instances/ins-" + str(instance_id) + ".dzn")

    # Create an Instance of model for the solver
    instance = Instance(solver, model)
    
    timeout = datetime.timedelta(seconds=300)
    return instance.solve(timeout=timeout, nogood="true")
    

def solve_all(max_instance):
    results = []
    times = []
    for i in range(1, max_instance+1):
        print("Solving: ", i)
        res = solve_instance(i)
        results.append({"obj": res["objective"], "sol": res["sol"], "dims": res["dims_v"], "w": res["w_v"]})
        times.append(res.statistics["solveTime"].total_seconds())
        print("Solved: ", i)
        
    times = np.array(times)
    plt.plot(range(1, max_instance+1), times, 'o')
    plt.show() 


if __name__ == "__main__":
    solve_all(20)

