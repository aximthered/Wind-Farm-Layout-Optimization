# Wind-Farm-Layout-Optimization
Hybrid Genetic Algorithm and Simulated Annealing implementation to search for the global optimum configuration and thereby to maximize the Annual Energy Production (AEP).

## Introduction 

An optimized layout of wind turbines improves the AEP by 5-15% from that of an unoptimized random arrangement. In our case of 50 turbines, with effective diameters of 400 m, are to be arranged optimally in a square-shaped farm of side 4 km. This optimization problem poses many exciting challenges due to its high dimensionality (100 DOF), Complex multimodality and discontinuous search space. Any analytic or gradient search is impossible for this 100-dimensional optimization problem. Hence, a grid-based hybrid genetic algorithm generates optimum grid solution, and we follow that up with a simulated annealing optimizer to refine the solution towards the global optimum.

![ensembel](/images/ensemble.gif)

## Hybrid Grid based Genetic Algorithm Optimizer (GA)

For reducing the complexity of the problem, the field is subdivided into grids, and the turbines are constrained to be on grid points only. Doing this has two main advantages, 

1. Reduces the complexity of the problem from infinity possible configurations to a combination of the number of grid points and the number of turbines (which is still large).
2. The grid points are chosen so that they implicitly satisfy the proximity constrain on any two turbines. This grid makes the generation of random configurations a mere choosing problem from an extremely time consuming random Markov chain sampling problem. (Note: The rejection rate of the 2D case is prohibitive when approaching the last 20 turbines this cannot be fixed as there are no rejection free sampling algorithms for more than 1D)

Using this grid, we can sample the positions of 50 turbines randomly. These random configurations form the initial population to the GA. Selection, crossover, and mutation operations are applied to the population based on their AEP. In addition to the global search, local search is implemented, which is why the Hybrid label is added. This local search leads to improvement in performance.

![ensembel](/images/GA2.gif)



![ensembel](/images/average_energy_step.png)

## Simulated Annealing Optimizer

Based on the metallurgical method of annealing in which metal at high temperature is slowly cooled to improve its properties. In this case, the virtual temperature is gradually reduced, and the configurations in the neighbourhood of the current optima are sampled. The advantage of this algorithm is its ability to escape local optima which is invaluable in our multimodal optimization problem. 

![ensembel](/images/simann2.gif)

![ensembel](/images/simann_energy_step.png)

## Winner of the north zone in special university edition and 10th on Global leaderboard.

![Certificate](/images/award.PNG)