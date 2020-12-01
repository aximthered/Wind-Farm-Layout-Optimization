# Wind-Farm-Layout-Optimization
Hybrid Genetic Algorithm and Simulated Annealing implementation to search for the global optimum configuration and thereby to maximize the Annual Energy Production (AEP). 

```
The Shell.ai Hackathon for Sustainable and Affordable Energy kicked-off on September 14th, focusing on a Windfarm 

Layout Optimisation coding challenge. The participants were invited to optimise the placement of 50 wind turbines 

of ‘100 m rotor diameter to help maximise the AEP (Annual Energy Production), each on a hypothetical offshore wind 

farm area. One of the key problems of an unoptimized layout is the combined effect wind turbines can have on the 

wind speed distribution in a windfarm. As a wind turbine extracts energy from incoming wind, it creates a region 

behind it downstream where the wind speed is decreased- this is called a wake region. Note that wind turbines 

automatically orient their rotors, to face incoming wind from any direction. Due to the induced speed deficit, a 

turbine placed inside the wake region of an upstream turbine will naturally generate reduced electrical power. This 

inter-turbine interference is known as a wake effect. An optimal windfarm layout is important to ensure a minimum 

loss of power during this combined wake effect.


The contestants faced challenges such as a high dimensionality, complex multimodality and the discontinuous nature 

of the search space. This made optimising the layout analytics difficult. But, armed with optimisation strategies 

and computer algorithms, around 5000 teams signed up to compete in this challenge.
```

More information at [Hackerearth](https://www.hackerearth.com/challenges/competitive/shell-hackathon/), [Shell](https://www.shell.in/energy-and-innovation/ai-hackathon.html#vanity-aHR0cHM6Ly93d3cuc2hlbGwuaW4vaGFja2F0aG9uLmh0bWw), [Slides](https://www.slideshare.net/Surajk15/windfarm-shell-hackathon)

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