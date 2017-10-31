# MathLib 


## Context
Written mostly between 2000-03 for my research work during my PhD. A number of main functions are developed and extended from [Numerical Recipes for C](http://numerical.recipes/). 

Written as a researcher long ago, I have posted the code mostly as is and have formatted some of the header files to make the functions more readable. From 2005-06, I began to re-write (re-factor) large chunks of this library as EconStack (also posted) but remains an unfinished *work-in-progress*.

#### Code modules and Style
Style is bad as all was written when I was a graduate student, and reflects the weak systems and engineering background of researchers.  

## Research Use Case
The primary use case is for model estimation using [simulated method of moments](https://en.wikipedia.org/wiki/Method_of_simulated_moments). The original Econometrica paper by Daniel McFadden is [here](https://www.jstor.org/stable/1913621) adapted to model [calibration](http://paulgomme.github.io/calibration.pdf). 

In my research work, I need to solve a model millions of times in a search in a large parameter space. The model is a typically an economic model with hundreds to millions of heterogenous economic agents (e.g. workers differentiated by their age, education endowment, unobserved ability to learn and ability to retain know-how). Finding a "solution" to each model required solving a (large) system of non-linear equations. 

One of the key computation problems to be reckoned with is that there is no guarantee that given an initial guess for the solution, that the algorthms (all variants of Newton's method) will converge to a solution.  

### DSGE (Dynamic Stochastic General Equilibrium) Economic Models 
One attribute of most DSGE models is that a closed solution to equations describing the economic equilibrium does not exist. Hence the solution must be computed. Below is a brief overview of the computation of solving one class of "complex" DSGE models. 

Solving the model essentially consists of finding a solution to a system of equations which corresponds to a steady state of a dynamic general equilibrium model. In my use case, such models had up to tens of millions of heterogeneous agents. It can be proven that the set of solutions to each agent's (constrained) optimization problem is isomorphic to a solution to a system of (non-linear) equations. In this way, it can be proven that, for a convex set in the parameter space, that a solution exists for the set of agents' optimization problems and therefore for a set of non-linear. 


## License
This library is licensed under the MIT license: 

Copyright 2000-2017, Ronald Leung

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE