# Econstack (v0) 


## Context
Set of functions written mostly between 2000-03 for my research work during my PhD. A number of standard functions are developed and extended from [Numerical Recipes for C](http://numerical.recipes/). 

Written as a researcher long ago, I have posted the code mostly as is and have formatted some of the header files to make the functions more readable. From 2005-06, I began to re-write (re-factor) large chunks of this library as EconStack (also posted) but remains an unfinished *work-in-progress*.

#### Code modules and Style
Code style is bad as all was written when I was a graduate student, and reflects the weak systems and engineering background of researchers.  

## Research Use Case
The primary use case is for model estimation using [simulated method of moments](https://en.wikipedia.org/wiki/Method_of_simulated_moments). The original Econometrica paper by Daniel McFadden is [here](https://www.jstor.org/stable/1913621) adapted to model [calibration](http://paulgomme.github.io/calibration.pdf). 

In much of my research work, I need to solve a model up to millions of times in a search of a large parameter space. The model is a typically an economic model with hundreds to millions of heterogenous economic agents (e.g. workers differentiated by their age, education endowment, and unobserved ability to learn). Finding a "solution" to each model required solving a (large) system of non-linear equations, a modest and manageble computational intensive cost if you have a small system of n=10 or n=20 equations. However it becomes very costly if you have to use simulated method of moments which would requires solving it thousands to millions of times. 

One of the key computation problems to be reckoned with is that there is no guarantee that given an initial guess for the solution, that the algorthms (all variants of Newton's method) will converge to a solution. This requires gluing together a series of perturbations from a known solution and exploiting properties of local continuity.  

### DSGE (Dynamic Stochastic General Equilibrium) Economic Models 
One attribute of most DSGE models is that a closed solution to equations describing the economic equilibrium does not exist. Hence the solution must be computed. Below is a brief overview of the computation of solving one class of "complex" DSGE models. 

Solving the model essentially consists of finding a solution to a system of equations which corresponds to a steady state of a dynamic general equilibrium model. In my use case, such models had up to tens of millions of heterogeneous agents. It can be proven that the set of solutions to each agent's (constrained) optimization problem is isomorphic to a solution to a system of (non-linear) equations. In this way, it can be proven that, for a convex set in the parameter space, that a solution exists for the set of agents' optimization problems and therefore for a set of non-linear. 

### Estimation 
Estimation of the model using SMM (simulated method of moments) means choosing *N(p)* parameters to best fit *N(m)* moments of data. This means solving a optimization problem with *N(p)* decision variables. I exploit some mathematical properties of the models so that I can deploy a variant of Newton's Method in some appropriately constructed metric space. This then requires solving for a Hessian of dimension *N(p)* at each iteration of the optimization algorithm which means solving the model *3N(p)^2* times. 

### Local Continuity
**Fmin** and **calibration** both allow the user to exploit a local continuity by allowing for an initial guess. In most applications, it is almost certainly necessary to start with a solution *x(p0)* to a known model parameterization *p0*. For an arbitrary parameterization *p*, starting from arbitrary init guess for solution *x*, the algo may not converge to solution *x(p)*. In practice, this is the case most of the time without a "good" init guess. For this reason, we start from a known solution to a known parameterization. We limit the size of the next step in the outer Newton method for the estimation so that *p(n+1)* is not too far from *p(n)*


## License
This library is licensed under the MIT license: 

Copyright 2000-2017, Ronald Leung

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE