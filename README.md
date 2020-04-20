## GeomScale Challenges 


![logo](https://user-images.githubusercontent.com/60981694/77258473-76ed4c00-6c83-11ea-9740-c2eba0b008f8.png)
![logo2](https://user-images.githubusercontent.com/60981694/77258157-24923800-6c4f-11ea-8258-6539ca614569.png)

:rocket:Geomscale 2020 gsoc Randomised LP Solver Challenges

•	**Easy Task:**   
  This challenge aims at ensuring the successful compilation and run of Volesti package.   
  To test the installation of the package I sampled with sample_points()   
  from a 100 dimensional hypercube and visualized the result using the R extension.  
  
The resulting plot is the following:
![plotpng](https://user-images.githubusercontent.com/60981694/77974408-bb21c180-72ff-11ea-948a-79bdf6edaf15.png)

•	**Medium Task**  
This challenge includes the implementation of randomized cutting plane method   
for Linear Programming problems according to [1].   
Benchmarks reporting run-time and performance of the algorithm were also conducted.  

**Result:** a first naive implementation of the algorithm gives really encouraging results. Heuristics and further improvements can be added in order to make the algorithm more competitive against current state-of-the-art linear programming solvers.

 A sample of run-time benchmarks is provided below:      
![plotimg2](https://user-images.githubusercontent.com/60981694/78061132-6af53e80-7395-11ea-8a02-41939f7d4e94.png)  

 
*References:*  
**[1]** Dabbene, Fabrizio, Pavel S. Shcherbakov, Boris T. Polyak. *A randomized cutting plane method with probabilistic geometric convergence* SIAM Journal on Optimization 20.6 (2010): 3185-3207.
