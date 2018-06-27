# MKL (Matlab)

Multiple kernel learning refers to a set of machine learning methods that use a predefined set of kernels and learn an optimal linear or non-linear combination of kernels as part of the algorithm. Reasons to use multiple kernel learning include a) the ability to select for an optimal kernel and parameters from a larger set of kernels, reducing bias due to kernel selection while allowing for more automated machine learning methods, and b) combining data from different sources (e.g. sound and images from a video) that have different notions of similarity and thus require different kernels. Instead of creating a new kernel, multiple kernel algorithms can be used to combine kernels already established for each individual data source.

## Database

Synthetic left ventricular myocardial velocities that emulate the span of cardiac abnormalities that may be observed in a HFPEF population, ranging from completely normal subjects (Group 1) to subjects with a severely impaired cardiac function (Group 5). Four features of the velocity traces, extracted from physiological knowledge about impaired cardiac function, have been modified to create the synthetic curves, namely:

1) Diminished systolic peak velocity 
2) Delay of the systolic peak velocity
3) Appearance of a post-systolic event
4) Fusion during diastole of the negative peaks corresponding to the left ventricular suction and the atrial contraction. 

Twenty subjects have been created for each group, making a total of 100 subjects. The velocity curves are split in 4 segments, as depicted in the figure above. Each of these segments will be a feature to be used as input to the MKL algorithm. 


| Segment number   |  Number of samples per segment  |
| :--------------: | :-----------------------------: |
| #1               | 500                             |
| #2               | 500                             |
| #3               | 250                             |
| #4               | 2000                            |

## Code

Clean version of the unsupervised Multiple Kernel Learning for dimensionality reduction code. A few remarks: 

### Requirements

It is necessary to install the CVX toolbox to be able to run the algorithm. See the details in the CVX web.

### Execution

The C scripts (computeENERGY.c, computeSWA.c and computeSWB.c) need to be compiled in Matlab. To do so just write in the Matlab command line: 

```javascript
mex computeENERGY
```
```javascript
mex computeSWA
```
```javascript
mex computeSWB
```

If any problem occurs during compilation, try:

```javascript
 mex -setup
```

to change the configuration.

After this is done, run **“Launch_MKL.m”** function

## Demo


