# CellAggregationAnalysis
This code was developed to analyze images of cell clusters formed through aggregation.  A Laplacian of Gaussian (LoG) filter is used to find objects in a fluorescent microscope image.  The number of clusters as well as statistics on the cluster areas are reported as output.

## Dependencies
This code includes the function `circfit.m` from the MATLAB file exchange:
> Izhak Bucher (2020). Circle fit (https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit), MATLAB Central File Exchange. Retrieved May 21, 2020.

The function `DecomposedMexiHat.m` was developed by Leonard Campanello ([@ljcamp1624](https://github.com/ljcamp1624)) in the lab of Wolfgang Losert ([@losertlab](https://github.com/losertlab)).

## Using the Code
### Automatic Analysis
The batch script `clusterFind_TiledImage.m` will analyze folders of .tif format images.  The user is required to specify the location and type of images, as well as the image resolution.  Image analysis parameters are robust to the same type of imaging across different experimental conditions, but need to be changed if the imaging set up or parameters change significantly.

### Manual Correction
The script `clusterFind_manualCorrect.m` allows the user to manually select clusters to remove from the automatic analysis (for example, pieces of dust that were incorrectly detected as cell clusters).  The user specifies the image to be edited and runs the script.  The user then interacts with command line prompts and figures to manually select clusters for deletion.