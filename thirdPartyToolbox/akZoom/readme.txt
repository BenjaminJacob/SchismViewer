___ Using akZoom automatically for all 2D-plots in Matlab___
The standard usage of akZoom is "on demand", this means you have to call akZoom explicitly after plotting something like shown in the akZoom_examples-file. However, if you are satisfied with the functionality, you can implement akZoom such that it is called automatically when calling plot, plotyy, etc..

To use akZoom automatically, you just have to include the folder "Wrapper of Matlab functions" to the Matlab search patch (here is how you do that: http://www.mathworks.de/de/help/matlab/ref/path.html). In this folder Matlab will find wrapper-functions for plot, plotyy, etc.. These wrapper-functions call the original functions and automatically add the akZoom call afterwards.

Note: since we overload Matlabs own functions, you will get warnings like the following on Matlab startup:
"Warning: Function ...\plot.m has the same name as a MATLAB builtin. We suggest you rename the function to avoid a potential name conflict."

This can either just be ignored or, if it really bugs you, turned off with the following command:
warning off MATLAB:dispatcher:nameConflict
However, sometimes these kind of warnings can be usefull if you accidentally overload a Matlab function. So I recommend not to turn them off.

___ Going back to default ___
If you do not want to call akZoom automatically, just remove the "Wrapper of Matlab functions" folder from the Matlab search path.

___ Using akZoom in GUIs ___
One important remark for GUI-programmers: When building an executable, Matlab somehow does not like it, when its built-in functions are overwritten. Therefore I recommend to call akZoom explicitly in the m-file of your GUI and remove the "Wrapper of Matlab functions" folder from the search path before building your application.


If you have any questions or remarks just send me an email: alexander.kessel(at)mpq.mpg.de