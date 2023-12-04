# Ray Tracer
This assignment is a ray tracer that can render spheres. It impliments phong
lighting, shadows, and reflections. 

## Summary of Work Done
I was able to successfully implement all the features required for each test case. ADS lighting is calculated 
correctly, shadows are rendered correctly, and reflections are rendered correctly. The program handles edge cases such 
as a sphere intersecting the image plane, or a light source being inside an intersected sphere. As a result, all of my
images render as expected and appear identical to the test keys. 

## Usage
To compile the program, run `make` in the root directory.
Alternatively, run `g++ -std=c++11 raytrace.cpp -o RayTracer.exe` to compile

To run the program, run `./RayTracer.exe <input file>`.

The program will output a ppm file with the name specified in the input file in the folder `output/`.
This folder must exist before running the program.

To run all test cases consecutively, run `make test` in the root directory.

## File Syntax
NEAR &lt;n&gt;  
LEFT &lt;l&gt;  
RIGHT &lt;r&gt;  
BOTTOM &lt;b&gt;  
TOP &lt;t&gt;  
RES &lt;x&gt; &lt;y&gt;  
SPHERE &lt;name&gt; &lt;pos x&gt;= &lt;pos y&gt;= &lt;pos z&gt; &lt;scl x&gt; &lt;scl y&gt; &lt;scl z&gt; &lt;r&gt; &lt;g&gt; &lt;b&gt; &lt;ka&gt; &lt;kd&gt; &lt;ks&gt; &lt;kr&gt; &lt;n&gt;  
// up to 14 additional sphere specifications

LIGHT &lt;name&gt; &lt;pos x&gt; &lt;pos y&gt; &lt;pos z&gt; &lt;ir&gt; &lt;ig&gt; &lt;ib&gt;  
// up to 9 additional light specifications

BACK &lt;r&gt; &lt;g&gt; &lt;b&gt;  
AMBIENT &lt;ir&gt; &lt;ig&gt; &lt;ib&gt;  
OUTPUT &lt;name&gt;  
  


