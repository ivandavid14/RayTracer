This is a ray tracer that has soft shadows, anti-aliasing, and reflection among some other features.
I intend to work on it some since it was quite fun to make. My next goal is to try to make it run faster 
by implementing a better search algorithm for actually detecting collisions as I currently use a linear search
to detect collisions and it is extremely slow for large objects. 

This was a homework assignment, so some of the dlls and pic.h library were given to us so we didn't have to worry
about figuring out how to extract an image from a jpeg image. Its also why the project is called "assign3." 

I plan to swap out the jpeg dlls and use an open-source library that i've used in the past which is smaller and 
easy to use.

specifically this one
https://code.google.com/p/jpeg-compressor/

NOTES
-----
-If you are testing out the program and are using the table image, don't using reflection (ie, set the number of bounces to 0)
 This is because the program will go slowly, also for some reason the triangles on the table are very reflective, so the will look
 a bit strange
-If you are using the siggraph image, don't use antialiasing or too many shadow rays, the program is simply too slow to handle these things
 since there is a lot of iteration going on.
-The pictures are located in the pic folder that is in the highest level directory, the same directory as the .sln file. 
-Also, in the triangle_ball picture, the triangle DOES cast a shadow on the ground, its just that because I used light attentuation, the area that 
 the shadow is in is a little dark, however you can notice the triangle casting the shadow if you look closely. 

INSTRUCTIONS:
-How to run. Take a file from the sample pictures file and copy paste it into the screenfile.txt file in the assign3 folder if you are running it 
from within visual studio. Then just run the program.
-This program is intended to be run on windows and was made in visual studio. If you are trying to run the program from the command line, to run it
 in the command line type "assign3" then the file you wish to use as input (has to be in the same format as the screenfile.txt) and if you want, after
 specifying the input file, you can specify an output jpeg file. 
-You can find sample input files in "samples" in assign3/s

When you run the program, you can set some features
-It will ask you to change the amount of bounces. The range is between 0 and 5. The values will be clamped if you go outside that range
-It will ask you for antialiasing
-It will ask you for amount of shadow rays, between 1 and 100. It will be clamped. The higher than better quality shadows.