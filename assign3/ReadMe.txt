Assignment #3: Ray tracing



FULL NAME: David Ivan
ID: 1912037993




MANDATORY FEATURES

------------------



<Under "Status" please indicate whether it has been implemented and is

functioning correctly.  If not, please explain the current status.>



Feature:                                 Status: finish? (yes/no)

-------------------------------------    -------------------------

1) Ray tracing triangles                  Yes


2) Ray tracing sphere                     Yes


3) Triangle Phong Shading                 Yes


4) Sphere Phong Shading                   Yes


5) Shadows rays                           Yes


6) Still images                           Yes
   

7) Extra Credit (up to 20 points)
   
I added soft shadow, anti-aliasing, and reflection.

SUPER DUPER IMPORTANT!!!
-If you are testing out the program and are using the table image, don't using reflection (ie, set the number of bounces to 0)
 This is because the program will go slowly, also for some reason the triangles on the table are very reflective, so the will look
 a bit strange
-If you are using the siggraph image, don't use antialiasing or too many shadow rays, the program is simply too slow to handle these things
 since there is a lot of iteration going on.
-The pictures are located in the pic folder that is in the highest level directory, the same directory as the .sln file. 
-Also, in the triangle_ball picture, the triangle DOES cast a shadow on the ground, its just that because I used light attentuation, the area that 
 the shadow is in is a little dark, however you can notice the triangle casting the shadow if you look closely. 

INSTRUCTIONS:

When you run the program, you can set some features
-It will ask you to change the amount of bounces. The range is between 0 and 5. The values will be clamped if you go outside that range
-It will ask you for antialiasing
-It will ask you for amount of shadow rays, between 1 and 100. It will be clamped. The higher than better quality shadows.