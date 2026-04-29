Follow Fusion 360 API instructions for installation. Fusion 360 may update and make GUI worse with time. This has been updated (May 2026) to be usable with the current update.

General Notes:
Individual source code is held under commands > "CommandName" > entry.py
Control Point Extractor pulls the control points from a 3D control point spline and writes them to an output file.
Airfoil Builder draws any airfoil onto a 2D plane. Input file is a tab delimited, 3 column, list of coordinates.
Propeller Builder is an early version of a propeller builder. It should still do everything, although likely fails with open trailing edge airfoils and other more complex requests.
VP1304 and its respective batch files were built for batch runs of screw propellers. 
NACA66 Builder and Toroidal Builder Batch are the most advanced files as the rest of the files stopped being updated as their use fell out.

NACA66 Builder:
For use only with airfoils whose starting point and ending point are different--the code connects them.
Input file is the output of the MatLab files.
Can build both toroidal and screw propellers.
Screw propellers can have the tip condition that the guide curve is tangent at the tip. Check said box if this is the case. 
Fixed splines helps prevent accidental moving of the airfoils.

To use, select the correct input file. If you are confident that the propeller will generate and you would like the STEP file, check said box, select a folder to place the file in, and give the file a name in "Output File ID". Run the command. If the propeller generates, congrats. If it doesn't, the program will stop at the last command that works, generally this is right before the loft. You'll be able to see all of the airfoils and guide curves. If you perform a surface loft, you'll be able to find the self-intersection points. 

Toroidal Builder Batch:
Input file is the output of the MatLab files.
Select point type as control point if you've figured out how to alter the MatLab files to output control point splines (honestly not too hard).
Fixed splines helps prevent accidental moving of the airfoils.
Connected airfoil TEs (trailing edges) allows this to build airfoils whose starting and ending point are different (checked) and whose are the same (unchecked).
Hub direction can be changed if you've altered the input to be on a different axis.
Expand hub radius by? is a multiplier to make booleaning the surfaces easier.
Rotate to sim? rotates and translates the propeller to the correct location and orientation of our sim.
STEP export? gives you a step file of each propeller you generate in the folder specified.
Save files? saves the files into the same Fusion 360 folder in which the command was run in (see below).

To use, create a Fusion 360 folder you would like to build your batch of propellers in. Save a blank file in this folder, this is needed to house the command. Select the input files you would like to run. Select any other options, although the defaults are quite nice. Any propellers that generate will be saved in the folder. Propellers that don't generate won't be saved. But, you can rebuild them with NACA66 Builder to see where they fail.


If you have any questions, shoot us an email at torus.gemstone@gmail.com, but try to read the source code first. I swear this works really well.