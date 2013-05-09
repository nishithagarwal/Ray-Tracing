Nishith Agarwal
ID: 1624-574-357
CSCI 480 HW3


Description:

The program uses Ray Tracing to display objects (triangles and spheres) on the a 640x480 window. It also uses the Phong Illumination Model to calculate the ambient, diffuse and specular component of the object color. It then uses shadow rays to cast shadows on the objects.

Execute Instructions:
> tar -xvf nishith_agarwal_csci480_hw3.zip
> cd pic
> make
> cd ..
> cd assign3
> make
> ./assign3 screenfile1.txt

Extra Credit Features:

1. Displays soft shadows by casting multiple shadow rays towards the light source and taking an average (for diffuse color).

2. Uses anti-aliasing techniques to make the image sharper and remove aliasing issues. The technique used casts multiple rays through each pixel to better approximate the pixel color.

Images:

Pictures 000.jpg - 004.jpg show focus on anti-aliasing
Pictures 005.jpg - 009.jpg focus on soft shadows
