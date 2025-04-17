http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=8002&objectType=file


%%%%%%%%%%%%%%%%%%%%%%%
normalize() function was not included in the original package.  But it should be a simple one, so I created it in the same folder: normalize.m  (Songbai Ji, 2/6/2006).
%%%%%%%%%%%%%%%%%%%%%%%%%



geom3d
  	Author: 	David Legland
  	Summary: 	functions to create, transform, intersect and draw 3D geometrical shapes
  	MATLAB Release: 	R14SP1
  	Description: 	A library to handles and manipulate 3D geometrical primitives, such as points, lines, planes, rays, 3D circles, spheres.

Library currently allows :

shape generation :
---------------
- line (from 2 points, parallel to another one)
- planes (from 3 point, from a point and a normal, normalize a plane)
- surface of revolution
- various polyhedra (cubeocathedron, octaheron,

shapes intersections & projections :
---------------------
- line and plane
- 2 planes
- line and sphere
- plane and sphere
- project point on plane

Measures between sets :
----------------------
- distances between points, and between points and a set of points
- distance from a point to a 3D line, or to a plane
- position of a point on a 3D circle, a 3D line, or a plane
- angles between 3 points
- spherical angles
- dihedral angle between 2 planes

Transforms :
----------
rotations, translation, and computation of transformed points .

Display :
-------
Draw 3D points, planes, lines, rays, circles, circles arcs, spheres, spherical triangles, polyhedra
Infinite shapes are automatically clipped with the current axis window

A detailed contents is given in file geom3D/Contents.m 