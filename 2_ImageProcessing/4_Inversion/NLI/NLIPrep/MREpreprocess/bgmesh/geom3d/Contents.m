% Geometry Toolbox
% Version 0.1 - 2005, March 04.
%
%   Creation, transformations, algorithms and visualization of geometrical
%   3D primitives, such as points, lines, planes, polyhedra, circles and
%   spheres.
%
%
%Primitives creation
%   createLine3d      create a line with various inputs.
%   createPlane       create a plane in parametrized form
%
%Primitive extraction from other primitives
%   circle3dOrigin    return the first point of a 3D circle
%   projPointOnPlane  return the projection of a point on a plane
%   planeNormal       compute the normal to a plane
%   normalizePlane    normalize parametric form of a plane
%   intersectLineSphere   return intersection between a line and a sphere
%   intersectPlaneLine    return intersection between a plane and a line
%   intersectPlanes       return intersection between 2 planes in space
%   intersectPlaneSphere  return intersection between a plane and a sphere
%
%Measurements
%   distancePointPlane    euclidean distance between 3D point and plane
%   cart2sph2         convert cartesian coordinate to spherical coordinate
%   sph2cart2         convert spherical coordinate to cartesian coordinate
%   distancePoints3d  compute euclidean distance between 3D Points
%   sphericalAngle    compute angle on the sphere
%   anglePoints3d     compute angle between 2 3D points
%   dihedralAngle     compute dihedral angle between 2 planes
%   linePosition3d    return position of a 3D point on a 3D line
%   circle3dPosition  return the angular position of a point on a 3D circle
%   planePosition     compute position of a point on a plane
%
%Polyhedra
%   createCube        create a Cube     
%   createOctahedron  create a graph representing this mesh
%   createCubeOctahedron        create a graph representing this mesh
%   createTetrakaidecahedron    create a graph representing this mesh
%   createRhombododecahedron    create a graph representing this mesh   
%
%Manage 3D affine transforms
%   TransformPoint3d  transform a point with a 3D affine transform
%   translation3d     return 4x4 matrix of a 3D translation
%   rotationOx        create 4x4 matrix of a rotation around x-axis
%   rotationOy        return 4x4 matrix of a rotation around x-axis
%   rotationOz        return 4x4 matrix of a rotation around z-axis
%
%Drawing procedures
%   drawPoint3d       draw 3D point on the current axis.
%   drawLine3d        clip and draw line in the current Window
%   drawEdge3d        draw the edge in the current Window
%   drawPlane3d       draw a plane clipped in the current window
%   drawCurve3d       draw a 3D curve specified by a list of points
%   drawCircle3d      draw a 3D circle
%   drawCircleArc3d   draw a 3D circle arc
%   drawSphere        draw a sphere as a mesh
%   drawSphericalTriangle draw a triangle on a sphere
%   drawPolyhedra     draw a Polyhedra defined by vertices and faces
%
%
%   See also : 
%   geom
%
%   -----
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/12/2003.
%

