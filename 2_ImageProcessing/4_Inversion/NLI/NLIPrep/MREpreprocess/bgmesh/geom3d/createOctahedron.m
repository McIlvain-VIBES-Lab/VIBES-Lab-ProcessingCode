function [nodes, edges, faces] = createOtahedron(varargin)
%createOtahedron create a 3D cube
%
%   c = createOtahedron create a unit cube, as a polyhedra representation.
%   c has the form [n, e, f], where n is a 6*3 array with vertices
%   coordinate, e is a 12*2 array containing indices of neighbour vertices,
%   and f is a 8*3 array containing vertices array of each face.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 10/02/2005.
%

x0 = 0; dx= 1;
y0 = 0; dy= 1;
z0 = 0; dz= 1;

nodes = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];


edges = [1 2;1 4;1 5; 1 6;2 3;2 5;2 6;3 4;3 5;3 6;4 5;4 6];

faces = [1 2 5;2 3 5;3 4 5;4 1 5;1 2 6;2 3 6;3 4 6;4 1 6];

