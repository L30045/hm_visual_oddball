function h = plt_cube(center, len_cube, facecolor)

% Define the coordinates of the vertices of the cube
v = [-len_cube -len_cube -len_cube;  % vertex 1
      len_cube -len_cube -len_cube;  % vertex 2
      len_cube  len_cube -len_cube;  % vertex 3
     -len_cube  len_cube -len_cube;  % vertex 4
     -len_cube -len_cube  len_cube;  % vertex 5
      len_cube -len_cube  len_cube;  % vertex 6
      len_cube  len_cube  len_cube;  % vertex 7
     -len_cube  len_cube  len_cube]; % vertex 8
 
% Move to Center 
v = v + reshape(center,1,[]);

% Define the faces of the cube by specifying the vertices
f = [1 2 3 4;  % bottom face
     2 6 7 3;  % front face
     6 5 8 7;  % top face
     5 1 4 8;  % back face
     4 3 7 8;  % right face
     1 5 6 2]; % left face

% Plot the cube
h = patch('Vertices',v,'Faces',f,'FaceColor',facecolor,'EdgeColor','k','FaceAlpha',0.5);
axis equal;  % make the axes have equal scales

end