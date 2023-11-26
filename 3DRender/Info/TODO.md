# Notes
- Video Command: ffmpeg -framerate 24 -i %d.bmp -c:v libx264rgb -r 24 output.mp4

# ToDo:
- [ ] Make canvas point scalar not a magic number
- [ ] Fix memory leaks

# ToDone:  
- [X] Add wireframe rendering to the 3D renderer so it utilise movement controls
- [X] Add raytracing
    - [X] Ray-triangle intersection
    - [X] Validate intersection
    - [X] Occlusion Shadows
    - [X] Sort out acne
    - [X] Fix colours
    - [X] Fix orientation
    
- [X] Add controls for switching between wireframe, raytraced and rasterized
- [X] Raytracing desperatly needs optimisation, potentially triangles need to only store vertex indices, remove pass by values, prevent copying where possible
- [X] Fix gourand and phong shaders
- [X] Fix ray origin to be the actual position of the camera
- [X] Add Dieletrics
- [X] Add mesh property transitions
- [X] Add soft shadows
- [X] Experiment with methods for shadows behind dieletrics
- [X] Add soft shadows
- [X] Fix texture mapping
- [X] Add mormal mapping aligned with textures
- [X] Parallelize Renderer
- [X] Add Bezier curves
- [X] Add stanford bunny
- [X] Add mesh rotations
