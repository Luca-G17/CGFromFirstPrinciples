# Notes
- Video Command: ffmpeg -framerate 24 -i %d.bmp -c:v libx264rgb -r 24 output.mp4

# ToDo:
- [ ] Make canvas point scalar not a magic number
- [ ] Fix memory leaks
- [ ] Experiment with methods for shadows behind dieletrics
- [ ] Investigate Photon Mapping
- [ ] Investigate Normal Maps
- [ ] Investigate Environment maps

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
