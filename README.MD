# How to run
$ make run input=<path/to/xml>

Example:
$ make run input=inputs/empty_box_cd.xml



# Implementation Steps

## 1. Initialize the Project
- Use the provided helper files (`Scene.cpp` and `Helpers.cpp`) and study their structures and functions.
- Set up the `Makefile` to compile your program into the `rasterizer` executable.

## 2. Implement Modeling Transformations
- Apply transformations (translation, rotation, scaling) to all vertices of the triangles in the order specified in the input file.
- Use matrix multiplication for these transformations, as discussed in the "Modeling Transformations" lecture.

## 3. Implement Viewing Transformations
- Transform vertices from world coordinates to camera (view) coordinates.
- Implement the camera transformation matrix based on the camera position, gaze vector, and up vector from the input file.

## 4. Projection Transformations
- Implement both orthographic and perspective projections based on the camera type specified in the input file.

## 5. Backface Culling
- Calculate normal vectors for each triangle.
- Cull triangles whose backfaces face the viewer, based on the dot product between the normal and the view direction.

## 6. Clipping
- Choose a clipping algorithm (Cohen-Sutherland or Liang-Barsky) to clip triangles in **wireframe mode only**.
- Skip clipping in **solid mode**, assuming all meshes are fully visible.

## 7. Rasterization
- **Line Rasterization**: Use the midpoint algorithm for triangle edges in wireframe mode.
- **Triangle Rasterization**: Use barycentric coordinates for solid mode to fill the triangle surfaces.

## 8. Depth Buffering
- Maintain a z-buffer to track the depth of fragments.
- Ensure that only the closest triangle's fragments are drawn at each pixel.

## 9. Output Rendering
- Write the final image to a `.ppm` file in the resolution specified in the input file.





# Structure

-Implement forwardRenderingPipeline(camera) method of scene:
            // do forward rendering pipeline operations
            scene->forwardRenderingPipeline(scene->cameras[i]);



-Scene members:
	
    // Objects
    std::vector<std::vector<Color> > image;
    std::vector<Vec3 *> vertices;
    std::vector<Mesh *> meshes;

    // Coloring
    Color backgroundColor;
    std::vector<Color *> colorsOfVertices;
    int makeBetweenZeroAnd255(double value);

    // Culling
    bool cullingEnabled;
	
    // Modelling Transformations
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	
    
- forwardRenderingPipeline:
    - Take camera as input
    - Render each object according to scene members
    - Write on image according to given camera





# Inputs
empty_box_cd.xml                
flag_czechia_alternative_ce.xml 
flag_turkey_alternative_cd.xml  
horse_and_mug_cd.xml
empty_box_ce.xml                
flag_czechia_cd.xml             
flag_turkey_alternative_ce.xml  
horse_and_mug_ce.xml
empty_box_clipped.xml           
flag_czechia_ce.xml             
flag_turkey_cd.xml              
horse_and_mug_orthographic.xml
filled_box_cd.xml               
flag_germany_cd.xml             
flag_turkey_ce.xml              
horse_and_mug_perspective.xml
filled_box_ce.xml               
flag_germany_ce.xml             
flag_turkey_orthographic.xml    
sample_cd.xml
flag_czechia_alternative_cd.xml 
flag_iceland_ce.xml             
flag_turkey_perspective.xml     
sample_ce.xml