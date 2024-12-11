# Software Rasterizer: Forward Rendering Pipeline Implementation

This project is an implementation of a software rasterizer in C++, focusing on the forward rendering pipeline. It processes 3D scenes composed of triangles to produce 2D images, supporting various rendering and transformation functionalities.

## Features

- **Modeling and Viewing Transformations**:
  - Translation, Rotation, and Scaling transformations.
  - Perspective and Orthographic projections.
- **Rasterization**:
  - Line drawing using the midpoint algorithm.
  - Triangle rasterization using barycentric coordinates.
- **Rendering Modes**:
  - Wireframe and Solid modes.
  - Backface culling (enabled/disabled via input).
- **Clipping**:
  - Implements Cohen-Sutherland or Liang-Barsky algorithm for wireframe mode.
- **Depth Buffering**:
  - Ensures proper rendering order of overlapping triangles.
- **Multiple Cameras**:
  - Supports multiple camera configurations with adjustable parameters.

## Input File Format

The input file is an XML document containing:

- **Scene Settings**:
  - Background color, backface culling settings, and camera configurations.
- **Geometric Data**:
  - Vertex positions and colors.
  - Transformation definitions (translations, rotations, scalings).
  - Mesh definitions and triangle connectivity.

For a detailed structure, see the provided `Description.pdf`.

## Output

The rasterizer generates 2D images in the `.ppm` format. Use tools like [ImageMagick](https://imagemagick.org) to convert `.ppm` files to viewable formats (e.g., PNG).

## Compilation and Execution

To build and run the project:

1. Compile the code using the provided `Makefile`:
   ```bash
   make rasterizer


   Here’s the full updated README.md file with the Makefile commands section included:

# Software Rasterizer: Forward Rendering Pipeline Implementation

This project is an implementation of a software rasterizer in C++, focusing on the forward rendering pipeline. It processes 3D scenes composed of triangles to produce 2D images, supporting various rendering and transformation functionalities.

## Features

- **Modeling and Viewing Transformations**:
  - Translation, Rotation, and Scaling transformations.
  - Perspective and Orthographic projections.
- **Rasterization**:
  - Line drawing using the midpoint algorithm.
  - Triangle rasterization using barycentric coordinates.
- **Rendering Modes**:
  - Wireframe and Solid modes.
  - Backface culling (enabled/disabled via input).
- **Clipping**:
  - Implements Cohen-Sutherland or Liang-Barsky algorithm for wireframe mode.
- **Depth Buffering**:
  - Ensures proper rendering order of overlapping triangles.
- **Multiple Cameras**:
  - Supports multiple camera configurations with adjustable parameters.

## Input File Format

The input file is an XML document containing:

- **Scene Settings**:
  - Background color, backface culling settings, and camera configurations.
- **Geometric Data**:
  - Vertex positions and colors.
  - Transformation definitions (translations, rotations, scalings).
  - Mesh definitions and triangle connectivity.

For a detailed structure, see the provided `Description.pdf`.

## Output

The rasterizer generates 2D images in the `.ppm` format. Use tools like [ImageMagick](https://imagemagick.org) to convert `.ppm` files to viewable formats (e.g., PNG).
