# Tetrapal

![tetrapal_banner_f](https://github.com/matejlou/Tetrapal/assets/120740455/08e48e4b-1620-474e-9a7b-09b198ada4c2)

Creates a tetrahedral tessellation from a colour palette via 3D Delaunay triangulation. Intended for use in arbitrary-palette colour image dithering.

# What is Tetrapal?

Tetrapal is shorthand for _tetrahedral palette_. It is a utility that takes a set of points in colour space and forms a [triangulated irregular network](https://en.wikipedia.org/wiki/Triangulated_irregular_network) of tetrahedra. As a result, any colour can be represented as the weighted sum of up to 4 existing palette colours.

The main motivation behind this library is to enable an efficient implementation of colour image [ordered dithering](https://en.wikipedia.org/wiki/Ordered_dithering) for irregular palettes. Typically, ordered dithering is only optimal when the palette contains colours that are equally distributed in colour space, i.e. that they form a regular grid. However, the algorithm can be modified to accommodate irregular or arbitrary palettes by representing input colours as the weighted sum of a number of existing palette colours. Tetrapal provides a data structure that can determine the necessary weights much faster and with more precision than existing implementations.

| Typical Algorithm | Tetrapal Algorithm |
|-|-|
|![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/def812c7-a231-432e-9906-96edab438bbf)|![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/c1c69d07-069f-4558-862a-b19783f89401)|

The idea to use a 3D triangulated irregular network as a means to dither colour images is not a new one. The earliest source I could find that describes such an implementation is the 1988 article "_Using tetrahedrons for dithering color pictures_" by Eduard Gröller and Werner Purgathofer[^1]. However, the technique remains relatively unknown and unused in common practice, at least on the web. To my knowledge, this repository is the first and only public implementation available online.

# How to Use

## New Triangulation

```c
Tetrapal* tetrapal_new(const float *points, const int size);
```

Creates a new triangulation from an existing palette. The parameter `*points` should be a pointer to a buffer of 3D coordinates in the following format:

$$\Huge{[x_0, y_0, z_0, x_1, y_1, z_1, ... x_{size-1}, y_{size-1}, z_{size-1]}}$$

Where `size` is the number of points represented in the buffer. Internally, points are indexed according to the order they appear in the buffer, where the starting index is 0. If successful this function will return an opaque pointer to the Tetrapal data structure, otherwise it will return `NULL`.

Tetrapal expects coordinate values to be in the range 0.0 to 1.0; values beyond this range will be clamped. Internally, coordinates are transformed and represented as integers with 16 bits of precision.

## Free Triangulation

```c
void tetrapal_free(Tetrapal* tetrapal);
```

Triangulation data can be safely freed by passing the `Tetrapal` pointer to the above function.

## Interpolation
### Barycentric Interpolation

```c
int tetrapal_interpolate(const Tetrapal* tetrapal, const float point[3], int* indices, float* weights);
```

Performs barycentric interpolation within a triangulation. This will return an `int` between 1 and 4 depending on the number of points contributing to the interpolant given by `point[3]`. The indices of the points will be written to the values at `*indices`, and their respective weights written to `*weights`. In the case where the number of points $N$ is less than 4, only the first $N$ values in each output array will be written. 

Naturally, the output arrays should be large enough to hold the maximum expected number of points, which is equivalent to `tetrapal_element_size`.

This is the recommended interpolation function as it provides the best combination of quality and performance.

### Natural Neighbour Interpolation

```c
int tetrapal_natural_neighbour(const Tetrapal* tetrapal, const float point[3], int* indices, float* weights, const int size);
```

In addition to standard barycentric interplation, Tetrapal is able to perform [natural neighbour interpolation](https://en.wikipedia.org/wiki/Natural_neighbor_interpolation) as well. This function returns an `int` specifiying the number of natural neighbours contributing to the interpolant given by `point[3]`. The size of the output arrays `*indices` and `*weights` should be given by the parameter `size`. It is possible for this function to fail, either due to a lack of available system memory or because the number of natural neighbours exceeds the size of the output arrays. In both cases, the function will return 0. 

In theory number of natural neighbours for a given input point can range from 1 to the total number of vertices in the triangulation. In practice the maximum number is much lower for most point sets. However, to guarantee success it is advised that the output arrays are at least as large as the number of vertices in the triangulation.

Because the function may fail, it is recommended to use barycentric interpolation for most cases. Natural neighbour interpolation is much slower and the resulting dither quality may not be much better than barycentric interpolation, perceptually speaking.

### Nearest Neighbour Interpolation

```c
int tetrapal_nearest_neighbour(const Tetrapal* tetrapal, const float point[3]);
```

Tetrapal also supports nearest-neighbour queries within the triangulation. This can be faster than a standard linear search under the right circumstances. It is included for convenience.

Returns the index of the nearest neighbour to the query point defined by `point[3]`.

## Utility Functions
### Number of Elements

```c
int tetrapal_number_of_elements(const Tetrapal* tetrapal);
```

Returns the number of real elements (simplices) in the triangulation, i.e. does not count symbolic 'infinite' elements.

### Number of Dimensions

```c
int tetrapal_number_of_dimensions(const Tetrapal* tetrapal);
```

Returns the number of dimensions spanned by the triangulation. This can range anywhere from 0 to 3. A value of -1 indicates a null triangulation.

### Element Size

```c
int tetrapal_element_size(const Tetrapal* tetrapal);
```

Gets the size of an element in the triangulation, or in other words, the number of vertices that make up a simplex within the triangulation. This is the same as `tetrapal_number_of_dimensions` + 1.

### Get Elements

```c
int tetrapal_get_elements(const Tetrapal* tetrapal, int* buffer);
```

Writes raw geometry data directly into `*buffer`. Every element is represented by a set of vertex indices which are packed contiguously into the buffer. Thus, the size of the buffer must be _no less_ than the number of elements present in the triangulation multiplied by the number of indices representing each element, or `tetrapal_number_of_elements` × `tetrapal_element_size`.

# Example Usage

Below is an example C implementation of a function that takes an input image and uses Tetrapal to perform ordered dithering with an arbitrary palette, producing an indexed image. The buffers `*input` and `*palette` are expected to contain RGB values from 0.0 to 1.0. The function `sort_by_luminance` simply orders the candidates and their corresponding weights in order of ascending/descending [luminance](https://en.wikipedia.org/wiki/Luma_(video)).

```c

// Normalised threshold matrix
const double bayer_matrix_4x4[4][4] =
{
  { 0.0  / 16.0,  12.0 / 16.0,  3.0  / 16.0,  15.0 / 16.0 },
  { 8.0  / 16.0,  4.0  / 16.0,  11.0 / 16.0,  7.0  / 16.0 },
  { 2.0  / 16.0,  14.0 / 16.0,  1.0  / 16.0,  13.0 / 16.0 },
  { 10.0 / 16.0,  6.0  / 16.0,  9.0  / 16.0,  5.0  / 16.0 }
};

void dither_image(const float *input, const int image_width, const int image_height,
  unsigned char *output, const float *palette, const int palette_size)
{
  // Candidate arrays
  int candidates[4];
  double weights[4];

  // Triangulate the palette
  Tetrapal* tetrapal = tetrapal_new(palette, palette_size);

  // Iterate over pixels in the input image
  for (int y = 0; y < image_height; y++)
  {
    for (int x = 0; x < image_width; x++)
    {
      // Get the current pixel from the input buffer
      const int image_index = x + y * image_width;
      const float* pixel = &input[image_index * 3];

      // Interpolate within the triangulation to get the candidates for the current pixel and their weights
      tetrapal_interpolate(tetrapal, pixel, candidates, weights);

      // Sort the candidates by luminance
      sort_by_luminance(candidates, weights, palette);
     
      // Use the value in the threshold matrix to select a candidate
      const double threshold = bayer_matrix_4x4[y % 4][x % 4];
      double sum = 0.0;

      // Accumulate the sum of weights until we pass the threshold
      for (int i = 0; i < 4; i++)
      {
        sum += weights[i];
       
        if (threshold < sum)
        {
          output[image_index] = candidates[i];
          break;
        }
      }
    }
  }

  // Remember to free the triangulation data
  tetrapal_free(tetrapal);
}
```

# Implementation

Tetrapal implements 3D Delaunay triangulation using the [Bowyer-Watson incremental construction algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm). It borrows ideas from computational geometry libraries such as [CGAL](https://www.cgal.org/) and [Geogram](https://github.com/BrunoLevy/geogram) to ensure accuracy and correctness. This includes the use of symbolic 'infinite' vertices to guarantee convexity, as well as the ability to gracefully handle degenerate inputs in the form of lower-dimensional triangulations (2D, 1D, 0D). Coincident points are ignored during triangulation but are still stored internally to maintain consistency when indexing vertices. 

A combination of extended precision integer arithmetic and static filtering via predetermined error bounds[^2] is used to provide robustness for geometric predicates in cases where standard arithmetic may fail to give an accurate result.

Tetrapal supports extrapolation of points outside the convex hull via projection onto the nearest surface element (which may be a triangle, line segment, or single vertex) and returning the barycentric coordinates with respect to this element. Natural-neighbour interpolation with Sibson weights is supported in 2D and 3D, but only for points that lie within the convex hull.

Point location is performed via remembering stochastic walk[^3], where a kd-tree approximate nearest neighbour search is used to quickly locate a starting tetrahedron close to the query point.

As its main purpose is to process colour palettes, some assumptions have been made to simplify the algorithm. However, Tetrapal can still be used as a general-purpose Delaunay triangulation library provided its limitations are observed. These include:
* That it expects the input to be normalised between 0.0 and 1.0.
* That the desired precision of the input does not exceed 1 / 65535.

# Performance & Comparison

## Time

The three tables below compare the running time of a Tetrapal-based ordered dithering algorithm against the more well-known algorithms of Thomas Knoll[^4] and Joel Yliluoma[^5]. A different independent variable was chosen for each test (palette size, threshold matrix size, and input image size, respectively). The "_Yliluoma's ordered dithering algorithm 2_" variant of Yliluoma's algorithms was implemented. The construction of the Tetrapal data structure itself is included in the timings, and barycentric interpolaton was used for all tests. All image/palette colours were transformed to linearised sRGB space prior to processing.

| Palette Size | Tetrapal   | Knoll   | Yliluoma | | Matrix Size | Tetrapal | Knoll   | Yliluoma  | | Image Size | Tetrapal | Knoll   | Yliluoma  |
| :--          | :-:        | :-:     | :-:      |-| :--         | :-:      | :-:     | :-:       |-| :--        | :-:      | :-:     | :-:       |
| 8            | 0.091s     | 0.830s  | 6.923s   | | 2x2         | 0.127s   | 0.100s  | 0.364s    | | 128x128    | 0.010s   | 0.092s  | 1.089s    |
| 16           | 0.133s     | 1.457s  | 14.837s  | | 4x4         | 0.118s   | 0.392s  | 2.882s    | | 256x256    | 0.038s   | 0.365s  | 4.226s    |
| 32           | 0.146s     | 2.733s  | 30.529s  | | 8x8         | 0.122s   | 1.494s  | 17.032s   | | 512x512    | 0.146s   | 1.455s  | 16.587s   |
| 64           | 0.166s     | 5.240s  | 61.521s  | | 16x16       | 0.128s   | 5.676s  | 92.412s   | | 1024x1024  | 0.523s   | 5.658s  | 65.638s   |
| 128          | 0.178s     | 10.194s | 126.047s | | 32x32       | 0.124s   | 22.421s | 470.400s  | | 2048x2048  | 2.039s   | 22.222s | 262.450s  |
| 256          | 0.186s     | 20.072s | 257.878s | | 64x64       | 0.129s   | 88.046s | 2297.365s | | 4096x4096  | 8.272s   | 87.802s | 1056.591s |

Tetrapal is faster in almost all cases. This is because both Knoll and Yliluoma are iterative algorithms whose time complexity is a factor of both the palette size $n$ as well as the number of candidates $m$ per pixel, which for optimal results is typically proportional to the size of the threshold matrix. Tetrapal's time complexity is bounded by its point location routine, which is shown to have an expected running time of $O(n^{1/4})$[^6]. The table below shows the theoretical worst-case time complexity for each algorithm.

| Algorithm      | Tetrapal       | Knoll         | Yliluoma     |
| :--            | :-:            | :-:           | :-:          |
| **Complexity** | $O(n^{1/4})$   | $O(nm)$       | $O(nmlogn)$  |

## Quality

This table shows the [peak signal-to-noise ratio](https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio) (PSNR) for the output of each algorithm as a rough estimate of the dither quality (higher is better). A Gaussian blur was applied to the dithered output images before measuring the PSNR. Two different 16-colour palettes were tested for each image; a custom palette by Andrew Kensler[^7] that remained the same for all images, and an [adaptive palette](https://en.wikipedia.org/wiki/List_of_software_palettes#Adaptive_palettes)  generated using a variance-based colour quantisation algorithm by Xiaolin Wu[^8].

| Image              | Tetrapal | Knoll    | Yliluoma | Palette |
| :--                | :-:      | :-:      | :-:      |-|
| [Peppers](https://www.hlevkin.com/hlevkin/TestImages/pepper.bmp) (Custom)   | 27.17dB  | 27.03dB  | 27.03dB  | ![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/1140a0d4-f710-4e19-980e-2c0634bbc6f5)|
| [Peppers](https://www.hlevkin.com/hlevkin/TestImages/pepper.bmp) (Adaptive) | 24.33dB  | 24.37dB  | 24.37dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/ee06b360-2686-444b-83b2-03dd671e2e29)|
| [Baboon](https://www.hlevkin.com/hlevkin/TestImages/baboon.bmp) (Custom)    | 21.14dB  | 21.10dB  | 21.10dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/bbcf7aba-fa59-42c7-bfc7-c55179d8f51e)|
| [Baboon](https://www.hlevkin.com/hlevkin/TestImages/baboon.bmp) (Adaptive)  | 20.58dB  | 20.61dB  | 20.61dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/dc6883cd-565e-4617-b93e-46f5422e6dc1)|
| [Parrots](https://r0k.us/graphics/kodak/kodak/kodim23.png) (Custom)   | 27.77dB  | 27.65dB  | 27.65dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/0c17aae1-8d36-49fc-a418-efc701f9dc6e)|
| [Parrots](https://r0k.us/graphics/kodak/kodak/kodim23.png) (Adaptive) | 20.42dB  | 20.42dB  | 20.42dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/a8549e82-9409-457e-b605-c134041b755f)|

Knoll and Yliluoma score exactly the same for each image, while the Tetrapal algorithm scores slightly better than the other two using the custom palette and slightly worse when using an adaptive palette. All algorithms scored worse using an adaptive palette, suggesting that a palette generation algorithm specifically catered towards dithering should be preferred in general. It is possible that Tetrapal is better suited to varied colour palettes with a good spread of colours, unlike the generated palettes which tended to contain many similar colours.

## Memory

This table records the size in memory of the Tetrapal data structure for various palette sizes, with random colours generated uniformly inside the unit cube. The memory consumption of the Knoll and Yliluoma algorithms that were implemented were trivial (equivalent to the size $N$ of the threshold matrix in the worst case, which is usually <1KB).

| Palette Size | 8     | 16    | 32    | 64    | 128   | 256   |
| :--          | :-:   | :-:   | :-:   | :-:   | :-:   | :-:   | 
| **Memory**   | 2KB   | 4KB   | 8KB   | 16KB  | 32KB  | 64KB  |

## Visuals

Here is a visual comparison between the dithered output of Tetrapal, Knoll, and a standard implementation of ordered dithering. Yliluoma has been omitted as the output is virtually identical to Knoll. The 16-colour [CGA](https://en.wikipedia.org/wiki/Color_Graphics_Adapter) palette was used.

| Algorithm     | Test Image 1 | Test Image 2 |
| :-:           | :-:      | :-: |
| **Tetrapal**  |![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/f3770dee-6aab-4e01-865c-20c98172656d)|![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/4958488d-ee14-434d-8a5b-b4c71b410638)|
| **Knoll**     |![output_KNOLL](https://github.com/matejlou/Tetrapal/assets/120740455/ec6dd924-d308-47e8-851a-23a4e8775bcf)|![output_KNOLL](https://github.com/matejlou/Tetrapal/assets/120740455/26de6284-d730-42cc-812a-f871872d2ef0)|
| **Standard**  |![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/280ea9c6-d8df-42df-b86d-88e83a472663)|![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/5fc44c0f-18b5-4b6c-b3e0-2acfc5b9ba8d)|

## Interpolation

These final images are intended to show the visual difference between dithering via barycentric interpolation and natural neighbour interpolation. A 256x256 void-and-cluster threshold matrix was used for both images. In general, natural neighbour dithering produces perceptually smoother gradations but introduces more high-frequency noise as a result of considering a greater number of candidate colours per pixel. This can best be seen in the appearance of the sky in each of the example images.

| Barycentric | Natural Neighbour |
| :-:         | :-:               |
|![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/d85c74e4-db01-4188-9c94-59827db79399)|![output_NATURAL_NEIGHBOUR](https://github.com/matejlou/Tetrapal/assets/120740455/8d3fc6a3-4556-47c5-a7ab-8d10c22d972a)|

[^1]: E. Gröller and W. Purgathofer, "_Using tetrahedrons for dithering color pictures_" (1988).

[^2]: S. Fortune and C. J. Van Wyk, "_Static Analysis Yields Efficient Exact Integer Arithmetic for Computational Geometry_" (1996).

[^3]: O. Devillers, S. Pion and M. Teillaud, "[_Walking in a Triangulation_"](https://inria.hal.science/inria-00102194/document) (2006).

[^4]: T. Knoll, ["_Pattern Dithering_"](https://patents.google.com/patent/US6606166B1/en) (1999).

[^5]: J. Yliluoma, ["_Joel Yliluoma's arbitrary-palette positional dithering algorithm_"](https://bisqwit.iki.fi/story/howto/dither/jy/) (2011).

[^6]: E. P. Mücke, I. Saias and B. Zhu, "_Fast randomized point location without preprocessing in two- and
three-dimensional Delaunay triangulations_" (1998).

[^7]: A. Kensler, ["_Pixel Art Palettes for Free_"](http://eastfarthing.com/blog/2016-05-06-palette/) (2016).

[^8]: X. Wu, "_Efficient Statistical Computations For Optimal Colour Quantisation_" (1995).

