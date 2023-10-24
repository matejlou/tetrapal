# Tetrapal

![tetrapal_banner_e](https://github.com/matejlou/Tetrapal/assets/120740455/0b82367d-eb12-4c6f-aa5c-705c9171e51a)

Creates a tetrahedral tessellation from a colour palette via 3D Delaunay triangulation for use in arbitrary-palette colour image dithering.

# What is Tetrapal?

Tetrapal is shorthand for _tetrahedral palette_. It is a utility that takes a set of points in colour space and forms a [triangulated irregular network](https://en.wikipedia.org/wiki/Triangulated_irregular_network) of tetrahedra. As a result, any other point in colour space can be represented by the weighted sum of up to 4 existing palette colours.

The main motivation behind this library is to enable an efficient implementation of colour image [ordered dithering](https://en.wikipedia.org/wiki/Ordered_dithering) for irregular palettes. Typically, ordered dithering is only optimal when the palette contains colours that are equally distributed in colour space, i.e. that they form a regular grid. However, the algorithm can be modified to accommodate irregular or arbitrary palettes by representing input colours as the weighted sum of a number of existing palette colours. Tetrapal provides a data structure that can determine the necessary weights much faster and with more precision than existing implementations.

| Palette | Original Image | Typical Algorithm | Tetrapal Algorithm |
|-|-|-|-|
|![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/56246eaf-62ab-4484-82be-44b6c2df3700)|![dog256](https://github.com/matejlou/Tetrapal/assets/120740455/bd7d0c22-e946-49e4-bf01-2133687cc553)|![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/def812c7-a231-432e-9906-96edab438bbf)|![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/c1c69d07-069f-4558-862a-b19783f89401)|

The idea to use a 3D triangulated irregular network as a means to dither colour images is not a new one. The earliest source I could find that describes such an implementation is the 1988 article "_Using tetrahedrons for dithering color pictures_" by Eduard Gröller and Werner Purgathofer[^1]. However, the technique remains relatively unknown and unused in common practice, at least on the web. To my knowledge, this repository is the first and only public implementation available online.

# How to Use

To create a new tessellation from an existing palette simply call the following function:

```c
TetrapalData* tetrapal_new(const float *points, const int size);
```

The parameter `*points` should be a pointer to a buffer of 3D coordinates in the following format:

$$\Huge{[x_0, y_0, z_0, x_1, y_1, z_1, ... x_{size-1}, y_{size-1}, z_{size-1]}}$$

Where `size` is the number of 3D points represented in the buffer. Internally, points are indexed according to their order in the buffer, where the starting index is 0. If successful this function will return an opaque pointer to the Tetrapal data structure, otherwise it will return `NULL`. The values in `*points` should be normalised between 0.0 to 1.0; values beyond this range will be clamped.

Remember to free the triangulation data when you are done via:

```c
void tetrapal_free(TetrapalData* tetrapal);
```

---

To interpolate within a tessellation, pass the `TetrapalData` pointer to the function below:

```c
int tetrapal_interpolate(const TetrapalData* tetrapal, const float point[3], int* a, int* b, int* c, int* d, double* u, double* v, double* w, double* x);
```

This will return an `int` between 1 and 4 depending on the number of points contributing to the interpolant defined by `point[3]`. The indices of the points will be written to the values at `*a`, `*b`, `*c`, and `*d`, and their weights written to `*u`, `*v`, `*w`, `*x`, respectively. In the case where the number of points $N$ is less than 4, the unused indices/weights will be set to 0 and the first $N$ return values will contain the valid indices. For example, a return value of 2 means that the values at `*a`, `*b`, `*u`, and `*v` will contain the contributing indices and their weights, while the values at `*c`, `*d`, `*w`, and `*x` will be set to 0.

---

Tetrapal also supports nearest-neighbour queries within the tessellation. This can be faster than a standard linear search under the right circumstances. It is included for convenience:

```c
int tetrapal_nearest_neighbour(const TetrapalData* tetrapal, const float point[3]);
```

Returns the index of the nearest neighbour to the query point defined by `point[3]`.

---

Also included are a number of functions that can provide useful information about the tessellation itself:

```c
int tetrapal_number_of_dimensions(const TetrapalData* tetrapal);
int tetrapal_number_of_vertices(const TetrapalData* tetrapal);
int tetrapal_number_of_segments(const TetrapalData* tetrapal);
int tetrapal_number_of_triangles(const TetrapalData* tetrapal);
int tetrapal_number_of_tetrahedra(const TetrapalData* tetrapal);
```

Important to note is that segments shared across different triangles are only counted once, i.e. the edges [_a_, _b_] and [_b_, _a_] count as a single segment. Likewise, the triangles [_a_, _b_, _c_] and [_c_, _b_, _a_] count as one triangle. This information is cached immediately after triangulation so there is no processing overhead in calling these functions.

---

It is possible to write geometry data directly into buffers provided by the user, where each element is represented by a set of vertex indices (or in the case of `tetrapal_get_vertices`, a single index):

```c
int tetrapal_get_vertices(const TetrapalData* tetrapal, int* buffer);
int tetrapal_get_segments(const TetrapalData* tetrapal, int* buffer);
int tetrapal_get_triangles(const TetrapalData* tetrapal, int* buffer);
int tetrapal_get_tetrahedra(const TetrapalData* tetrapal, int* buffer);
```

For all functions, the size of `*buffer` must be _no less_ than the number of elements present in the triangulation multiplied by the number of indices representing each element. For example, to use `tetrapal_get_triangles` correctly the size of `*buffer` must be no smaller than `tetrapal_number_of_triangles` × 3. All functions will return non-zero on failure.

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
  TetrapalData* tetrapal = tetrapal_new(palette, palette_size);

  // Iterate over pixels in the input image
  for (int y = 0; y < image_height; y++)
  {
    for (int x = 0; x < image_width; x++)
    {
      // Get the current pixel from the input buffer
      const int image_index = x + y * image_width;
      const float* pixel = &input[image_index * 3];

      // Interpolate within the triangulation to get the candidates for the current pixel and their weights
      tetrapal_interpolate(tetrapal, pixel,
        &candidates[0], &candidates[1], &candidates[2], &candidates[3],
        &weights[0], &weights[1], &weights[2], &weights[3]);

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

At its core, Tetrapal is an implementation of 3D Delaunay triangulation using the [Bowyer-Watson incremental construction algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm). Tetrapal borrows ideas from computational geometry libraries such as [CGAL](https://www.cgal.org/) and [Geogram](https://github.com/BrunoLevy/geogram) to ensure accuracy and correctness. This includes the use of 'infinite' vertices to guarantee convexity, as well as the ability to gracefully handle degenerate inputs in the form of lower-dimensional triangulations (2D, 1D, 0D). A combination of extended precision integer arithmetic and static filtering via [interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic) is used to provide robustness for geometric predicates in cases where standard floating point arithmetic can fail to give an accurate result.

Tetrapal supports linear interpolation of points inside the convex hull by locating the enclosing tetrahedron and computing the 3D barycentric coordinates with respect to the query point (in lower-dimensional triangulations the barycentric coordinates are taken from the enclosing line segment or triangle). Points that lie outside the convex hull of the triangulation are projected onto the closest triangle/segment/vertex on the surface and the barycentric coordinates are taken from it instead.

Point location is performed via stochastic walk[^4], and [spatial indexing](https://en.wikipedia.org/wiki/Grid_(spatial_index)) is used to quickly locate a starting tetrahedron that is close enough to the query point.

As its main purpose is to process colour palettes, some assumptions have been made to simplify the algorithm. However, Tetrapal can still be used as a general-purpose Delaunay triangulation library provided its limitations are observed. These include:
* That it is unoptimised compared to more mature Delaunay triangulation libraries, especially for large numbers of points.
* That it expects the input to be normalised between 0.0 and 1.0.
* That the desired precision of the input does not exceed 1 / 16,777,215.

# Performance & Comparison

The three tables below compare the running time of a Tetrapal-based ordered dithering algorithm against the more well-known algorithms of Thomas Knoll[^2] and Joel Yliluoma[^3]. A different independent variable was chosen for each test (palette size, threshold matrix size, and input image size, respectively). The "_Yliluoma's ordered dithering algorithm 2_" variant of Yliluoma's algorithms was implemented for all tests. The construction of the Tetrapal data structure itself is included in the timings. All image/palette colours were transformed to linearised sRGB space prior to processing.

| Palette Size | Tetrapal   | Knoll   | Yliluoma | | Matrix Size | Tetrapal | Knoll   | Yliluoma  | | Image Size | Tetrapal | Knoll   | Yliluoma  |
| :--          | :-:        | :-:     | :-:      |-| :--         | :-:      | :-:     | :-:       |-| :--        | :-:      | :-:     | :-:       |
| 8            | 0.246s     | 0.830s  | 6.923s   | | 2x2         | 0.238s   | 0.100s  | 0.364s    | | 128x128    | 0.021s   | 0.092s  | 1.089s    |
| 16           | 0.246s     | 1.457s  | 14.837s  | | 4x4         | 0.211s   | 0.392s  | 2.882s    | | 256x256    | 0.049s   | 0.365s  | 4.226s    |
| 32           | 0.266s     | 2.733s  | 30.529s  | | 8x8         | 0.207s   | 1.494s  | 17.032s   | | 512x512    | 0.144s   | 1.455s  | 16.587s   |
| 64           | 0.318s     | 5.240s  | 61.521s  | | 16x16       | 0.205s   | 5.676s  | 92.412s   | | 1024x1024  | 0.523s   | 5.658s  | 65.638s   |
| 128          | 0.408s     | 10.194s | 126.047s | | 32x32       | 0.204s   | 22.421s | 470.400s  | | 2048x2048  | 2.039s   | 22.222s | 262.450s  |
| 256          | 0.467s     | 20.072s | 257.878s | | 64x64       | 0.200s   | 88.046s | 2297.365s | | 4096x4096  | 8.272s   | 87.802s | 1056.591s |

Tetrapal is faster in almost all cases and scales far better with both the palette size and threshold matrix size. This is because both Knoll and Yliluoma repeatedly iterate over the entire palette in order to fill an array of size $N$ with candidates, which is then sampled according to the value in the threshold matrix for the current pixel. For optimal results it is necessary that $N$ be equal to the number of elements in the threshold matrix, owing to the decrease in performance. Because Tetrapal determines candidates geometrically, the size of the threshold matrix does not matter. The palette size also has less of an impact on performance due to the nature of the triangulation and point location routines.

---

This table shows the [peak signal-to-noise ratio](https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio) (PSNR) for the output of each algorithm as a rough estimate of the dither quality (higher is better). A Gaussian blur was applied to the dithered output images before measuring the PSNR. Two different 16-colour palettes were tested for each image; a custom palette by Andrew Kensler[^5] that remained the same for all images, and an [adaptive palette](https://en.wikipedia.org/wiki/List_of_software_palettes#Adaptive_palettes)  generated using a variance-based colour quantisation algorithm by Xiaolin Wu[^6].

| Image              | Tetrapal | Knoll    | Yliluoma | Palette |
| :--                | :-:      | :-:      | :-:      |-|
| [Peppers](https://www.hlevkin.com/hlevkin/TestImages/pepper.bmp) (Custom)   | 27.17dB  | 27.03dB  | 27.03dB  | ![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/1140a0d4-f710-4e19-980e-2c0634bbc6f5)|
| [Peppers](https://www.hlevkin.com/hlevkin/TestImages/pepper.bmp) (Adaptive) | 24.33dB  | 24.37dB  | 24.37dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/ee06b360-2686-444b-83b2-03dd671e2e29)|
| [Baboon](https://www.hlevkin.com/hlevkin/TestImages/baboon.bmp) (Custom)    | 21.14dB  | 21.10dB  | 21.10dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/bbcf7aba-fa59-42c7-bfc7-c55179d8f51e)|
| [Baboon](https://www.hlevkin.com/hlevkin/TestImages/baboon.bmp) (Adaptive)  | 20.58dB  | 20.61dB  | 20.61dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/dc6883cd-565e-4617-b93e-46f5422e6dc1)|
| [Parrots](https://r0k.us/graphics/kodak/kodak/kodim23.png) (Custom)   | 27.77dB  | 27.65dB  | 27.65dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/0c17aae1-8d36-49fc-a418-efc701f9dc6e)|
| [Parrots](https://r0k.us/graphics/kodak/kodak/kodim23.png) (Adaptive) | 20.42dB  | 20.42dB  | 20.42dB  |![output_PALETTE](https://github.com/matejlou/Tetrapal/assets/120740455/a8549e82-9409-457e-b605-c134041b755f)|

Knoll and Yliluoma score exactly the same for each image, while the Tetrapal algorithm scores slightly better than the other two using the custom palette and slightly worse when using an adaptive palette. All algorithms scored worse using an adaptive palette, suggesting that a palette generation algorithm specifically catered towards dithering should be preferred in general. It is possible that Tetrapal is better suited to varied colour palettes with a good spread of colours, unlike the generated palettes which tended to contain many similar colours.

---

This next table records the size in memory of the Tetrapal data structure for various palette sizes, with random colours generated uniformly inside the unit cube. The memory consumption of the Knoll and Yliluoma algorithms that were implemented were trivial (equivalent to the size $N$ of the threshold matrix, which is usually <1KB).

| Palette Size | 8     | 16    | 32    | 64    | 128   | 256   |
| :--          | :-:   | :-:   | :-:   | :-:   | :-:   | :-:   | 
| **Memory**   | 2KB  | 4KB  | 8KB | 16KB | 32KB | 64KB |

---

Finally, here is a visual comparison between the dithered output of Tetrapal, Knoll, and a standard implementation of ordered dithering. Yliluoma has been omitted as the output is virtually identical to Knoll. The 16-colour [CGA](https://en.wikipedia.org/wiki/Color_Graphics_Adapter) palette was used.

| Algorithm | Test Image 1 | Test Image 2 |
| :-:       | :-:      | :-: |
| Tetrapal  |![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/f3770dee-6aab-4e01-865c-20c98172656d)|![output_DELAUNAY](https://github.com/matejlou/Tetrapal/assets/120740455/4958488d-ee14-434d-8a5b-b4c71b410638)|
| Knoll     |![output_KNOLL](https://github.com/matejlou/Tetrapal/assets/120740455/ec6dd924-d308-47e8-851a-23a4e8775bcf)|![output_KNOLL](https://github.com/matejlou/Tetrapal/assets/120740455/26de6284-d730-42cc-812a-f871872d2ef0)|
| Standard  |![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/280ea9c6-d8df-42df-b86d-88e83a472663)|![output_THRESHOLD](https://github.com/matejlou/Tetrapal/assets/120740455/5fc44c0f-18b5-4b6c-b3e0-2acfc5b9ba8d)|

[^1]: E. Gröller and W. Purgathofer, "_Using tetrahedrons for dithering color pictures_" (1988).

[^2]: T. Knoll, ["_Pattern Dithering_"](https://patents.google.com/patent/US6606166B1/en) (1999).

[^3]: J. Yliluoma, ["_Joel Yliluoma's arbitrary-palette positional dithering algorithm_"](https://bisqwit.iki.fi/story/howto/dither/jy/) (2011).

[^4]: O. Devillers, S. Pion and M. Teillaud, "[_Walking in a Triangulation_"](https://inria.hal.science/inria-00102194/document) (2006).

[^5]: A. Kensler, ["_Pixel Art Palettes for Free_"](http://eastfarthing.com/blog/2016-05-06-palette/) (2016).

[^6]: X. Wu, "_Efficient Statistical Computations For Optimal Colour Quantisation_" (1995).

[^7]: J. Shewchuk, ["_Adaptive Precision Floating-Point Arithmetic
and Fast Robust Geometric Predicates_"](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf) (1997).
