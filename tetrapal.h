#ifndef TETRAPAL_H
#define TETRAPAL_H

#ifdef __cplusplus
extern "C" {
#endif

	// =====================================================================
	//		Tetrapal Colour Palette Triangulator (Delaunay Triangulation)
	// =====================================================================

	/*
		The triangulation is generated using a variant of the Bowyer-Watson incremental insertion algorithm[1].
		Instead of an initial bounding 'super-tetrahedron', this algorithm makes use of an 'infinite' vertex
		that is connected to every face on the convex hull, which is typical of many Delaunay triangulation
		libraries. The infinite vertex simplifies insphere assessments for points that are inserted outside the
		convex hull of the current triangulation, and in turn ensures that the resulting triangulation is truly
		convex (it is common for a bounding tetrahedron consisting of finite vertices to produce concavities in
		the final triangulation). It also allows for convenient traversal along the convex hull itself.

		For the infinite vertex technique to work correctly, an initial starting tetrahedron must be created in
		advance. This poses a problem for point sets which are 1 or 2-dimensional, which, although uncommon, should
		still be considered valid input. To account for this, I took inspiration from CGAL's triangulation routines,
		and included support for 1D and 2D triangulations embedded in 3D. The algorithm will successively attempt to
		create the first d-simplex starting from d = 3, reducing the dimension by 1 for every failed attempt. Once
		a d-simplex is found, the point set will be triangulated using the appropriate d-triangulation algorithm.
		The Delaunay struct contains a variable specifying the dimensions of the point set, which is used to determine
		the appropriate functions to call after the initial triangulation is completed.

		To support adjecency queries within the triangulation, a hash table data structure is used, where a tetrahedron's
		vertex [d] is keyed by a tuple of incident vertices [a, b, c], such that [a, b, c, d] form a positively-oriented
		tetrahedron. This idea is borrowed from [2][pp.46], and I have found it to be extremely useful in providing a
		fast solution to adjecency queries while avoiding the added code and often headache that comes with manually
		keeping track of a tetrahedron's neighbours, such as in pointer-based methods.

		In addition, a graph is implemented in the form of an adjacency list in order to keep track of vertex connections
		along the convex hull. This data structure makes it possible to walk from face, to edge, to vertex, etc. and is
		in turn extremely helpful when determining the orthogonal projection of external points onto the hull.

		To determine the enclosing tetrahedra during point insertion, I utilise a 'stochastic walk' algorithm[4][pp.8].
		It consists of performing a 3D orientation test for every face of the current tetrahedron, in random order. If a
		test returns negative, then the point exists somewhere behind the plane, so we move into the adjacent tetrahedron
		and repeat the process. The procedure terminates when all tests return positive (point is inside) or when an
		infinite tetrahedron is reached.

		Spatial indexing is implemented to keep a record of one finite tetrahedron for every grid unit. Whenever a new point
		is inserted, the walking algorithm will start from an appropriate tetrahedron according to its location in the grid. 

		For numerical robustness in the computation of geometric tests, a static filter is used in combination with
		extended-precision exact arithmetic in the case of an uncertain result. First, the predicates are evaluated via
		interval arithmetic[3] using double precision floats with accumulating error bounds. If the signs of the upper
		and lower bounds are in agreement, then the sign of the result can be guaranteed and the algorithm can continue
		as normal. If the signs are not in agreement, then we resort to computing the exact value using the less performant
		extended arithmetic method. The filter tries to ensure that we only compute the more computationally intensive exact 
		value if it is necessary to do so.

		[1] https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
		[2] https://people.eecs.berkeley.edu/~jrs/meshpapers/delnotes.pdf
		[3] https://en.wikipedia.org/wiki/Interval_arithmetic
		[4] https://hal.inria.fr/inria-00102194/document
	*/

	/* Opaque handle to the Tetrapal data structure. */
	typedef struct TetrapalData TetrapalData;

	/* Triangulate a given point set, allocating memory for a new triangulation. */
	TetrapalData* tetrapal_new(const float *points, const int size);

	/* Free the memory representing an existing triangulation. */
	void tetrapal_free(TetrapalData* tetrapal);

	/* Interpolate a given input point as the weighted sum of up to four exisiting points in the triangulation.
		Returns an int from 1 - 4 signifying the number of points contributing to the interpolated value.
		Valid weights/indices are always given first. */
	int tetrapal_interpolate(const TetrapalData* tetrapal, const float point[3], int* a, int* b, int* c, int* d, double* u, double* v, double* w, double* x);

	/* Return the index of the nearest palette colour to the given input point. */
	int tetrapal_nearest_neighbour(const TetrapalData* tetrapal, const float point[3]);

	/* Returns the number of dimensions spanned by the triangulation. */
	int tetrapal_number_of_dimensions(const TetrapalData* tetrapal);

	/* Returns the number of vertices in the triangulation. */
	int tetrapal_number_of_vertices(const TetrapalData* tetrapal);

	/* Returns the number of segments in the triangulation. */
	int tetrapal_number_of_segments(const TetrapalData* tetrapal);

	/* Returns the number of triangles in the triangulation. */
	int tetrapal_number_of_triangles(const TetrapalData* tetrapal);

	/* Returns the number of tetrahedra in the triangulation. */
	int tetrapal_number_of_tetrahedra(const TetrapalData* tetrapal);

	/* Get the indices of every vertex in the triangulation.
		The buffer must be large enough to hold (tetrapal_number_of_vertices) indices. 
		A non-zero return value indicates that the function has failed to fill the buffer. */
	int tetrapal_get_vertices(const TetrapalData* tetrapal, int* buffer);

	/* Get the indices of every segment in the triangulation, where a single segment is represented by 2 indices at a time.
		The buffer must be large enough to hold (tetrapal_number_of_segments * 2) indices. 
		A non-zero return value indicates that the function has failed to fill the buffer. */
	int tetrapal_get_segments(const TetrapalData* tetrapal, int* buffer);

	/* Get the indices of every triangle in the triangulation, where a single triangle is represented by 3 indices at a time.
		The buffer must be large enough to hold (tetrapal_number_of_triangles * 3) indices. 
		A non-zero return value indicates that the function has failed to fill the buffer. */
	int tetrapal_get_triangles(const TetrapalData* tetrapal, int* buffer);

	/* Get the indices of every tetrahedron in the triangulation, where a single tetrahedron is represented by 4 indices at a time.
		The buffer must be large enough to hold (tetrapal_number_of_tetrahedra * 4) indices. 
		A non-zero return value indicates that the function has failed to fill the buffer. */
	int tetrapal_get_tetrahedra(const TetrapalData* tetrapal, int* buffer);

#ifdef __cplusplus
}
#endif

#endif // !TETRAPAL_H

/*
	MIT License

	Copyright (c) 2023 matejlou

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
*/
