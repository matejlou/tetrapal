
#define TETRAPAL_IMPLEMENTATION
#include "tetrapal.h"

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>

#define TETRAPAL_MALLOC(size) malloc(size)
#define TETRAPAL_REALLOC(ptr, size) realloc(ptr, size)
#define TETRAPAL_FREE(ptr) free(ptr)

//#define TETRAPAL_DEBUG
#ifdef TETRAPAL_DEBUG
	#include <time.h>
	#include <stdio.h>
	clock_t TIME_LOCATE;
	clock_t TIME_STELLATE;
#endif

// =========================================================
//		Geometric Predicates
// =========================================================

/* Returns the squared distance between [a] and [b]. */
static double distance_squared(const int64_t a[3], const int64_t b[3]);

/* Returns 1 if [a] and [b] are coincident, returns 0 otherwise. */
static int is_coincident(const int64_t a[3], const int64_t b[3]);

/* Determine whether the projection of point [c] on the line [a, b] lies before or after point [a].
	Returns 0 if the projection of [c] is coincident with [a]. */
static double orient1d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3]);
static int orient1d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3]);
static int orient1d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3]);

/* Determine whether the projection of point [c] on the plane defined by [d, e, f] lies on the positive side
	of the line [a, b]. Assumes [a, b] and [d, e, f] are coplanar. */
static double orient2d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3]);
static int orient2d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3]);
static int orient2d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3]);

/* Determine whether the points [a, b, c, d] form a positively-oriented tetrahedron such that the triangle
	[a, b, c] has a counterclockwise winding order when seen from [d]. */
static double orient3d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);
static int orient3d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);
static int orient3d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);

/* Determine whether the point [c] lies on the line segment defined by the points [a, b].
	Assumes [a, b, c] are colinear. Returns positive if it is on the segment, negative if outside,
	and 0 if it lies exactly on an end point. */
// static double insegment_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3]);
static int insegment_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3]);
static int insegment_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3]);

/* Determine whether the point [d] lies within the circle passing through [a, b, c]. All
	points are given in 3D space, and [d] is assumed to be coplanar with [a, b, c]. There
	is no guarantee of correct behaviour if [d] is not coplanar! */
// static double incircle_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);
static int incircle_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);
static int incircle_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3]);

/* Determine whether the point [e] lies in, on, or outside the sphere passing through [a, b, c, d]. */
// static double insphere_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3]);
static int insphere_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3]);
static int insphere_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3]);

// =========================================================
//		Interval Arithmetic
// =========================================================

/*
	Interval representation for floating point values. Used as a static filter
	for geometric predicates.
*/

typedef struct
{
	double lower, upper;
} Interval;

static Interval interval(double value);
static Interval interval_add(const Interval lhs, const Interval rhs);
static Interval interval_sub(const Interval lhs, const Interval rhs);
static Interval interval_mul(const Interval lhs, const Interval rhs);
static int interval_sign(const Interval interval);
static inline void interval_swap_bounds(Interval* interval);
static inline double interval_min(const double a, const double b);
static inline double interval_max(const double a, const double b);
static inline void interval_add_error(Interval* interval);

#ifdef TETRAPAL_DEBUG
static void interval_print(const Interval interval);
#endif

// =========================================================
//		Extended Precision Integer Arithmetic (BigInt)
// =========================================================

/*
	Struct to represent large signed integers. Digits are represented in base 2^32,
	stored using 64-bit unsigned integer types to handle the carry from arithmetic operations.
	No dynamic allocation, so not technically 'arbitrary precision', more like 'extended
	precision'.

	For now, every digit in the array gets evaluated regardless of the length of the actual number.
	Improving this would be a welcome change, but it's not high priority yet.
*/

enum { BIGINT_SIZE = 8 };

typedef struct
{
	uint64_t digits[BIGINT_SIZE];
	int sign;
} BigInt;

static BigInt bigint(const int64_t val);
static BigInt bigint_zero();
static BigInt bigint_abs(const BigInt bigint);
static BigInt bigint_neg(const BigInt bigint);
static BigInt bigint_add(const BigInt lhs, const BigInt rhs);
static BigInt bigint_sub(const BigInt lhs, const BigInt rhs);
static BigInt bigint_mul(const BigInt lhs, const BigInt rhs);
static int bigint_compare(const BigInt lhs, const BigInt rhs);
static void bigint_shift_l(BigInt* bigint, unsigned int shift);

#ifdef TETRAPAL_DEBUG
static int bigint_is_zero(const BigInt bigint);
static void bigint_print_digits(const BigInt bigint);
static void bigint_print_decimal(const BigInt bigint);
#endif

// =========================================================
//		Barycentric Interpolation
// =========================================================

/*  Get the barycentric coordinates [u, v] for the projection of point [q] onto the line segment [a, b]. */
static void barycentric_1d(const int64_t q[3], const int64_t a[3], const int64_t b[3], double* u, double* v);

/* Get the barycentric coordinates [u, v, w] for the projection of point [q] onto the positively oriented triangle [a, b, c]. */
static void barycentric_2d(const int64_t q[3], const int64_t a[3], const int64_t b[3], const int64_t c[3], double* u, double* v, double* w);

/* Get the barycentric coordinates [u, v, w, x] for point [q] wrt the positively oriented tetrahedron [a, b, c, d]. */
static void barycentric_3d(const int64_t q[3], const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], double* u, double* v, double* w, double* x);

// =========================================================
//		Spatial Grid
// =========================================================

/* Resolution of the spatial grid. */
enum { SPATIALGRID_N = 16 }; 

typedef struct
{
	int grid[SPATIALGRID_N][SPATIALGRID_N][SPATIALGRID_N][4];
} SpatialGrid;

/* Build the spatial grid. */
static void spatialgrid_build(TetrapalData* tetrapal);

/* Locate the simplex associated with the grid index for a given point. */
static void spatialgrid_locate(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, int* d);

// =========================================================
//		Vertex Graph (Adjacency List)
// =========================================================

static const size_t GRAPH_LIST_RESERVE = 4;

typedef struct
{
	size_t capacity, size;
	int* data;
} GraphList;

typedef struct
{
	size_t size;
	GraphList* list;
} Graph;

static Graph* graph_new(size_t size);
static void graph_free(Graph* graph);

static void graph_list_insert(GraphList* list, int a);
// static void graph_list_remove(GraphList* list, int a);

/* Connect [a] and [b] in the graph. */
static void graph_insert(Graph* graph, int a, int b);

// /* Disconnect [a] and [b] in the graph. */
// static void graph_remove(Graph* graph, int a, int b);

/* Return the size of a given vertex list. */
static size_t graph_list_size(const Graph* graph, const size_t vertex);

/* Return a connected vertex at the index in the given list. */
static int graph_get(const Graph* graph, const size_t vertex, const size_t index);

/* Iterate through the graph and return the number of unique (undirected) edges. */
static int graph_find_number_of_unique_edges(const Graph* graph);

// ===============================================================
//	Hash Map 
// ===============================================================	

	/*
		An unordered map structure where values are associated with keys via a hash table. This implementation uses separate
		chaining via contiguous arrays rather than a linked list. Additionally, keys and values are stored a separate
		contiguous array for fast iteration. The hash table simply holds an index into these arrays for the appropriate
		element.

		The map is ordered until a removal is called, in which case the element to be removed is swapped out with the element
		at the back of the array.
	*/

typedef size_t(*HashMap_Func_Hash)(const void*);
typedef int(*HashMap_Func_Equal)(const void*, const void*);

typedef struct
{
	size_t size, capacity;
	size_t data;
} HashMap_Bucket;

typedef struct
{
	size_t capacity;
	HashMap_Bucket** buckets;
} HashMap_Table;

typedef struct
{
	size_t size, capacity;
	size_t stride_key, stride_value;
	void* data;
	HashMap_Table table;
	HashMap_Func_Hash func_hash;
	HashMap_Func_Equal func_equal;
} HashMap;

static const double HASHMAP_LOAD_FACTOR = 0.7;
static const double HASHMAP_GROWTH_FACTOR = 1.618;
static const size_t HASHMAP_BUCKET_RESERVE = 2;

static HashMap* hashmap_new(const size_t stride_key, const size_t stride_value, size_t reserve, const HashMap_Func_Hash func_hash, const HashMap_Func_Equal func_equal);

/* Free all data in the map. */
static void hashmap_free(HashMap* hashmap);

/* Insert an element into the map. */
static void hashmap_insert(HashMap* hashmap, const void* key, const void* value);

/* Remove an element from the map. Does not maintain the order of elements. */
static void hashmap_remove(HashMap* hashmap, const void* key);

// /* Check if the element exists in the map. Returns 0 if false, non-zero if true. */
// static int hashmap_has(const HashMap* hashmap, const void* key);

/* Lookup a value in the map via its key. */
static void* hashmap_find(const HashMap* hashmap, const void* key);

/* Get the key at the specified index. Elements are not guaranteed to be ordered! */
static void* hashmap_get_key(const HashMap* hashmap, const size_t index);

/* Get the value at the specified index. Elements are not guaranteed to be ordered! */
static void* hashmap_get_value(const HashMap* hashmap, const size_t index);

/* Return the number of elements in the map. */
static size_t hashmap_size(const HashMap* hashmap);

/* Increase the capacity of the hash table, rehashing all elements */
static void hashmap_table_grow(HashMap* hashmap);

// ===============================================================
//	Hash Set 
// ===============================================================	

/*
	An unordered set data structure that stores elements in a contiguous array. A hash table is used to map values
	to their position in the array. This data structure supports O(1) insertion and lookup, as well as fast iteration
	of all elements in the flat array. O(1) removal is also supported, as the structure simply replaces the element
	to be removed with the last element before decreasing the size. The hash table is updated in a similar manner.
	Because of this, the elements are only ordered as long as no removals are performed.
*/

typedef size_t(*HashSet_Func_Hash)(const void*);
typedef int(*HashSet_Func_Equal)(const void*, const void*);

typedef struct HashSet_Bucket {
	size_t size, capacity;
	size_t data;
} HashSet_Bucket;

typedef struct HashSet_Map {
	size_t capacity;
	HashSet_Bucket** buckets;
} HashSet_Map;

typedef struct HashSet {
	size_t size, capacity, stride;
	void* data;
	HashSet_Map map;
	HashSet_Func_Hash func_hash;
	HashSet_Func_Equal func_equal;
} HashSet;

static const double HASHSET_LOAD_FACTOR = 0.7;
static const double HASHSET_GROWTH_FACTOR = 1.618;
static const size_t HASHSET_BUCKET_RESERVE = 2;

/* Allocate memory for a new set. */
static HashSet* hashset_new(const size_t stride, size_t reserve, const HashSet_Func_Hash func_hash, const HashSet_Func_Equal func_equal);

/* Free all data in the set. */
static void hashset_free(HashSet* hashset);

/* Insert an element into the set. */
static void hashset_insert(HashSet* hashset, const void* element);

// /* Remove an element from the set. Does not maintain the order of elements. */
// static void hashset_remove(HashSet* hashset, const void* element);

// /* Check if the element exists in the set. Returns 0 if false, non-zero if true. */
// static int hashset_has(const HashSet* hashset, const void* element); 

/* Get the element at the specified index. Elements are not guaranteed to be ordered! */
static const void* hashset_get(const HashSet* hashset, const size_t index);

/* Return a pointer to the first element in the set. */
static const void* hashset_begin(const HashSet* hashset);

/* Return a pointer to the last element in the set. */
// static const void* hashset_end(const HashSet* hashset);

/* Return the number of elements in the set. */
static size_t hashset_size(const HashSet* hashset);

/* Increase the capacity of the map, rehashing all elements */
static void hashset_map_grow(HashSet* hashset);

// ===============================================================
//	Tetrapal Core
// ===============================================================

typedef enum
{
	VERTEX_INFINITE = -1,
	VERTEX_NULL = -2,
	VERTEX_NEXT = -3,
	VERTEX_PREV = -4
} TetrapalEnum;

typedef struct
{
	int64_t x, y, z;
} TetrapalVertex;

typedef struct
{
	int a, b, c;
} TetrapalAdjacencyKey;

struct TetrapalData
{
	size_t size;
	size_t dimensions;
	int last[4];
	int number_of_segments;
	int number_of_triangles;
	int number_of_tetrahedra;
	int64_t orient_a[3], orient_b[3], orient_c[3];
	TetrapalVertex* vertices;
	HashMap* adjacency;
	Graph* graph;
	SpatialGrid grid;
};

/* Internal precision for input data.
	Input float values from 0.0 to 1.0 are scaled up by this value and then cast to int64. */
static const unsigned int TETRAPAL_PRECISION = 16777216;
static const uint32_t TETRAPAL_RANDOM_MAX = 0xffff;

/* Generate a random integer. */
static uint32_t random_int(uint32_t* state);

/* Generate a random integer from 0 to (max - 1). */
static uint32_t random_range_int(uint32_t* state, const int32_t range);

/* Swap the values of two ints, given their pointers. */
static void swap_int(int* a, int* b);

/* Clamp a given float within the range spcified by min, max. */
static void clamp_float(float* value, float min, float max);

/* Given a vertex at index i, get its 3D coordinates.
If the vertex is infinite or doesn't exist, the function will return {-1, -1, -1}. */
static void get_coordinates(const TetrapalData* tetrapal, const int i, int64_t coords[3]);

/* Hash function for the adjacency map. */
static size_t adjacency_hash(const void* ptr_key);

/* Comparison function for adjacency keys. */
static int adjacency_compare_equal(const void* a, const void* b);

/* Hash function for a line segment defined by two indices. */
static size_t segment_hash(const void* ptr_segment);

/* Comparison function for a non-directed line segment defined by two indices. */
static int segment_compare_equal(const void* ptr_a, const void* ptr_b);

/* Arranges triangle indices to be in canonical order. */
static void triangle_make_canonical(int* a, int* b, int* c);

/* Hash function for a triangle defined by three indices. */
static size_t triangle_hash(const void* ptr_triangle);

/* Comparison function for a triangle defined by three indices. */
static int triangle_compare_equal(const void* ptr_a, const void* ptr_b);

/* Arranges tetrahedron indices to be in canonical order. */
static void tetrahedron_make_canonical(int* a, int* b, int* c, int *d);

/* Hash function for a tetrahedron defined by four indices. */
static size_t tetrahedron_hash(const void* ptr_tetrahedron);

/* Comparison function for a tetrahedron defined by four indices. */
static int tetrahedron_compare_equal(const void* ptr_a, const void* ptr_b);

// ===============================================================
//		Interpolation
// ===============================================================	

static int interpolate_3d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, int* d, double* u, double* v, double* w, double* x);

static int interpolate_2d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, double* u, double* v, double* w);

static int interpolate_1d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, double* u, double* v);

// ===============================================================
//		3D Triagulation
// ===============================================================	

static int triangulate_3d(TetrapalData* tetrapal);

static void rotate_tetrahedron(int* a, int* b, int* c, int* d, int rotations);

static void add_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d);

static void delete_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d);

static int adjacent_tetrahedron(const TetrapalData* tetrapal, int a, int b, int c);

static int check_conflict_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d, int e);

static void consider_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int e);

static int locate_tetrahedron(TetrapalData* tetrapal, int64_t p[3], int* a, int* b, int* c, int* d);

static int is_infinite_tetrahedron(int a, int b, int c, int d);

// ===============================================================
//		2D Triangulation
// ===============================================================	

static int triangulate_2d(TetrapalData* tetrapal);

static int locate_triangle(TetrapalData* tetrapal, int64_t p[3], int* a, int* b, int* c);

static void consider_triangle(TetrapalData* tetrapal, int a, int b, int e);

static int check_conflict_triangle(TetrapalData* tetrapal, int a, int b, int c, int e);

static void add_triangle(TetrapalData* tetrapal, int a, int b, int c);

static void delete_triangle(TetrapalData* tetrapal, int a, int b, int c);

static int adjacent_triangle(const TetrapalData* tetrapal, int a, int b);

static void rotate_triangle(int* a, int* b, int* c, int rotations);

static int is_infinite_triangle(int a, int b, int c);

// ===============================================================
//		1D Triangulation
// ===============================================================	

static int triangulate_1d(TetrapalData* tetrapal);

static int locate_segment(TetrapalData* tetrapal, int64_t p[3], int* a, int* b);

static void add_segment(TetrapalData* tetrapal, int a, int b);

static void delete_segment(TetrapalData* tetrapal, int a, int b);

static int adjacent_segment(const TetrapalData* tetrapal, int a, TetrapalEnum direction);

static int is_infinite_segment(int a, int b);

//*******************************************************************************************************************************
//		IMPLEMENTATION
//*******************************************************************************************************************************

// =========================================================
//		Geometric Predicates
// =========================================================

static double distance_squared(const int64_t a[3], const int64_t b[3])
{
	double ab[3] =
	{
		(double)a[0] - (double)b[0],
		(double)a[1] - (double)b[1],
		(double)a[2] - (double)b[2]
	};

	return ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
}

static int is_coincident(const int64_t a[3], const int64_t b[3])
{
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

static double orient1d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const double ba[3] = { (double)b[0] - (double)a[0], (double)b[1] - (double)a[1], (double)b[2] - (double)a[2] };
	const double ca[3] = { (double)c[0] - (double)a[0], (double)c[1] - (double)a[1], (double)c[2] - (double)a[2] };

	const double result =
		ca[0] * ba[0] +
		ca[1] * ba[1] +
		ca[2] * ba[2];

	return result;
}

static int orient1d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };

	const BigInt ba[3] = { bigint_sub(pb[0], pa[0]), bigint_sub(pb[1], pa[1]), bigint_sub(pb[2], pa[2]) };
	const BigInt ca[3] = { bigint_sub(pc[0], pa[0]), bigint_sub(pc[1], pa[1]), bigint_sub(pc[2], pa[2]) };

	const BigInt temp[3] = {
		bigint_mul(ca[0], ba[0]),
		bigint_mul(ca[1], ba[1]),
		bigint_mul(ca[2], ba[2]) };

	const BigInt result = bigint_add(bigint_add(temp[0], temp[1]), temp[2]);

	return result.sign;
}

static int orient1d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };

	const Interval ba[3] = { interval_sub(pb[0], pa[0]), interval_sub(pb[1], pa[1]), interval_sub(pb[2], pa[2]) };
	const Interval ca[3] = { interval_sub(pc[0], pa[0]), interval_sub(pc[1], pa[1]), interval_sub(pc[2], pa[2]) };

	const Interval temp[3] = {
		interval_mul(ca[0], ba[0]),
		interval_mul(ca[1], ba[1]),
		interval_mul(ca[2], ba[2]) };

	const Interval result = interval_add(interval_add(temp[0], temp[1]), temp[2]);

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return orient1d_exact(a, b, c);
}

static double orient2d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3])
{
	const double ba[3] = { (double)b[0] - (double)a[0], (double)b[1] - (double)a[1], (double)b[2] - (double)a[2] };
	const double ca[3] = { (double)c[0] - (double)a[0], (double)c[1] - (double)a[1], (double)c[2] - (double)a[2] };
	const double ed[3] = { (double)e[0] - (double)d[0], (double)e[1] - (double)d[1], (double)e[2] - (double)d[2] };
	const double fd[3] = { (double)f[0] - (double)d[0], (double)f[1] - (double)d[1], (double)f[2] - (double)d[2] };

	const double norm[3] = {
		ed[1] * fd[2] - ed[2] * fd[1],
		ed[2] * fd[0] - ed[0] * fd[2],
		ed[0] * fd[1] - ed[1] * fd[0] };

	const double cross[3] = {
		ba[1] * ca[2] - ba[2] * ca[1],
		ba[2] * ca[0] - ba[0] * ca[2],
		ba[0] * ca[1] - ba[1] * ca[0] };

	const double result =
		cross[0] * norm[0] +
		cross[1] * norm[1] +
		cross[2] * norm[2];

	return result;
}

static int orient2d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };
	const BigInt pd[3] = { bigint(d[0]), bigint(d[1]), bigint(d[2]) };
	const BigInt pe[3] = { bigint(e[0]), bigint(e[1]), bigint(e[2]) };
	const BigInt pf[3] = { bigint(f[0]), bigint(f[1]), bigint(f[2]) };

	const BigInt ba[3] = { bigint_sub(pb[0], pa[0]), bigint_sub(pb[1], pa[1]), bigint_sub(pb[2], pa[2]) };
	const BigInt ca[3] = { bigint_sub(pc[0], pa[0]), bigint_sub(pc[1], pa[1]), bigint_sub(pc[2], pa[2]) };
	const BigInt ed[3] = { bigint_sub(pe[0], pd[0]), bigint_sub(pe[1], pd[1]), bigint_sub(pe[2], pd[2]) };
	const BigInt fd[3] = { bigint_sub(pf[0], pd[0]), bigint_sub(pf[1], pd[1]), bigint_sub(pf[2], pd[2]) };

	const BigInt norm[3] = {
		bigint_sub(bigint_mul(ed[1], fd[2]), bigint_mul(ed[2], fd[1])),
		bigint_sub(bigint_mul(ed[2], fd[0]), bigint_mul(ed[0], fd[2])),
		bigint_sub(bigint_mul(ed[0], fd[1]), bigint_mul(ed[1], fd[0])) };

	const BigInt cross[3] = {
		bigint_sub(bigint_mul(ba[1], ca[2]), bigint_mul(ba[2], ca[1])),
		bigint_sub(bigint_mul(ba[2], ca[0]), bigint_mul(ba[0], ca[2])),
		bigint_sub(bigint_mul(ba[0], ca[1]), bigint_mul(ba[1], ca[0])) };

	const BigInt temp[3] = {
		bigint_mul(cross[0], norm[0]),
		bigint_mul(cross[1], norm[1]),
		bigint_mul(cross[2], norm[2]) };

	const BigInt result = bigint_add(bigint_add(temp[0], temp[1]), temp[2]);

	return result.sign;
}

static int orient2d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3], const int64_t f[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };
	const Interval pd[3] = { interval(d[0]), interval(d[1]), interval(d[2]) };
	const Interval pe[3] = { interval(e[0]), interval(e[1]), interval(e[2]) };
	const Interval pf[3] = { interval(f[0]), interval(f[1]), interval(f[2]) };

	const Interval ba[3] = { interval_sub(pb[0], pa[0]), interval_sub(pb[1], pa[1]), interval_sub(pb[2], pa[2]) };
	const Interval ca[3] = { interval_sub(pc[0], pa[0]), interval_sub(pc[1], pa[1]), interval_sub(pc[2], pa[2]) };
	const Interval ed[3] = { interval_sub(pe[0], pd[0]), interval_sub(pe[1], pd[1]), interval_sub(pe[2], pd[2]) };
	const Interval fd[3] = { interval_sub(pf[0], pd[0]), interval_sub(pf[1], pd[1]), interval_sub(pf[2], pd[2]) };

	const Interval norm[3] = {
		interval_sub(interval_mul(ed[1], fd[2]), interval_mul(ed[2], fd[1])),
		interval_sub(interval_mul(ed[2], fd[0]), interval_mul(ed[0], fd[2])),
		interval_sub(interval_mul(ed[0], fd[1]), interval_mul(ed[1], fd[0])) };

	const Interval cross[3] = {
		interval_sub(interval_mul(ba[1], ca[2]), interval_mul(ba[2], ca[1])),
		interval_sub(interval_mul(ba[2], ca[0]), interval_mul(ba[0], ca[2])),
		interval_sub(interval_mul(ba[0], ca[1]), interval_mul(ba[1], ca[0])) };

	const Interval temp[3] = {
		interval_mul(cross[0], norm[0]),
		interval_mul(cross[1], norm[1]),
		interval_mul(cross[2], norm[2]) };

	const Interval result = interval_add(interval_add(temp[0], temp[1]), temp[2]);

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return orient2d_exact(a, b, c, d, e, f);
}

static double orient3d_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const double ad[3] = { (double)a[0] - (double)d[0], (double)a[1] - (double)d[1], (double)a[2] - (double)d[2] };
	const double bd[3] = { (double)b[0] - (double)d[0], (double)b[1] - (double)d[1], (double)b[2] - (double)d[2] };
	const double cd[3] = { (double)c[0] - (double)d[0], (double)c[1] - (double)d[1], (double)c[2] - (double)d[2] };

	return
		ad[0] * (bd[1] * cd[2] - bd[2] * cd[1]) +
		bd[0] * (cd[1] * ad[2] - cd[2] * ad[1]) +
		cd[0] * (ad[1] * bd[2] - ad[2] * bd[1]);
}

static int orient3d_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };
	const BigInt pd[3] = { bigint(d[0]), bigint(d[1]), bigint(d[2]) };

	const BigInt ad[3] = { bigint_sub(pa[0], pd[0]), bigint_sub(pa[1], pd[1]), bigint_sub(pa[2], pd[2]) };
	const BigInt bd[3] = { bigint_sub(pb[0], pd[0]), bigint_sub(pb[1], pd[1]), bigint_sub(pb[2], pd[2]) };
	const BigInt cd[3] = { bigint_sub(pc[0], pd[0]), bigint_sub(pc[1], pd[1]), bigint_sub(pc[2], pd[2]) };

	const BigInt temp[3] = {
		bigint_mul(ad[0], bigint_sub(bigint_mul(bd[1], cd[2]), bigint_mul(bd[2], cd[1]))),
		bigint_mul(bd[0], bigint_sub(bigint_mul(cd[1], ad[2]), bigint_mul(cd[2], ad[1]))),
		bigint_mul(cd[0], bigint_sub(bigint_mul(ad[1], bd[2]), bigint_mul(ad[2], bd[1]))) };

	const BigInt result = bigint_add(bigint_add(temp[0], temp[1]), temp[2]);

	return result.sign;
}

static int orient3d_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };
	const Interval pd[3] = { interval(d[0]), interval(d[1]), interval(d[2]) };

	const Interval ad[3] = { interval_sub(pa[0], pd[0]), interval_sub(pa[1], pd[1]), interval_sub(pa[2], pd[2]) };
	const Interval bd[3] = { interval_sub(pb[0], pd[0]), interval_sub(pb[1], pd[1]), interval_sub(pb[2], pd[2]) };
	const Interval cd[3] = { interval_sub(pc[0], pd[0]), interval_sub(pc[1], pd[1]), interval_sub(pc[2], pd[2]) };

	const Interval temp[3] = {
		interval_mul(ad[0], interval_sub(interval_mul(bd[1], cd[2]), interval_mul(bd[2], cd[1]))),
		interval_mul(bd[0], interval_sub(interval_mul(cd[1], ad[2]), interval_mul(cd[2], ad[1]))),
		interval_mul(cd[0], interval_sub(interval_mul(ad[1], bd[2]), interval_mul(ad[2], bd[1]))) };

	const Interval result = interval_add(interval_add(temp[0], temp[1]), temp[2]);

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return orient3d_exact(a, b, c, d);
}
/* 
static double insegment_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const double ac[3] = { (double)a[0] - (double)c[0], (double)a[1] - (double)c[1], (double)a[2] - (double)c[2] };
	const double bc[3] = { (double)b[0] - (double)c[0], (double)b[1] - (double)c[1], (double)b[2] - (double)c[2] };

	const double result = ac[0] * bc[0] - ac[1] * bc[1] - ac[2] * bc[2];

	return result;
}
 */
static int insegment_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };

	const BigInt ac[3] = { bigint_sub(pa[0], pc[0]), bigint_sub(pa[1], pc[1]), bigint_sub(pa[2], pc[2]) };
	const BigInt bc[3] = { bigint_sub(pb[0], pc[0]), bigint_sub(pb[1], pc[1]), bigint_sub(pb[2], pc[2]) };

	const BigInt temp[3] = {
		bigint_mul(ac[0], bc[0]),
		bigint_mul(ac[1], bc[1]),
		bigint_mul(ac[2], bc[2]) };

	const BigInt result = bigint_sub(bigint_sub(temp[0], temp[1]), temp[2]);

	return result.sign;
}

static int insegment_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };

	const Interval ac[3] = { interval_sub(pa[0], pc[0]), interval_sub(pa[1], pc[1]), interval_sub(pa[2], pc[2]) };
	const Interval bc[3] = { interval_sub(pb[0], pc[0]), interval_sub(pb[1], pc[1]), interval_sub(pb[2], pc[2]) };

	const Interval temp[3] = {
		interval_mul(ac[0], bc[0]),
		interval_mul(ac[1], bc[1]),
		interval_mul(ac[2], bc[2]) };

	const Interval result = interval_sub(interval_sub(temp[0], temp[1]), temp[2]);

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return insegment_exact(a, b, c);
}
/* 
static double incircle_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const double ad[3] = { (double)a[0] - (double)d[0], (double)a[1] - (double)d[1], (double)a[2] - (double)d[2] };
	const double bd[3] = { (double)b[0] - (double)d[0], (double)b[1] - (double)d[1], (double)b[2] - (double)d[2] };
	const double cd[3] = { (double)c[0] - (double)d[0], (double)c[1] - (double)d[1], (double)c[2] - (double)d[2] };
	const double ba[3] = { (double)b[0] - (double)a[0], (double)b[1] - (double)a[1], (double)b[2] - (double)a[2] };
	const double ca[3] = { (double)c[0] - (double)a[0], (double)c[1] - (double)a[1], (double)c[2] - (double)a[2] };

	const double norm[3] = {
		ba[1] * ca[2] - ba[2] * ca[1],
		ba[2] * ca[0] - ba[0] * ca[2],
		ba[0] * ca[1] - ba[1] * ca[0] };

	const double ad_dot = ad[0] * ad[0] + ad[1] * ad[1] + ad[2] * ad[2];
	const double bd_dot = bd[0] * bd[0] + bd[1] * bd[1] + bd[2] * bd[2];
	const double cd_dot = cd[0] * cd[0] + cd[1] * cd[1] + cd[2] * cd[2];

	const double a_lift[3] = { ad[0] + norm[0] * ad_dot, ad[1] + norm[1] * ad_dot, ad[2] + norm[2] * ad_dot };
	const double b_lift[3] = { bd[0] + norm[0] * bd_dot, bd[1] + norm[1] * bd_dot, bd[2] + norm[2] * bd_dot };
	const double c_lift[3] = { cd[0] + norm[0] * cd_dot, cd[1] + norm[1] * cd_dot, cd[2] + norm[2] * cd_dot };

	const double result =
		c_lift[0] * (a_lift[1] * b_lift[2] - a_lift[2] * b_lift[1]) +
		c_lift[1] * (a_lift[2] * b_lift[0] - a_lift[0] * b_lift[2]) +
		c_lift[2] * (a_lift[0] * b_lift[1] - a_lift[1] * b_lift[0]);

	return result;
}
 */
static int incircle_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };
	const BigInt pd[3] = { bigint(d[0]), bigint(d[1]), bigint(d[2]) };

	const BigInt ad[3] = { bigint_sub(pa[0], pd[0]), bigint_sub(pa[1], pd[1]), bigint_sub(pa[2], pd[2]) };
	const BigInt bd[3] = { bigint_sub(pb[0], pd[0]), bigint_sub(pb[1], pd[1]), bigint_sub(pb[2], pd[2]) };
	const BigInt cd[3] = { bigint_sub(pc[0], pd[0]), bigint_sub(pc[1], pd[1]), bigint_sub(pc[2], pd[2]) };
	const BigInt ba[3] = { bigint_sub(pb[0], pa[0]), bigint_sub(pb[1], pa[1]), bigint_sub(pb[2], pa[2]) };
	const BigInt ca[3] = { bigint_sub(pc[0], pa[0]), bigint_sub(pc[1], pa[1]), bigint_sub(pc[2], pa[2]) };

	const BigInt norm[3] = {
		bigint_sub(bigint_mul(ba[1], ca[2]), bigint_mul(ba[2], ca[1])),
		bigint_sub(bigint_mul(ba[2], ca[0]), bigint_mul(ba[0], ca[2])),
		bigint_sub(bigint_mul(ba[0], ca[1]), bigint_mul(ba[1], ca[0])) };

	const BigInt dot[3] = {
		bigint_add(bigint_add(bigint_mul(ad[0], ad[0]), bigint_mul(ad[1], ad[1])), bigint_mul(ad[2], ad[2])),
		bigint_add(bigint_add(bigint_mul(bd[0], bd[0]), bigint_mul(bd[1], bd[1])), bigint_mul(bd[2], bd[2])),
		bigint_add(bigint_add(bigint_mul(cd[0], cd[0]), bigint_mul(cd[1], cd[1])), bigint_mul(cd[2], cd[2])) };

	const BigInt a_lift[3] = {
		bigint_add(bigint_mul(dot[0], norm[0]), ad[0]),
		bigint_add(bigint_mul(dot[0], norm[1]), ad[1]),
		bigint_add(bigint_mul(dot[0], norm[2]), ad[2]) };

	const BigInt b_lift[3] = {
		bigint_add(bigint_mul(dot[1], norm[0]), bd[0]),
		bigint_add(bigint_mul(dot[1], norm[1]), bd[1]),
		bigint_add(bigint_mul(dot[1], norm[2]), bd[2]) };

	const BigInt c_lift[3] = {
		bigint_add(bigint_mul(dot[2], norm[0]), cd[0]),
		bigint_add(bigint_mul(dot[2], norm[1]), cd[1]),
		bigint_add(bigint_mul(dot[2], norm[2]), cd[2]) };

	const BigInt temp[3] = {
		bigint_mul(c_lift[0], bigint_sub(bigint_mul(a_lift[1], b_lift[2]), bigint_mul(a_lift[2], b_lift[1]))),
		bigint_mul(c_lift[1], bigint_sub(bigint_mul(a_lift[2], b_lift[0]), bigint_mul(a_lift[0], b_lift[2]))),
		bigint_mul(c_lift[2], bigint_sub(bigint_mul(a_lift[0], b_lift[1]), bigint_mul(a_lift[1], b_lift[0]))) };

	const BigInt result = bigint_add(bigint_add(temp[0], temp[1]), temp[2]);

	return result.sign;
}

static int incircle_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };
	const Interval pd[3] = { interval(d[0]), interval(d[1]), interval(d[2]) };

	const Interval ad[3] = { interval_sub(pa[0], pd[0]), interval_sub(pa[1], pd[1]), interval_sub(pa[2], pd[2]) };
	const Interval bd[3] = { interval_sub(pb[0], pd[0]), interval_sub(pb[1], pd[1]), interval_sub(pb[2], pd[2]) };
	const Interval cd[3] = { interval_sub(pc[0], pd[0]), interval_sub(pc[1], pd[1]), interval_sub(pc[2], pd[2]) };
	const Interval ba[3] = { interval_sub(pb[0], pa[0]), interval_sub(pb[1], pa[1]), interval_sub(pb[2], pa[2]) };
	const Interval ca[3] = { interval_sub(pc[0], pa[0]), interval_sub(pc[1], pa[1]), interval_sub(pc[2], pa[2]) };

	const Interval norm[3] = {
		interval_sub(interval_mul(ba[1], ca[2]), interval_mul(ba[2], ca[1])),
		interval_sub(interval_mul(ba[2], ca[0]), interval_mul(ba[0], ca[2])),
		interval_sub(interval_mul(ba[0], ca[1]), interval_mul(ba[1], ca[0])) };

	const Interval dot[3] = {
		interval_add(interval_add(interval_mul(ad[0], ad[0]), interval_mul(ad[1], ad[1])), interval_mul(ad[2], ad[2])),
		interval_add(interval_add(interval_mul(bd[0], bd[0]), interval_mul(bd[1], bd[1])), interval_mul(bd[2], bd[2])),
		interval_add(interval_add(interval_mul(cd[0], cd[0]), interval_mul(cd[1], cd[1])), interval_mul(cd[2], cd[2])) };

	const Interval a_lift[3] = {
		interval_add(interval_mul(dot[0], norm[0]), ad[0]),
		interval_add(interval_mul(dot[0], norm[1]), ad[1]),
		interval_add(interval_mul(dot[0], norm[2]), ad[2]) };

	const Interval b_lift[3] = {
		interval_add(interval_mul(dot[1], norm[0]), bd[0]),
		interval_add(interval_mul(dot[1], norm[1]), bd[1]),
		interval_add(interval_mul(dot[1], norm[2]), bd[2]) };

	const Interval c_lift[3] = {
		interval_add(interval_mul(dot[2], norm[0]), cd[0]),
		interval_add(interval_mul(dot[2], norm[1]), cd[1]),
		interval_add(interval_mul(dot[2], norm[2]), cd[2]) };

	const Interval temp[3] = {
		interval_mul(c_lift[0], interval_sub(interval_mul(a_lift[1], b_lift[2]), interval_mul(a_lift[2], b_lift[1]))),
		interval_mul(c_lift[1], interval_sub(interval_mul(a_lift[2], b_lift[0]), interval_mul(a_lift[0], b_lift[2]))),
		interval_mul(c_lift[2], interval_sub(interval_mul(a_lift[0], b_lift[1]), interval_mul(a_lift[1], b_lift[0]))) };

	const Interval result = interval_add(interval_add(temp[0], temp[1]), temp[2]);

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return incircle_exact(a, b, c, d);
}
/* 
static double insphere_fast(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3])
{
	const double ea[3] = { (double)a[0] - (double)e[0], (double)a[1] - (double)e[1], (double)a[2] - (double)e[2] };
	const double eb[3] = { (double)b[0] - (double)e[0], (double)b[1] - (double)e[1], (double)b[2] - (double)e[2] };
	const double ec[3] = { (double)c[0] - (double)e[0], (double)c[1] - (double)e[1], (double)c[2] - (double)e[2] };
	const double ed[3] = { (double)d[0] - (double)e[0], (double)d[1] - (double)e[1], (double)d[2] - (double)e[2] };

	const double ab = ea[0] * eb[1] - eb[0] * ea[1];
	const double bc = eb[0] * ec[1] - ec[0] * eb[1];
	const double cd = ec[0] * ed[1] - ed[0] * ec[1];
	const double da = ed[0] * ea[1] - ea[0] * ed[1];
	const double ac = ea[0] * ec[1] - ec[0] * ea[1];
	const double bd = eb[0] * ed[1] - ed[0] * eb[1];

	const double abc = ea[2] * bc - eb[2] * ac + ec[2] * ab;
	const double bcd = eb[2] * cd - ec[2] * bd + ed[2] * bc;
	const double cda = ec[2] * da + ed[2] * ac + ea[2] * cd;
	const double dab = ed[2] * ab + ea[2] * bd + eb[2] * da;

	const double a_lift = ea[0] * ea[0] + ea[1] * ea[1] + ea[2] * ea[2];
	const double b_lift = eb[0] * eb[0] + eb[1] * eb[1] + eb[2] * eb[2];
	const double c_lift = ec[0] * ec[0] + ec[1] * ec[1] + ec[2] * ec[2];
	const double d_lift = ed[0] * ed[0] + ed[1] * ed[1] + ed[2] * ed[2];

	return (d_lift * abc - c_lift * dab) + (b_lift * cda - a_lift * bcd);
}
 */
static int insphere_exact(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3])
{
	const BigInt pa[3] = { bigint(a[0]), bigint(a[1]), bigint(a[2]) };
	const BigInt pb[3] = { bigint(b[0]), bigint(b[1]), bigint(b[2]) };
	const BigInt pc[3] = { bigint(c[0]), bigint(c[1]), bigint(c[2]) };
	const BigInt pd[3] = { bigint(d[0]), bigint(d[1]), bigint(d[2]) };
	const BigInt pe[3] = { bigint(e[0]), bigint(e[1]), bigint(e[2]) };

	const BigInt ae[3] = { bigint_sub(pa[0], pe[0]), bigint_sub(pa[1], pe[1]), bigint_sub(pa[2], pe[2]) };
	const BigInt be[3] = { bigint_sub(pb[0], pe[0]), bigint_sub(pb[1], pe[1]), bigint_sub(pb[2], pe[2]) };
	const BigInt ce[3] = { bigint_sub(pc[0], pe[0]), bigint_sub(pc[1], pe[1]), bigint_sub(pc[2], pe[2]) };
	const BigInt de[3] = { bigint_sub(pd[0], pe[0]), bigint_sub(pd[1], pe[1]), bigint_sub(pd[2], pe[2]) };

	const BigInt ab = bigint_sub(bigint_mul(ae[0], be[1]), bigint_mul(be[0], ae[1]));
	const BigInt bc = bigint_sub(bigint_mul(be[0], ce[1]), bigint_mul(ce[0], be[1]));
	const BigInt cd = bigint_sub(bigint_mul(ce[0], de[1]), bigint_mul(de[0], ce[1]));
	const BigInt da = bigint_sub(bigint_mul(de[0], ae[1]), bigint_mul(ae[0], de[1]));
	const BigInt ac = bigint_sub(bigint_mul(ae[0], ce[1]), bigint_mul(ce[0], ae[1]));
	const BigInt bd = bigint_sub(bigint_mul(be[0], de[1]), bigint_mul(de[0], be[1]));

	const BigInt abc = bigint_add(bigint_sub(bigint_mul(ae[2], bc), bigint_mul(be[2], ac)), bigint_mul(ce[2], ab));
	const BigInt bcd = bigint_add(bigint_sub(bigint_mul(be[2], cd), bigint_mul(ce[2], bd)), bigint_mul(de[2], bc));
	const BigInt cda = bigint_add(bigint_add(bigint_mul(ce[2], da), bigint_mul(de[2], ac)), bigint_mul(ae[2], cd));
	const BigInt dab = bigint_add(bigint_add(bigint_mul(de[2], ab), bigint_mul(ae[2], bd)), bigint_mul(be[2], da));

	const BigInt alift = bigint_add(bigint_add(bigint_mul(ae[0], ae[0]), bigint_mul(ae[1], ae[1])), bigint_mul(ae[2], ae[2]));
	const BigInt blift = bigint_add(bigint_add(bigint_mul(be[0], be[0]), bigint_mul(be[1], be[1])), bigint_mul(be[2], be[2]));
	const BigInt clift = bigint_add(bigint_add(bigint_mul(ce[0], ce[0]), bigint_mul(ce[1], ce[1])), bigint_mul(ce[2], ce[2]));
	const BigInt dlift = bigint_add(bigint_add(bigint_mul(de[0], de[0]), bigint_mul(de[1], de[1])), bigint_mul(de[2], de[2]));

	const BigInt result = bigint_add(
		bigint_sub(bigint_mul(dlift, abc), bigint_mul(clift, dab)),
		bigint_sub(bigint_mul(blift, cda), bigint_mul(alift, bcd)));

	return result.sign;
}

static int insphere_filtered(const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], const int64_t e[3])
{
	const Interval pa[3] = { interval(a[0]), interval(a[1]), interval(a[2]) };
	const Interval pb[3] = { interval(b[0]), interval(b[1]), interval(b[2]) };
	const Interval pc[3] = { interval(c[0]), interval(c[1]), interval(c[2]) };
	const Interval pd[3] = { interval(d[0]), interval(d[1]), interval(d[2]) };
	const Interval pe[3] = { interval(e[0]), interval(e[1]), interval(e[2]) };

	const Interval ae[3] = { interval_sub(pa[0], pe[0]), interval_sub(pa[1], pe[1]), interval_sub(pa[2], pe[2]) };
	const Interval be[3] = { interval_sub(pb[0], pe[0]), interval_sub(pb[1], pe[1]), interval_sub(pb[2], pe[2]) };
	const Interval ce[3] = { interval_sub(pc[0], pe[0]), interval_sub(pc[1], pe[1]), interval_sub(pc[2], pe[2]) };
	const Interval de[3] = { interval_sub(pd[0], pe[0]), interval_sub(pd[1], pe[1]), interval_sub(pd[2], pe[2]) };

	const Interval ab = interval_sub(interval_mul(ae[0], be[1]), interval_mul(be[0], ae[1]));
	const Interval bc = interval_sub(interval_mul(be[0], ce[1]), interval_mul(ce[0], be[1]));
	const Interval cd = interval_sub(interval_mul(ce[0], de[1]), interval_mul(de[0], ce[1]));
	const Interval da = interval_sub(interval_mul(de[0], ae[1]), interval_mul(ae[0], de[1]));
	const Interval ac = interval_sub(interval_mul(ae[0], ce[1]), interval_mul(ce[0], ae[1]));
	const Interval bd = interval_sub(interval_mul(be[0], de[1]), interval_mul(de[0], be[1]));

	const Interval abc = interval_add(interval_sub(interval_mul(ae[2], bc), interval_mul(be[2], ac)), interval_mul(ce[2], ab));
	const Interval bcd = interval_add(interval_sub(interval_mul(be[2], cd), interval_mul(ce[2], bd)), interval_mul(de[2], bc));
	const Interval cda = interval_add(interval_add(interval_mul(ce[2], da), interval_mul(de[2], ac)), interval_mul(ae[2], cd));
	const Interval dab = interval_add(interval_add(interval_mul(de[2], ab), interval_mul(ae[2], bd)), interval_mul(be[2], da));

	const Interval alift = interval_add(interval_add(interval_mul(ae[0], ae[0]), interval_mul(ae[1], ae[1])), interval_mul(ae[2], ae[2]));
	const Interval blift = interval_add(interval_add(interval_mul(be[0], be[0]), interval_mul(be[1], be[1])), interval_mul(be[2], be[2]));
	const Interval clift = interval_add(interval_add(interval_mul(ce[0], ce[0]), interval_mul(ce[1], ce[1])), interval_mul(ce[2], ce[2]));
	const Interval dlift = interval_add(interval_add(interval_mul(de[0], de[0]), interval_mul(de[1], de[1])), interval_mul(de[2], de[2]));

	const Interval result = interval_add(
		interval_sub(interval_mul(dlift, abc), interval_mul(clift, dab)),
		interval_sub(interval_mul(blift, cda), interval_mul(alift, bcd)));

	const int sign = interval_sign(result);

	if (sign != 0) return sign;
	return insphere_exact(a, b, c, d, e);
}

// =========================================================
//		Interval Arithmetic
// =========================================================

static inline void interval_swap_bounds(Interval* interval)
{
	double temp = interval->lower;
	interval->lower = interval->upper;
	interval->upper = temp;
}

static inline double interval_min(const double a, const double b)
{
	return a < b ? a : b;
}

static inline double interval_max(const double a, const double b)
{
	return a > b ? a : b;
}

static inline void interval_add_error(Interval* interval)
{
	//if (isnan(interval->lower) || isnan(interval->upper)) 
	//{
	//	interval->lower = interval->upper = nan(""); return;
	//}

	if (interval->lower > interval->upper)
		interval_swap_bounds(interval);

	if 		(interval->lower < DBL_MIN) interval->lower *= 1.0 + DBL_EPSILON;
	else if (interval->lower > DBL_MIN) interval->lower *= 1.0 - DBL_EPSILON;
	else interval->lower -= DBL_MIN * DBL_EPSILON;

	if 		(interval->upper < DBL_MIN) interval->upper *= 1.0 - DBL_EPSILON;
	else if (interval->upper > DBL_MIN) interval->upper *= 1.0 + DBL_EPSILON;
	else interval->upper += DBL_MIN * DBL_EPSILON;
}

static Interval interval(double value)
{
	Interval interval = { value, value };
	interval_add_error(&interval);
	return interval;
}

static Interval interval_add(const Interval lhs, const Interval rhs)
{
	Interval result;
	result.lower = lhs.lower + rhs.lower;
	result.upper = lhs.upper + rhs.upper;
	interval_add_error(&result);
	return result;
}

static Interval interval_sub(const Interval lhs, const Interval rhs)
{
	Interval result;
	result.lower = lhs.lower - rhs.upper;
	result.upper = lhs.upper - rhs.lower;
	interval_add_error(&result);
	return result;
}

static Interval interval_mul(const Interval lhs, const Interval rhs)
{
	Interval result;
	double aa = lhs.lower * rhs.lower;
	double ab = lhs.lower * rhs.upper;
	double ba = lhs.upper * rhs.lower;
	double bb = lhs.upper * rhs.upper;
	result.lower = interval_min(interval_min(aa, ab), interval_min(ba, bb));
	result.upper = interval_max(interval_max(aa, ab), interval_max(ba, bb));
	interval_add_error(&result);
	return result;
}

static int interval_sign(const Interval interval)
{
	if (interval.lower > 0.0 && interval.upper > 0.0) return 1;
	if (interval.lower < 0.0 && interval.upper < 0.0) return -1;
	return 0;
}

#ifdef TETRAPAL_DEBUG
static void interval_print(const Interval interval)
{

	printf("INTERVAL: %.17g : %.17g\n", interval.lower, interval.upper);
}
#endif


// =========================================================
//		Extended Precision Integer Arithmetic (BigInt)
// =========================================================

static BigInt bigint(int64_t number)
{
	BigInt bigint = bigint_zero();
	if (number == 0) return bigint;
	bigint.sign = number < 0 ? -1 : 1;
	bigint.digits[BIGINT_SIZE - 1] = (uint64_t)llabs(number) & 0x00000000ffffffff;
	bigint.digits[BIGINT_SIZE - 2] = ((uint64_t)llabs(number) & 0xffffffff00000000) >> 32;
	return bigint;
}

static BigInt bigint_zero()
{
	return (BigInt) { { 0 }, 0 };
}

static BigInt bigint_abs(const BigInt bigint)
{
	BigInt absolute = bigint;
	absolute.sign = absolute.sign != 0 ? 1 : 0;
	return absolute;
}

static BigInt bigint_neg(const BigInt bigint)
{
	BigInt negative = bigint;
	negative.sign = negative.sign != 0 ? -1 : 0;
	return negative;
}

static BigInt bigint_add(const BigInt lhs, const BigInt rhs)
{
	if (rhs.sign == 0) return lhs;
	if (lhs.sign == 0) return rhs;
	if (lhs.sign < rhs.sign) return bigint_sub(rhs, bigint_abs(lhs));
	if (lhs.sign > rhs.sign) return bigint_sub(lhs, bigint_abs(rhs));

	BigInt result = lhs;
	uint64_t sum, carry = 0;

	for (int i = BIGINT_SIZE - 1; i >= 0; i--)
	{
		sum = lhs.digits[i] + rhs.digits[i] + carry;
		carry = (sum & 0xffffffff00000000) >> 32;
		result.digits[i] = sum & 0x00000000ffffffff;
	}

	return result;
}

static BigInt bigint_sub(const BigInt lhs, const BigInt rhs)
{
	if (rhs.sign == 0) return lhs;
	if (lhs.sign < rhs.sign) return bigint_neg(bigint_add(rhs, bigint_abs(lhs)));
	if (lhs.sign > rhs.sign) return bigint_add(lhs, bigint_abs(rhs));
	if (lhs.sign < 0 && rhs.sign < 0) return bigint_sub(bigint_abs(rhs), bigint_abs(lhs));

	int comparison = bigint_compare(lhs, rhs);
	if (comparison < 0) return bigint_neg(bigint_sub(rhs, lhs));
	if (comparison == 0) return bigint_zero();

	BigInt result = lhs;
	uint64_t increment = 0, borrow;

	for (int i = BIGINT_SIZE - 1; i >= 0; i--)
	{
		borrow = 0;
		while (lhs.digits[i] + borrow < rhs.digits[i] + increment) borrow += (1uLL << 32);
		result.digits[i] = (lhs.digits[i] + borrow) - (rhs.digits[i] + increment);
		increment = borrow >> 32;
	}

	return result;
}

static BigInt bigint_mul(const BigInt lhs, const BigInt rhs)
{
	if (lhs.sign == 0 || rhs.sign == 0) return bigint_zero();

	BigInt addend, result = bigint_zero();
	uint64_t product, carry;

	for (int i = BIGINT_SIZE - 1; i >= 0; i--)
	{
		addend = (BigInt){ { 0 }, 1 };
		carry = 0;
		for (int j = BIGINT_SIZE - 1; j >= 0; j--)
		{
			product = lhs.digits[j] * rhs.digits[i] + carry;
			carry = (product & 0xffffffff00000000) >> 32;
			addend.digits[j] = product & 0x00000000ffffffff;
		}
		bigint_shift_l(&addend, BIGINT_SIZE - 1 - i);
		result = bigint_add(result, addend);
	}

	result.sign = lhs.sign == rhs.sign ? 1 : -1;
	return result;
}

static void bigint_shift_l(BigInt* bigint, unsigned int shift)
{
	if (shift == 0)
		return;

	if (shift > BIGINT_SIZE - 1) {
		*bigint = bigint_zero();
		return;
	}

	memmove(&bigint->digits[0], &bigint->digits[shift], sizeof(uint64_t) * (BIGINT_SIZE - shift));
	memset(&bigint->digits[BIGINT_SIZE - shift], 0, sizeof(uint64_t) * shift);
}

static int bigint_compare(const BigInt lhs, const BigInt rhs)
{
	if (lhs.sign == -1 && rhs.sign == -1) return bigint_compare(bigint_abs(rhs), bigint_abs(lhs));
	if (lhs.sign == 0 && rhs.sign == 0) return 0;
	if (lhs.sign < rhs.sign) return -1;
	if (lhs.sign > rhs.sign) return 1;

	int64_t difference;

	for (int i = 0; i < BIGINT_SIZE; i++)
	{
		difference = lhs.digits[i] - rhs.digits[i];
		if (difference != 0) return difference < 0 ? -1 * lhs.sign : 1 * lhs.sign;
	}

	return 0;
}

#ifdef TETRAPAL_DEBUG
static int bigint_is_zero(const BigInt bigint)
{
	for (int i = BIGINT_SIZE - 1; i >= 0; i--)
		if (bigint.digits[i] != 0)
			return 0;
	return 1;
}
#endif

#ifdef TETRAPAL_DEBUG
static void bigint_print_digits(const BigInt bigint)
{
	printf("BIG INT (BASE 2^32): ");
	if (bigint.sign < 0) printf("-");

	for (int i = 0; i < BIGINT_SIZE; i++) {
		printf("%llu|", bigint.digits[i]);
	}
	printf("\n");
}
#endif

#ifdef TETRAPAL_DEBUG
static void bigint_print_decimal(const BigInt bigint)
{
	printf("BIG INT (DECIMAL): ");
	if (bigint.sign < 0) printf("-");

	if (bigint.sign == 0) {
		printf("0\n"); return;
	}

	char remainder, buffer[BIGINT_SIZE * 16] = { 0 };
	BigInt dividend = bigint;
	uint64_t quotient, digit = 0;

	while (!bigint_is_zero(dividend))
		//while (dt_bigint_compare(dividend, bigint_zero()) != 0)
	{
		quotient = remainder = 0;
		for (int i = 0; i < BIGINT_SIZE; i++)
		{
			quotient = (dividend.digits[i] + remainder * 0x100000000) / 10;
			remainder = (dividend.digits[i] + remainder * 0x100000000) % 10;
			dividend.digits[i] = quotient;
		}
		buffer[digit++] = remainder + '0';
	}
	buffer[digit] = '\0';

	// Reverse the buffer
	int length = strlen(buffer);
	for (int y = 0; y < length / 2; y++)
	{
		char temp = buffer[y];
		buffer[y] = buffer[length - y - 1];
		buffer[length - y - 1] = temp;
	}

	printf("%s\n", buffer);
}
#endif

// =========================================================
//		Barycentric Interpolation
// =========================================================

static void barycentric_1d(const int64_t q[3], const int64_t a[3], const int64_t b[3], double* u, double* v)
{
	const double pq[3] = { (double)q[0], (double)q[1], (double)q[2] };
	const double pa[3] = { (double)a[0], (double)a[1], (double)a[2] };
	const double pb[3] = { (double)b[0], (double)b[1], (double)b[2] };

	const double ba[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] };
	const double qa[3] = { pq[0] - pa[0], pq[1] - pa[1], pq[2] - pa[2] };

	const double total = ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2];

	*v = (qa[0] * ba[0] + qa[1] * ba[1] + qa[2] * ba[2]) / total;
	//*v = (*v < 0.0) ? 0.0 : (*v > 1.0 ? 1.0 : *v); // Clamp
	*u = 1.0 - *v;
}

static void barycentric_2d(const int64_t q[3], const int64_t a[3], const int64_t b[3], const int64_t c[3], double* u, double* v, double* w)
{
	const double pq[3] = { (double)q[0], (double)q[1], (double)q[2] };
	const double pa[3] = { (double)a[0], (double)a[1], (double)a[2] };
	const double pb[3] = { (double)b[0], (double)b[1], (double)b[2] };
	const double pc[3] = { (double)c[0], (double)c[1], (double)c[2] };

	const double ba[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] };
	const double ca[3] = { pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2] };
	const double qa[3] = { pq[0] - pa[0], pq[1] - pa[1], pq[2] - pa[2] };

	const double abc[3] = { ba[1] * ca[2] - ba[2] * ca[1], ba[2] * ca[0] - ba[0] * ca[2], ba[0] * ca[1] - ba[1] * ca[0] };
	const double abq[3] = { ba[1] * qa[2] - ba[2] * qa[1], ba[2] * qa[0] - ba[0] * qa[2], ba[0] * qa[1] - ba[1] * qa[0] };
	const double aqc[3] = { qa[1] * ca[2] - qa[2] * ca[1], qa[2] * ca[0] - qa[0] * ca[2], qa[0] * ca[1] - qa[1] * ca[0] };

	const double total = abc[0] * abc[0] + abc[1] * abc[1] + abc[2] * abc[2];
	const double inverse = 1.0 / total;

	*w = (abq[0] * abc[0] + abq[1] * abc[1] + abq[2] * abc[2]) * inverse;
	*v = (aqc[0] * abc[0] + aqc[1] * abc[1] + aqc[2] * abc[2]) * inverse;
	*u = 1.0 - *v - *w;
}

static void barycentric_3d(const int64_t q[3], const int64_t a[3], const int64_t b[3], const int64_t c[3], const int64_t d[3], double* u, double* v, double* w, double* x)
{
	const double pq[3] = { (double)q[0], (double)q[1], (double)q[2] };
	const double pa[3] = { (double)a[0], (double)a[1], (double)a[2] };
	const double pb[3] = { (double)b[0], (double)b[1], (double)b[2] };
	const double pc[3] = { (double)c[0], (double)c[1], (double)c[2] };
	const double pd[3] = { (double)d[0], (double)d[1], (double)d[2] };

	const double ba[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] };
	const double ca[3] = { pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2] };
	const double da[3] = { pd[0] - pa[0], pd[1] - pa[1], pd[2] - pa[2] };
	const double qa[3] = { pq[0] - pa[0], pq[1] - pa[1], pq[2] - pa[2] };

	const double abc[3] = { ba[1] * ca[2] - ba[2] * ca[1], ba[2] * ca[0] - ba[0] * ca[2], ba[0] * ca[1] - ba[1] * ca[0] };
	const double adb[3] = { da[1] * ba[2] - da[2] * ba[1], da[2] * ba[0] - da[0] * ba[2], da[0] * ba[1] - da[1] * ba[0] };
	const double acd[3] = { ca[1] * da[2] - ca[2] * da[1], ca[2] * da[0] - ca[0] * da[2], ca[0] * da[1] - ca[1] * da[0] };

	const double total = da[0] * abc[0] + da[1] * abc[1] + da[2] * abc[2];
	const double inverse = 1.0 / total;

	*x = (qa[0] * abc[0] + qa[1] * abc[1] + qa[2] * abc[2]) * inverse;
	*w = (qa[0] * adb[0] + qa[1] * adb[1] + qa[2] * adb[2]) * inverse;
	*v = (qa[0] * acd[0] + qa[1] * acd[1] + qa[2] * acd[2]) * inverse;
	*u = 1.0 - *v - *w - *x;
}

// =========================================================
//		Spatial Grid
// =========================================================

static void spatialgrid_build(TetrapalData* tetrapal)
{
	uint32_t random_state = 1;

	if (tetrapal->dimensions == 3)
	{
		int a, b, c, d;
		int64_t pa[3], pb[3], pc[3], pd[3];

		int last[4] =
		{
			tetrapal->last[0],
			tetrapal->last[1],
			tetrapal->last[2],
			tetrapal->last[3],
		};

		for (int x = 0; x < SPATIALGRID_N; x++)
			for (int y = 0; y < SPATIALGRID_N; y++)
				for (int z = 0; z < SPATIALGRID_N; z++)
				{
					const int64_t step = TETRAPAL_PRECISION / SPATIALGRID_N;
					const int64_t p[3] = { x * step, y * step, z * step };

					a = last[0];
					b = last[1];
					c = last[2];
					d = last[3];

					// Start walking from within the triangulation
					while (1) {

						// If we arrive at an infinite tetrahedron, we must perform a convex hull walk
						if (is_infinite_tetrahedron(a, b, c, d)) {
							goto WALK_FACE;
						}

						// Randomly rotate indices (stochastic walk)
						rotate_tetrahedron(&a, &b, &c, &d, random_range_int(&random_state, 4));
						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);
						get_coordinates(tetrapal, c, pc);
						get_coordinates(tetrapal, d, pd);

						// Test the query point against every face
						if (orient3d_fast(pa, pb, pc, p) < 0) { d = adjacent_tetrahedron(tetrapal, c, b, a); swap_int(&a, &c); continue; }
						if (orient3d_fast(pa, pc, pd, p) < 0) { b = adjacent_tetrahedron(tetrapal, d, c, a); swap_int(&b, &d); continue; }
						if (orient3d_fast(pa, pd, pb, p) < 0) { c = adjacent_tetrahedron(tetrapal, b, d, a); swap_int(&c, &d); continue; }
						if (orient3d_fast(pb, pd, pc, p) < 0) { a = adjacent_tetrahedron(tetrapal, c, d, b); swap_int(&a, &d); continue; }

						// If none of the faces return negative, then we have found the enclosing tetrahedron
						goto WALK_END_3D;
					}

				WALK_FACE:
					while (1)
					{
						// Randomly rotate indices (stochastic walk)
						rotate_triangle(&a, &b, &c, random_range_int(&random_state, 3));
						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);
						get_coordinates(tetrapal, c, pc);

						// Check if the point lies in an incident edge region
						if (orient2d_fast(pa, pb, p, pa, pb, pc) < 0) { goto WALK_EDGE; }
						if (orient2d_fast(pb, pc, p, pb, pc, pa) < 0) { a = b; b = c; goto WALK_EDGE; }
						if (orient2d_fast(pc, pa, p, pc, pa, pb) < 0) { b = a; a = c; goto WALK_EDGE; }

						// Get the tetrahedron
						d = adjacent_tetrahedron(tetrapal, c, b, a);
						swap_int(&a, &c);

						goto WALK_END_3D;
					}

				WALK_EDGE:
					while (1)
					{
						int j, k;
						j = adjacent_tetrahedron(tetrapal, a, b, VERTEX_INFINITE);
						k = adjacent_tetrahedron(tetrapal, b, a, VERTEX_INFINITE);

						int64_t pa[3], pb[3], pj[3], pk[3];
						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);
						get_coordinates(tetrapal, j, pj);
						get_coordinates(tetrapal, k, pk);

						// Check if point lies in an incident vertex region
						if (orient1d_fast(pa, pb, p) < 0) { goto WALK_VERTEX; }
						if (orient1d_fast(pb, pa, p) < 0) { a = b; goto WALK_VERTEX; }

						// Check if point lies in an incident face region
						if (orient2d_fast(pa, pb, p, pb, pa, pj) < 0) { swap_int(&a, &b); c = j; goto WALK_FACE; }
						if (orient2d_fast(pb, pa, p, pa, pb, pk) < 0) { c = k; goto WALK_FACE; }

						// Get the tetrahedron
						c = adjacent_tetrahedron(tetrapal, a, b, VERTEX_INFINITE);
						d = adjacent_tetrahedron(tetrapal, a, b, c);

						goto WALK_END_3D;
					}

				WALK_VERTEX:
					while (1)
					{
						// Get all the hull edges connected to this vertex and check if the orthogonal projection of [p] lies on any of these edges
						for (size_t i = 0; i < graph_list_size(tetrapal->graph, a); i++) {

							b = graph_get(tetrapal->graph, a, i);

							// Only check vertices that lie on the convex hull
							if (adjacent_tetrahedron(tetrapal, a, b, VERTEX_INFINITE) == VERTEX_NULL)
								continue;

							int64_t pa[3], pb[3];
							get_coordinates(tetrapal, a, pa);
							get_coordinates(tetrapal, b, pb);

							// Check if point lies in an incident edge region
							if (orient1d_fast(pa, pb, p) > 0) { goto WALK_EDGE; }
						}

						// Get the tetrahedron
						for (size_t i = 0; i < graph_list_size(tetrapal->graph, a); i++)
						{
							b = graph_get(tetrapal->graph, a, i);
							c = adjacent_tetrahedron(tetrapal, a, b, VERTEX_INFINITE);
							if (c == VERTEX_NULL) continue;
							break;
						}

						d = adjacent_tetrahedron(tetrapal, a, b, c);
						goto WALK_END_3D;
					}

				WALK_END_3D:

					last[0] = tetrapal->grid.grid[x][y][z][0] = a;
					last[1] = tetrapal->grid.grid[x][y][z][1] = b;
					last[2] = tetrapal->grid.grid[x][y][z][2] = c;
					last[3] = tetrapal->grid.grid[x][y][z][3] = d;
				}
	}
	else if (tetrapal->dimensions == 2)
	{
		int a, b, c;
		int64_t pa[3], pb[3], pc[3];

		int last[3] =
		{
			tetrapal->last[0],
			tetrapal->last[1],
			tetrapal->last[2],
		};

		get_coordinates(tetrapal, last[0], pa);
		get_coordinates(tetrapal, last[1], pb);
		get_coordinates(tetrapal, last[2], pc);

		for (int x = 0; x < SPATIALGRID_N; x++)
			for (int y = 0; y < SPATIALGRID_N; y++)
				for (int z = 0; z < SPATIALGRID_N; z++)
				{
					const int64_t step = TETRAPAL_PRECISION / SPATIALGRID_N;
					const int64_t p[3] = { x * step, y * step, z * step };

					a = last[0];
					b = last[1];
					c = last[2];

					// Walk as normal, but use the fast predicates (numerical stability is not that important)
					while (1)
					{
						// If we arrive at an infinite triangle, then our point lies outside the current triangulation
						if (is_infinite_triangle(a, b, c))
						{
							break;
						}

						// Randomly rotate indices (stochastic walk)
						rotate_triangle(&a, &b, &c, random_range_int(&random_state, 3));
						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);
						get_coordinates(tetrapal, c, pc);

						// Test the query point against every face
						if (orient2d_fast(pa, pb, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { c = adjacent_triangle(tetrapal, b, a); swap_int(&a, &b); continue; }
						if (orient2d_fast(pb, pc, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { a = adjacent_triangle(tetrapal, c, b); swap_int(&c, &a); continue; }
						if (orient2d_fast(pc, pa, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { b = adjacent_triangle(tetrapal, a, c); swap_int(&c, &b); continue; }

						// If none of the faces return negative, then we have found the enclosing triangle
						goto WALK_END_2D;
					}

					// Walk the hull
					while (1)
					{
						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);

						int adjacent[2] = {
							adjacent_triangle(tetrapal, a, VERTEX_INFINITE),
							adjacent_triangle(tetrapal, VERTEX_INFINITE, b) };

						int64_t pj[3], pk[3];
						get_coordinates(tetrapal, adjacent[0], pj);
						get_coordinates(tetrapal, adjacent[1], pk);

						if ((orient2d_fast(pj, pa, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) > 0) && (orient1d_fast(pa, pj, p) > 0)) { b = adjacent[0]; swap_int(&a, &b); continue; }
						if ((orient2d_fast(pb, pk, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) > 0) && (orient1d_fast(pb, pk, p) > 0)) { a = adjacent[1]; swap_int(&a, &b); continue; }

						// This point is closest to this edge, so get the full triangle
						c = adjacent_triangle(tetrapal, b, a);
						swap_int(&a, &b);

						goto WALK_END_2D;
					}

				WALK_END_2D:

					last[0] = tetrapal->grid.grid[x][y][z][0] = a;
					last[1] = tetrapal->grid.grid[x][y][z][1] = b;
					last[2] = tetrapal->grid.grid[x][y][z][2] = c;
					tetrapal->grid.grid[x][y][z][3] = VERTEX_NULL;
				}
	}
	else if (tetrapal->dimensions == 1)
	{
		int a, b;
		int64_t pa[3], pb[3];

		int last[2] =
		{
			tetrapal->last[0],
			tetrapal->last[1],
		};

		for (int x = 0; x < SPATIALGRID_N; x++)
			for (int y = 0; y < SPATIALGRID_N; y++)
				for (int z = 0; z < SPATIALGRID_N; z++)
				{
					const int64_t step = TETRAPAL_PRECISION / SPATIALGRID_N;
					const int64_t p[3] = { x * step, y * step, z * step };

					a = last[0];
					b = last[1];

					while (1)
					{
						// If we arrive at an infinite segment, then our point lies beyond the boundary vertex
						if (is_infinite_segment(a, b))
						{
							// Get the finite edge connected to the finite vertex
							if (a == VERTEX_INFINITE)
							{
								a = adjacent_segment(tetrapal, b, VERTEX_NEXT);
							}
							else
							{
								b = adjacent_segment(tetrapal, a, VERTEX_PREV);
							}

							break;
						}

						get_coordinates(tetrapal, a, pa);
						get_coordinates(tetrapal, b, pb);

						// Test the query point against both sides of the line
						if (orient1d_fast(pa, pb, p) < 0) { b = adjacent_segment(tetrapal, a, VERTEX_PREV); swap_int(&a, &b); continue; }
						if (orient1d_fast(pb, pa, p) < 0) { a = adjacent_segment(tetrapal, b, VERTEX_NEXT); swap_int(&a, &b); continue; }

						// If none of the faces return negative, then we have found the enclosing segment
						break;
					}

					// End of the walk
					last[0] = tetrapal->grid.grid[x][y][z][0] = a;
					last[1] = tetrapal->grid.grid[x][y][z][1] = b;
					tetrapal->grid.grid[x][y][z][2] = VERTEX_NULL;
					tetrapal->grid.grid[x][y][z][3] = VERTEX_NULL;
				}
	}
}

static void spatialgrid_locate(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, int* d)
{
	const int64_t step = TETRAPAL_PRECISION / (SPATIALGRID_N - 1);
	const int64_t x = p[0] / step;
	const int64_t y = p[1] / step;
	const int64_t z = p[2] / step;
	*a = tetrapal->grid.grid[x][y][z][0];
	*b = tetrapal->grid.grid[x][y][z][1];
	*c = tetrapal->grid.grid[x][y][z][2];
	*d = tetrapal->grid.grid[x][y][z][3];
}

// ===============================================================
//	Vertex Graph
// ===============================================================

static Graph* graph_new(size_t size)
{
	Graph* graph = TETRAPAL_MALLOC(sizeof(*graph));

	if (graph == NULL) 
	{ 
		perror("GRAPH: MALLOC FAILED"); 
		return NULL;
	}

	graph->list = TETRAPAL_MALLOC(sizeof(*graph->list) * size);

	if (graph->list == NULL) 
	{ 
		perror("GRAPH: MALLOC FAILED"); 
		TETRAPAL_FREE(graph);
		return NULL;
	}

	graph->size = size;

	for (size_t i = 0; i < graph->size; i++) {

		graph->list[i].data = TETRAPAL_MALLOC(sizeof(*(graph->list->data)) * GRAPH_LIST_RESERVE);

		if ((graph->list[i].data) == NULL) 
		{ 
			perror("GRAPH: MALLOC FAILED"); 
			graph_free(graph);
			return NULL;
		}

		graph->list[i].capacity = GRAPH_LIST_RESERVE;
		graph->list[i].size = 0;
	}

	return graph;
}

static void graph_free(Graph* graph)
{
	if (graph == NULL) return;

	for (size_t i = 0; i < graph->size; i++)
	{
		if (graph->list[i].data != NULL)
			TETRAPAL_FREE(graph->list[i].data);
	}

	if (graph->list != NULL)
		TETRAPAL_FREE(graph->list);

	TETRAPAL_FREE(graph);
}

static void graph_list_insert(GraphList* list, int a)
{
	// Resize array if we reach capacity
	if (list->size == list->capacity) {

		list->capacity *= 2;

		int* temp = list->data;
		list->data = TETRAPAL_REALLOC(list->data, sizeof(*list->data) * list->capacity);

		if (list->data == NULL) 
		{ 
			perror("GRAPH: REALLOC FAILED"); 
			list->data = temp;
			return;
		}
	}

	list->data[list->size++] = a;
}
/* 
static void graph_list_remove(GraphList* list, int a)
{
	for (size_t i = 0; i < list->size; i++) if (list->data[i] == a)
	{
		list->data[i] = list->data[--list->size]; return;
	}
}
 */
static void graph_insert(Graph* graph, int a, int b)
{
	// Validity check
	if ((size_t)a >= graph->size || (size_t)b >= graph->size || a < 0 || b < 0)
		return;

	// Ensure the vertices are not already linked
	for (size_t i = 0; i < graph->list[a].size; i++)
		if (graph->list[a].data[i] == b) return;

	graph_list_insert(&graph->list[a], b);
	graph_list_insert(&graph->list[b], a);
}
/* 
static void graph_remove(Graph* graph, int a, int b)
{
	// Validity check
	if ((size_t)a >= graph->size || (size_t)b >= graph->size || a < 0 || b < 0)
		return;

	graph_list_remove(&graph->list[a], b);
	graph_list_remove(&graph->list[b], a);
}
 */
static size_t graph_list_size(const Graph* graph, const size_t vertex)
{
	return graph->list[vertex].size;
}

static int graph_get(const Graph* graph, const size_t vertex, const size_t index)
{
	return graph->list[vertex].data[index];
}

static int graph_find_number_of_unique_edges(const Graph* graph)
{
	int count = 0;

	for (size_t i = 0; i < graph->size; i++)
	{
		count += graph->list[i].size;
	}

	// Each segment is recorded twice
	return count / 2;
}

// ===============================================================
//	Hash Map 
// ===============================================================	

static HashMap* hashmap_new(const size_t stride_key, const size_t stride_value, size_t reserve, const HashMap_Func_Hash func_hash, const HashMap_Func_Equal func_equal)
{
	HashMap* hashmap = TETRAPAL_MALLOC(sizeof(*hashmap));

	if (hashmap == NULL) 
		return NULL;

	if (reserve < 1) reserve = 1;

	hashmap->capacity = reserve;
	hashmap->stride_key = stride_key;
	hashmap->stride_value = stride_value;
	hashmap->size = 0;
	hashmap->func_hash = func_hash;
	hashmap->func_equal = func_equal;
	hashmap->table.capacity = (size_t)(reserve * HASHMAP_GROWTH_FACTOR) + 1;

	hashmap->data = TETRAPAL_MALLOC((stride_key + stride_value) * hashmap->capacity);

	if (!hashmap->data) 
	{ 
		TETRAPAL_FREE(hashmap);
		return NULL;
	}

	hashmap->table.buckets = TETRAPAL_MALLOC(sizeof(*hashmap->table.buckets) * hashmap->table.capacity);
	
	if (!hashmap->table.buckets) 
	{ 
		TETRAPAL_FREE(hashmap->data);
		TETRAPAL_FREE(hashmap);
		return NULL;
	}

	// Allocate buckets
	for (size_t i = 0; i < hashmap->table.capacity; i++) 
	{
		HashMap_Bucket* bucket = TETRAPAL_MALLOC(sizeof(*bucket) + sizeof(bucket->data) * HASHMAP_BUCKET_RESERVE);

		if (bucket == NULL) 
		{ 
			perror("HASHMAP: MALLOC FAILED"); 
			hashmap_free(hashmap);
			return NULL;
		}

		bucket->capacity = HASHMAP_BUCKET_RESERVE;
		bucket->size = 0;
		hashmap->table.buckets[i] = bucket;
	}

	return hashmap;
}

static void hashmap_free(HashMap* hashmap)
{
	if (!hashmap) return;
	
	for (size_t i = 0; i < hashmap->table.capacity; i++)
	{
		if (hashmap->table.buckets[i] != NULL)
			TETRAPAL_FREE(hashmap->table.buckets[i]);
	}

	if (hashmap->table.buckets != NULL)
		TETRAPAL_FREE(hashmap->table.buckets);
	
	if (hashmap->data != NULL)
		TETRAPAL_FREE(hashmap->data);

	TETRAPAL_FREE(hashmap);
}

static void hashmap_insert(HashMap* hashmap, const void* key, const void* value)
{
	// Grow and rehash the map if we hit the load factor
	if (hashmap->size / (double)hashmap->table.capacity > HASHMAP_LOAD_FACTOR)
		hashmap_table_grow(hashmap);

	const size_t hash = hashmap->func_hash(key);
	const size_t index = hash % hashmap->table.capacity;
	HashMap_Bucket* bucket = hashmap->table.buckets[index];

	// Iterate through the bucket until we get to the end, ensuring the element doesn't already exist.
	for (size_t i = 0; i < bucket->size; i++) 
	{
		const size_t elem_index = (&bucket->data)[i];
		const void* current_key = hashmap_get_key(hashmap, elem_index);

		if (hashmap->func_equal(key, current_key))
			return;
	}

	// Realloc element arrays if it is at capacity
	if (hashmap->size == hashmap->capacity) 
	{
		hashmap->capacity = (size_t)(hashmap->capacity * HASHMAP_GROWTH_FACTOR) + 1;

		void* temp = hashmap->data;
		hashmap->data = TETRAPAL_REALLOC(hashmap->data, (hashmap->stride_key + hashmap->stride_value) * hashmap->capacity);
		
		// If we're out of memory, do nothing
		if (hashmap->data == NULL) 
		{ 
			perror("HASHMAP: REALLOC FAILED"); 
			hashmap->data = temp; // Restore pointer to old data
			return; 
		}
	}

	// Realloc the bucket if it is at capacity
	if (bucket->size == bucket->capacity) 
	{
		bucket->capacity = (size_t)(bucket->capacity * HASHMAP_GROWTH_FACTOR) + 1;
		void* temp = bucket;
		bucket = TETRAPAL_REALLOC(bucket, sizeof(*bucket) + sizeof(bucket->data) * bucket->capacity);

		// If we're out of memory, do nothing
		if (bucket == NULL) 
		{ 
			perror("HASHMAP: REALLOC FAILED"); 
			bucket = temp; // Restore pointer to old data
			return;  
		}

		hashmap->table.buckets[index] = bucket; // Prevent dangling pointer!
	}

	// Add the element to the element arrays
	memcpy(((char*)hashmap->data) + (hashmap->stride_key + hashmap->stride_value) * hashmap->size, key, hashmap->stride_key);
	memcpy(((char*)hashmap->data) + (hashmap->stride_key + hashmap->stride_value) * hashmap->size + hashmap->stride_key, value, hashmap->stride_value);
	hashmap->size++;

	// Add index to the bucket
	(&bucket->data)[bucket->size] = hashmap->size - 1;
	bucket->size++;
}

static void hashmap_remove(HashMap* hashmap, const void* key)
{
	size_t hash = hashmap->func_hash(key);
	size_t index = hash % hashmap->table.capacity;
	HashMap_Bucket* bucket = hashmap->table.buckets[index];

	// Iterate through the bucket until we find the element, if it exists
	for (size_t i = 0; i < bucket->size; i++) {

		const size_t elem_index = (&bucket->data)[i];
		void* current_key = hashmap_get_key(hashmap, elem_index);

		if (hashmap->func_equal(key, current_key)) { // Element was found

			// Replace the elements with the last and decrease size
			if (elem_index != hashmap->size - 1)
				memcpy(current_key, hashmap_get_key(hashmap, hashmap->size - 1), hashmap->stride_key + hashmap->stride_value);
			hashmap->size--;

			// Do the same for the bucket entry
			(&bucket->data)[i] = (&bucket->data)[bucket->size - 1];
			bucket->size--;

			// Now find the element that was moved and update the index
			HashMap_Bucket* other_bucket = hashmap->table.buckets[hashmap->func_hash(current_key) % hashmap->table.capacity];

			for (size_t i = 0; i < other_bucket->size; i++) {

				if ((&other_bucket->data)[i] == hashmap->size) {
					(&other_bucket->data)[i] = elem_index;
					return;
				}

			}
		}
	}
}

/*
static int hashmap_has(const HashMap* hashmap, const void* key)
{
	size_t hash = hashmap->func_hash(key);
	size_t index = hash % hashmap->table.capacity;
	HashMap_Bucket* bucket = hashmap->table.buckets[index];

	for (size_t i = 0; i < bucket->size; i++) {
		const size_t elem_index = (&bucket->data)[i];
		const void* current_key = hashmap_get_key(hashmap, elem_index);
		if (hashmap->func_equal(key, current_key))
			return 1;
	}

	return 0;
}
*/

static void* hashmap_find(const HashMap* hashmap, const void* key)
{
	size_t hash = hashmap->func_hash(key);
	size_t index = hash % hashmap->table.capacity;
	HashMap_Bucket* bucket = hashmap->table.buckets[index];

	for (size_t i = 0; i < bucket->size; i++) {
		const size_t elem_index = (&bucket->data)[i];
		const void* current_key = hashmap_get_key(hashmap, elem_index);
		if (hashmap->func_equal(key, current_key))
			return hashmap_get_value(hashmap, elem_index);
	}

	return NULL;
}

static void* hashmap_get_key(const HashMap* hashmap, const size_t index)
{
	//if (index >= hashmap->size) return NULL;
	return ((char*)hashmap->data) + (hashmap->stride_key + hashmap->stride_value) * index;
}

static void* hashmap_get_value(const HashMap* hashmap, const size_t index)
{
	//if (index >= hashmap->size) return NULL;
	return ((char*)hashmap->data) + (hashmap->stride_key + hashmap->stride_value) * index + hashmap->stride_key;
}

static size_t hashmap_size(const HashMap* hashmap)
{
	return hashmap->size;
}

static void hashmap_table_grow(HashMap* hashmap)
{
	// Clear the buckets
	const size_t old_capacity = hashmap->table.capacity;

	for (size_t i = 0; i < old_capacity; i++)
		hashmap->table.buckets[i]->size = 0;

	// Realloc the bucket array
	void* old_buckets = hashmap->table.buckets;
	hashmap->table.capacity = (size_t)(hashmap->table.capacity * HASHMAP_GROWTH_FACTOR) + 1;
	hashmap->table.buckets = TETRAPAL_REALLOC(hashmap->table.buckets, sizeof(*hashmap->table.buckets) * hashmap->table.capacity);

	if (hashmap->table.buckets == NULL) 
	{ 
		perror("HASHMAP: REALLOC FAILED"); 
		hashmap->table.buckets = old_buckets; // Restore pointer to old data
		return;
	}

	// Allocate new buckets
	for (size_t i = old_capacity; i < hashmap->table.capacity; i++) 
	{
		HashMap_Bucket* bucket = TETRAPAL_MALLOC(sizeof(*bucket) + sizeof(bucket->data) * HASHMAP_BUCKET_RESERVE);

		if (!bucket) 
		{ 
			perror("HASHMAP: REALLOC FAILED"); 
			exit(1); 
		}

		bucket->capacity = HASHMAP_BUCKET_RESERVE;
		bucket->size = 0;
		hashmap->table.buckets[i] = bucket;
	}

	// Rehash all elements by iterating through the element array
	for (size_t i = 0; i < hashmap->size; i++) {

		const void* key = hashmap_get_key(hashmap, i);
		size_t hash = hashmap->func_hash(key);
		size_t index = hash % hashmap->table.capacity;
		HashMap_Bucket* bucket = hashmap->table.buckets[index];

		// Realloc the bucket if it is at capacity
		if (bucket->size == bucket->capacity) 
		{
			bucket->capacity = (size_t)(bucket->capacity * HASHMAP_GROWTH_FACTOR) + 1;
			HashMap_Bucket* temp = bucket;

			bucket = TETRAPAL_REALLOC(bucket, sizeof(*bucket) + sizeof(bucket->data) * bucket->capacity);

			if (bucket == NULL) 
			{ 
				perror("HASHMAP: REALLOC FAILED"); 
				bucket = temp; // Restore pointer to old data
				return;
			}

			hashmap->table.buckets[index] = bucket; // Prevent dangling pointer!
		}

		// Add index to the bucket
		(&bucket->data)[bucket->size] = i;
		bucket->size++;
	}
}

// ===============================================================
//	Hash Set 
// ===============================================================	

static HashSet* hashset_new(const size_t stride, size_t reserve, const HashSet_Func_Hash func_hash, const HashSet_Func_Equal func_equal)
{
	HashSet* hashset = TETRAPAL_MALLOC(sizeof(*hashset));
	if (!hashset) { perror("MAP: MALLOC ERROR"); exit(1); }

	if (reserve < 1) reserve = 1;

	hashset->capacity = reserve;
	hashset->stride = stride;
	hashset->size = 0;
	hashset->func_hash = func_hash;
	hashset->func_equal = func_equal;
	hashset->map.capacity = (size_t)(reserve * HASHSET_GROWTH_FACTOR) + 1;

	hashset->data = TETRAPAL_MALLOC(stride * hashset->capacity);
	if (!hashset->data) { perror("MAP: MALLOC ERROR"); exit(1); }

	hashset->map.buckets = TETRAPAL_MALLOC(sizeof(*hashset->map.buckets) * hashset->map.capacity);
	if (!hashset->map.buckets) { perror("MAP: MALLOC ERROR"); exit(1); }

	// Allocate buckets
	for (size_t i = 0; i < hashset->map.capacity; i++) {
		hashset->map.buckets[i] = TETRAPAL_MALLOC(sizeof(*hashset->map.buckets[i]) + sizeof(hashset->map.buckets[i]->data) * HASHSET_BUCKET_RESERVE);
		if (!hashset->map.buckets[i]) { perror("MAP: MALLOC ERROR"); exit(1); }
		hashset->map.buckets[i]->capacity = HASHSET_BUCKET_RESERVE;
		hashset->map.buckets[i]->size = 0;
	}

	return hashset;
}

static void hashset_free(HashSet* hashset)
{
	if (!hashset) return;

	for (size_t i = 0; i < hashset->map.capacity; i++)
		TETRAPAL_FREE(hashset->map.buckets[i]);

	TETRAPAL_FREE(hashset->map.buckets);
	TETRAPAL_FREE(hashset->data);
	TETRAPAL_FREE(hashset);
}

static void hashset_insert(HashSet* hashset, const void* element)
{
	size_t hash = hashset->func_hash(element);
	size_t index = hash % hashset->map.capacity;
	HashSet_Bucket* bucket = hashset->map.buckets[index];

	// Iterate through the bucket until we get to the end, ensuring the element doesn't already exist.
	for (size_t i = 0; i < bucket->size; i++) {
		const size_t elem_index = (&bucket->data)[i];
		const void* current = hashset_get(hashset, elem_index);
		if (hashset->func_equal(element, current))
			return;
	}

	// Realloc the element array if it is at capacity
	if (hashset->size == hashset->capacity) {
		hashset->capacity = (size_t)(hashset->capacity * HASHSET_GROWTH_FACTOR) + 1;
		void* temp = hashset->data;
		hashset->data = TETRAPAL_REALLOC(hashset->data, hashset->stride * hashset->capacity);
		if (!hashset->data) { perror("MAP: REALLOC ERROR"); TETRAPAL_FREE(temp); exit(1); }
	}

	// Realloc the bucket if it is at capacity
	if (bucket->size == bucket->capacity) {
		bucket->capacity = (size_t)(bucket->capacity * HASHSET_GROWTH_FACTOR) + 1;
		void* temp = bucket;
		bucket = TETRAPAL_REALLOC(bucket, sizeof(*bucket) + sizeof(bucket->data) * bucket->capacity);
		if (!bucket) { perror("MAP: REALLOC ERROR"); TETRAPAL_FREE(temp); exit(1); }
		hashset->map.buckets[index] = bucket; // Prevent dangling pointer!
	}

	// Add the element to the element array
	memcpy(((char*)hashset->data) + hashset->stride * hashset->size, element, hashset->stride);
	hashset->size++;

	// Add index to the bucket
	(&bucket->data)[bucket->size] = hashset->size - 1;
	bucket->size++;

	// Grow and rehash the map if we hit the load factor
	if (hashset->size / (double)hashset->map.capacity > HASHSET_LOAD_FACTOR)
		hashset_map_grow(hashset);
}

/* 
static void hashset_remove(HashSet* hashset, const void* element)
{
	size_t hash = hashset->func_hash(element);
	size_t index = hash % hashset->map.capacity;
	HashSet_Bucket* bucket = hashset->map.buckets[index];

	// Iterate through the bucket until we find the element, if it exists
	for (size_t i = 0; i < bucket->size; i++) {

		const size_t elem_index = (&bucket->data)[i];
		void* current = ((char*)hashset->data) + hashset->stride * elem_index;

		if (hashset->func_equal(element, current)) { // Element was found

			// Replacing the element with the last element and decrease size
			memcpy(current, hashset_end(hashset), hashset->stride);
			hashset->size--;

			// Do the same for the bucket entry
			(&bucket->data)[i] = (&bucket->data)[bucket->size - 1];
			bucket->size--;

			// Now find the element that was moved and update the index
			HashSet_Bucket* other_bucket = hashset->map.buckets[hashset->func_hash(current) % hashset->map.capacity];

			for (size_t i = 0; i < other_bucket->size; i++) {
				if ((&other_bucket->data)[i] == hashset->size) {
					(&other_bucket->data)[i] = elem_index;
					return;
				}
			}
		}
	}
}

static int hashset_has(const HashSet* hashset, const void* element)
{
	size_t hash = hashset->func_hash(element);
	size_t index = hash % hashset->map.capacity;
	HashSet_Bucket* bucket = hashset->map.buckets[index];

	for (size_t i = 0; i < bucket->size; i++) {
		const size_t elem_index = (&bucket->data)[i];
		const void* current = hashset_get(hashset, elem_index);
		if (hashset->func_equal(element, current))
			return 1;
	}

	return 0;
}
 */

static const void* hashset_get(const HashSet* hashset, const size_t index)
{
	//if (index >= hashset->size) return NULL;
	return ((char*)hashset->data) + hashset->stride * index;
}

static const void* hashset_begin(const HashSet* hashset)
{
	return ((char*)hashset->data);
}
/* 
static const void* hashset_end(const HashSet* hashset)
{
	return ((char*)hashset->data) + hashset->stride * (hashset->size - 1);
}
 */
static size_t hashset_size(const HashSet* hashset)
{
	return hashset->size;
}

static void hashset_map_grow(HashSet* hashset)
{
	// Clear the buckets
	const size_t old_capacity = hashset->map.capacity;
	hashset->map.capacity = (size_t)(hashset->map.capacity * HASHSET_GROWTH_FACTOR) + 1;

	for (size_t i = 0; i < old_capacity; i++)
		hashset->map.buckets[i]->size = 0;

	// Realloc the bucket array
	void* old_buckets = hashset->map.buckets;
	hashset->map.buckets = TETRAPAL_REALLOC(hashset->map.buckets, sizeof(*hashset->map.buckets) * hashset->map.capacity);
	if (!hashset->map.buckets) { perror("MAP: REALLOC ERROR"); TETRAPAL_FREE(old_buckets); exit(1); }

	// Allocate new buckets
	for (size_t i = old_capacity; i < hashset->map.capacity; i++) {
		hashset->map.buckets[i] = TETRAPAL_MALLOC(sizeof(*hashset->map.buckets[i]) + sizeof(hashset->map.buckets[i]->data) * HASHSET_BUCKET_RESERVE);
		if (!hashset->map.buckets[i]) { perror("MAP: MALLOC ERROR"); exit(1); }
		hashset->map.buckets[i]->capacity = HASHSET_BUCKET_RESERVE;
		hashset->map.buckets[i]->size = 0;
	}

	// Rehash all elements by iterating through the element array
	for (size_t i = 0; i < hashset->size; i++) {

		const void* element = ((char*)hashset->data) + hashset->stride * i;

		size_t hash = hashset->func_hash(element);
		size_t index = hash % hashset->map.capacity;
		HashSet_Bucket* bucket = hashset->map.buckets[index];

		// Realloc the bucket if it is at capacity
		if (bucket->size == bucket->capacity) {
			bucket->capacity = (size_t)(bucket->capacity * HASHSET_GROWTH_FACTOR) + 1;
			void* temp = bucket;
			bucket = TETRAPAL_REALLOC(bucket, sizeof(*bucket) + sizeof(bucket->data) * bucket->capacity);
			if (!bucket) { perror("MAP: REALLOC ERROR"); TETRAPAL_FREE(temp); exit(1); }
			hashset->map.buckets[index] = bucket; // Prevent dangling pointer!
		}

		// Add index to the bucket
		(&bucket->data)[bucket->size] = i;
		bucket->size++;
	}
}

// ===============================================================
//	Tetrapal Core
// ===============================================================

TetrapalData* tetrapal_new(const float *points, const int size)
{
	// If no points are given then fail gracefully
	if (size <= 0)
	{
		return NULL;
	}

	TetrapalData* tetrapal = TETRAPAL_MALLOC(sizeof(*tetrapal));

	if (tetrapal == NULL)
	{
		return NULL;
	}

	// Begin initialising struct data
	tetrapal->size = size;
	tetrapal->dimensions = 0;
	tetrapal->vertices = TETRAPAL_MALLOC(sizeof(*tetrapal->vertices) * size);
	tetrapal->adjacency = NULL;
	tetrapal->graph = NULL;
	tetrapal->number_of_segments = 0;
	tetrapal->number_of_triangles = 0;
	tetrapal->number_of_tetrahedra = 0;
	tetrapal->orient_a[0] = tetrapal->orient_a[1] = tetrapal->orient_a[2] = 0;
	tetrapal->orient_b[0] = tetrapal->orient_b[1] = tetrapal->orient_b[2] = 0;
	tetrapal->orient_c[0] = tetrapal->orient_c[1] = tetrapal->orient_c[2] = 0;

	if (tetrapal->vertices == NULL)
	{
		tetrapal_free(tetrapal);
		return NULL;
	}

	// Add all vertices
	for (int i = 0; i < size; i++)
	{
		// Clamp the input between 0.0 and 1.0
		float px = points[i * 3 + 0];
		float py = points[i * 3 + 1];
		float pz = points[i * 3 + 2];

		clamp_float(&px, 0.0f, 1.0f);
		clamp_float(&py, 0.0f, 1.0f);
		clamp_float(&pz, 0.0f, 1.0f);

		tetrapal->vertices[i].x = (double)px * TETRAPAL_PRECISION;
		tetrapal->vertices[i].y = (double)py * TETRAPAL_PRECISION;
		tetrapal->vertices[i].z = (double)pz * TETRAPAL_PRECISION;
	}

	// A 0-dimensional triangulation requires no further processing (only one vertex, spatial/adjacency structures are not needed)
	if (size == 1)
	{
		tetrapal->dimensions = 0;
		return tetrapal;
	}

	// Initialise adjacency structures
	tetrapal->adjacency = hashmap_new(sizeof(TetrapalAdjacencyKey), sizeof(int), size * 7, adjacency_hash, adjacency_compare_equal);
	tetrapal->graph = graph_new(size);

	if (tetrapal->adjacency == NULL || tetrapal->graph == NULL)
	{
		tetrapal_free(tetrapal);
		return NULL;
	}

	// Now we can start the main triangulation algorithm, repeatedly attempting decreasing dimensions until we find the appropriate one
	if (triangulate_3d(tetrapal) == 0)
	{
		if (triangulate_2d(tetrapal) == 0)
		{
			if (triangulate_1d(tetrapal) == 0)
			{
				// If all attempts failed for whatever reason, abort
				tetrapal_free(tetrapal);
				return NULL;
			}
		}
	}

	// Cache the number of segments
	tetrapal->number_of_segments = graph_find_number_of_unique_edges(tetrapal->graph);

	// Success!!
	return tetrapal;
}

void tetrapal_free(TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
	{
		return;
	}

	if (tetrapal->vertices != NULL)
	{
		TETRAPAL_FREE(tetrapal->vertices);
	}

	// These functions should have their own NULL checks
	graph_free(tetrapal->graph);
	hashmap_free(tetrapal->adjacency);

	TETRAPAL_FREE(tetrapal);
}

int tetrapal_interpolate(const TetrapalData* tetrapal, const float point[3],
	int* a, int* b, int* c, int* d,
	double* u, double* v, double* w, double* x)
{
	*a = *b = *c = *d = 0;
	*u = *v = *w = *x = 0.0;

	// Clamp the input between 0.0 and 1.0
	float px = point[0];
	float py = point[1];
	float pz = point[2];

	clamp_float(&px, 0.0f, 1.0f);
	clamp_float(&py, 0.0f, 1.0f);
	clamp_float(&pz, 0.0f, 1.0f);

	const int64_t p[3] = {
		(double)px * TETRAPAL_PRECISION,
		(double)py * TETRAPAL_PRECISION,
		(double)pz * TETRAPAL_PRECISION };

	if (!tetrapal) return 0;
	switch (tetrapal->dimensions) {
	case 0: *a = 1; *u = 1.0; return 1;
	case 1: return interpolate_1d(tetrapal, p, a, b, u, v);
	case 2: return interpolate_2d(tetrapal, p, a, b, c, u, v, w);
	case 3: return interpolate_3d(tetrapal, p, a, b, c, d, u, v, w, x);
	default: return 0;
	}
}

int tetrapal_nearest_neighbour(const TetrapalData* tetrapal, const float point[3])
{
	const int64_t p[3] =
	{
		(double)point[0] * TETRAPAL_PRECISION,
		(double)point[1] * TETRAPAL_PRECISION,
		(double)point[2] * TETRAPAL_PRECISION
	};

	if (!tetrapal) return 0;
	if (tetrapal->dimensions < 1) return 0;

	int64_t pa[3], pb[3];
	int a, null[3];

	// Use the spatial grid to get a sensible starting index
	spatialgrid_locate(tetrapal, p, &a, &null[0], &null[1], &null[2]);

	// Walk through the graph, moving towards the closest vertex
	while (1)
	{
		get_coordinates(tetrapal, a, pa);
		const double dist_a = distance_squared(p, pa);

		// Check the distance to every connected vertex and find the closest
		int closest = a;

		for (size_t i = 0; i < graph_list_size(tetrapal->graph, a); i++)
		{
			const int b = graph_get(tetrapal->graph, a, i);
			get_coordinates(tetrapal, b, pb);

			const double dist_b = distance_squared(p, pb);

			if (dist_b < dist_a)
				closest = b;
		}

		// If the closest vertex is the one we started with, we are done
		if (closest == a)
			break;

		// Otherwise update the current vertex and start again
		else a = closest;
		continue;
	}

	return a;
}

int tetrapal_number_of_dimensions(const TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
		return -1;

	return tetrapal->dimensions;
}

int tetrapal_number_of_vertices(const TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
		return 0;

	return tetrapal->size;
}

int tetrapal_number_of_segments(const TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
		return 0;

	if (tetrapal->dimensions < 1)
		return 0;

	return tetrapal->number_of_segments;
}

int tetrapal_number_of_triangles(const TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
		return 0;

	if (tetrapal->dimensions < 2)
		return 0;

	return tetrapal->number_of_triangles;
}

int tetrapal_number_of_tetrahedra(const TetrapalData* tetrapal)
{
	if (tetrapal == NULL)
		return 0;

	if (tetrapal->dimensions < 3)
		return 0;

	return tetrapal->number_of_tetrahedra;
}

int tetrapal_get_vertices(const TetrapalData* tetrapal, int* buffer)
{
	if (tetrapal == NULL)
		return EXIT_FAILURE;

	if (buffer == NULL)
		return EXIT_FAILURE;

	int count = tetrapal_number_of_vertices(tetrapal);

	for (int i = 0; i < count; i++)
	{
		buffer[i] = i;
	}

	return EXIT_SUCCESS;
}

int tetrapal_get_segments(const TetrapalData* tetrapal, int* buffer)
{
	if (tetrapal == NULL)
		return EXIT_FAILURE;

	if (tetrapal->dimensions < 1)
		return EXIT_FAILURE;

	if (buffer == NULL)
		return EXIT_FAILURE;

	int count = tetrapal_number_of_segments(tetrapal);
	HashSet* stack = hashset_new(sizeof(int[2]), count, segment_hash, segment_compare_equal);

	if (stack == NULL)
		return EXIT_FAILURE;

	for (size_t i = 0; i < tetrapal->graph->size; i++)
	{
		for (size_t j = 0; j < graph_list_size(tetrapal->graph, i); j++)
		{
			int segment[2] = {i, graph_get(tetrapal->graph, i, j)};
			hashset_insert(stack, segment);
		}
	}

	// Copy contents of stack to the output buffer
	memcpy(buffer, hashset_begin(stack), sizeof(*buffer) * 2 * hashset_size(stack));

	// Free hash set
	hashset_free(stack);

	return EXIT_SUCCESS;
}

int tetrapal_get_triangles(const TetrapalData* tetrapal, int* buffer)
{
	if (tetrapal == NULL)
		return EXIT_FAILURE;

	if (tetrapal->dimensions < 2)
		return EXIT_FAILURE;

	if (buffer == NULL)
		return EXIT_FAILURE;

	int count = tetrapal_number_of_triangles(tetrapal);
	HashSet* stack = hashset_new(sizeof(int[3]), count, triangle_hash, triangle_compare_equal);

	if (stack == NULL)
		return EXIT_FAILURE;

	if (tetrapal->dimensions == 3)
	{
		// Get finite triangles by iterating through the adjacency list
		for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
		{
			const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);

			// Ignore infinite triangles
			if (is_infinite_triangle(key->a, key->b, key->c))
				continue;

			int triangle[3] = { key->a, key->b, key->c };
			hashset_insert(stack, triangle);
		}
	}
	else if (tetrapal->dimensions == 2)
	{
		// Get finite triangles by iterating through the adjacency list
		for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
		{
			const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);
			const int* vertex = hashmap_get_value(tetrapal->adjacency, i);

			// Ignore infinite triangles
			if (is_infinite_triangle(key->a, key->b, *vertex))
				continue;

			int triangle[3] = { key->a, key->b, *vertex };
			hashset_insert(stack, triangle);
		}
	}

	// Copy contents of stack to the output buffer
	memcpy(buffer, hashset_begin(stack), sizeof(*buffer) * 3 * hashset_size(stack));

	// Free hash set
	hashset_free(stack);

	return EXIT_SUCCESS;
}

int tetrapal_get_tetrahedra(const TetrapalData* tetrapal, int* buffer)
{
	if (tetrapal == NULL)
		return EXIT_FAILURE;

	if (tetrapal->dimensions < 3)
		return EXIT_FAILURE;

	if (buffer == NULL)
		return EXIT_FAILURE;

	int count = tetrapal_number_of_tetrahedra(tetrapal);
	HashSet* stack = hashset_new(sizeof(int[4]), count, tetrahedron_hash, tetrahedron_compare_equal);

	if (stack == NULL)
		return EXIT_FAILURE;

	// Get finite triangles by iterating through the adjacency list
	for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
	{
		const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);
		const int* vertex = hashmap_get_value(tetrapal->adjacency, i);

		// Ignore infinite tetrahedra
		if (is_infinite_tetrahedron(key->a, key->b, key->c, *vertex))
			continue;

		int tetrahedron[4] = { key->a, key->b, key->c, *vertex };
		hashset_insert(stack, tetrahedron);
	}

	// Copy contents of stack to the output buffer
	memcpy(buffer, hashset_begin(stack), sizeof(*buffer) * 4 * hashset_size(stack));

	// Free hash set
	hashset_free(stack);

	return EXIT_SUCCESS;
}

static uint32_t random_int(uint32_t* state)
{
	*state = 214013u * *state + 2531011u;
    return (*state >> 16) & TETRAPAL_RANDOM_MAX;
}

static uint32_t random_range_int(uint32_t* state, const int32_t range)
{
	return random_int(state) / (TETRAPAL_RANDOM_MAX / range + 1);
}

static void swap_int(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

static void clamp_float(float* value, float min, float max)
{
	*value = *value < min ? min : (*value > max ? max : *value);
}

static void get_coordinates(const TetrapalData* tetrapal, const int i, int64_t coords[3])
{
	if (i == VERTEX_INFINITE || i == VERTEX_NULL || (size_t)i >= tetrapal->size)
	{
		coords[0] = coords[1] = coords[2] = -1;
		return;
	}

	coords[0] = tetrapal->vertices[i].x;
	coords[1] = tetrapal->vertices[i].y;
	coords[2] = tetrapal->vertices[i].z;
}

static size_t adjacency_hash(const void* ptr_key)
{
	TetrapalAdjacencyKey key = *(TetrapalAdjacencyKey*)ptr_key;
	return
		(10079 * key.a * key.a + 20047 * key.b + key.c) ^
		(10079 * key.b * key.b + 20047 * key.c + key.a) ^
		(10079 * key.c * key.c + 20047 * key.a + key.b);
}

static int adjacency_compare_equal(const void* a, const void* b)
{
	TetrapalAdjacencyKey first = *(TetrapalAdjacencyKey*)a;
	TetrapalAdjacencyKey second = *(TetrapalAdjacencyKey*)b;
	return (
		(first.a == second.a && first.b == second.b && first.c == second.c) ||
		(first.a == second.b && first.b == second.c && first.c == second.a) ||
		(first.a == second.c && first.b == second.a && first.c == second.b));
}

static size_t segment_hash(const void* ptr_segment)
{
	int a = ((int*)ptr_segment)[0];
	int b = ((int*)ptr_segment)[1];

	return (10079 * a * a + 20047 * b) ^ (10079 * b * b + 20047 * a);
}

static int segment_compare_equal(const void* ptr_a, const void* ptr_b)
{
	int a1 = ((int*)ptr_a)[0];
	int a2 = ((int*)ptr_a)[1];

	int b1 = ((int*)ptr_b)[0];
	int b2 = ((int*)ptr_b)[1];

	return (a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1);
}

static void triangle_make_canonical(int* a, int* b, int* c)
{
	if (*a > *b) swap_int(a, b);
	if (*b > *c) swap_int(b, c);
	if (*a > *b) swap_int(a, b);
}

static size_t triangle_hash(const void* ptr_triangle)
{
	int a = ((int*)ptr_triangle)[0];
	int b = ((int*)ptr_triangle)[1];
	int c = ((int*)ptr_triangle)[2];

	triangle_make_canonical(&a, &b, &c);

	return (10079 * a * a + 20047 * b + c);
}

static int triangle_compare_equal(const void* ptr_a, const void* ptr_b)
{
	int a1 = ((int*)ptr_a)[0];
	int a2 = ((int*)ptr_a)[1];
	int a3 = ((int*)ptr_a)[2];

	int b1 = ((int*)ptr_b)[0];
	int b2 = ((int*)ptr_b)[1];
	int b3 = ((int*)ptr_b)[2];

	triangle_make_canonical(&a1, &a2, &a3);
	triangle_make_canonical(&b1, &b2, &b3);

	return (a1 == b1 && a2 == b2 && a3 == b3);
}

static void tetrahedron_make_canonical(int* a, int* b, int* c, int *d)
{
	if (*a > *b) swap_int(a, b);
	if (*c > *d) swap_int(c, d);
	if (*a > *c) swap_int(a, c);
	if (*b > *d) swap_int(b, d);
	if (*b > *c) swap_int(b, c);
}

static size_t tetrahedron_hash(const void* ptr_tetrahedron)
{
	int a = ((int*)ptr_tetrahedron)[0];
	int b = ((int*)ptr_tetrahedron)[1];
	int c = ((int*)ptr_tetrahedron)[2];
	int d = ((int*)ptr_tetrahedron)[3];

	tetrahedron_make_canonical(&a, &b, &c, &d);

	return (10079 * a * a * a + 20047 * b * b + 12619 * c + d);
}

static int tetrahedron_compare_equal(const void* ptr_a, const void* ptr_b)
{
	int a1 = ((int*)ptr_a)[0];
	int a2 = ((int*)ptr_a)[1];
	int a3 = ((int*)ptr_a)[2];
	int a4 = ((int*)ptr_a)[3];

	int b1 = ((int*)ptr_b)[0];
	int b2 = ((int*)ptr_b)[1];
	int b3 = ((int*)ptr_b)[2];
	int b4 = ((int*)ptr_b)[3];

	tetrahedron_make_canonical(&a1, &a2, &a3, &a4);
	tetrahedron_make_canonical(&b1, &b2, &b3, &b4);

	return (a1 == b1 && a2 == b2 && a3 == b3 && a4 == b4);
}

// ===============================================================
//		Interpolation
// ===============================================================	

static int interpolate_3d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, int* d, double* u, double* v, double* w, double* x)
{
	int64_t pa[3], pb[3], pc[3], pd[3];
	uint32_t random_state = 1;

	spatialgrid_locate(tetrapal, p, a, b, c, d);

	//*a = tetrapal->last[0];
	//*b = tetrapal->last[1];
	//*c = tetrapal->last[2];
	//*d = tetrapal->last[3];

	// Start walking from within the triangulation
	while (1) {

		// If we arrive at an infinite tetrahedron, we must perform a convex hull walk
		if (is_infinite_tetrahedron(*a, *b, *c, *d))
			goto WALK_FACE;

		// Randomly rotate indices (stochastic walk)
		rotate_tetrahedron(a, b, c, d, random_range_int(&random_state, 4));
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, *c, pc);
		get_coordinates(tetrapal, *d, pd);

		// Test the query point against every face
		if (orient3d_fast(pa, pb, pc, p) < 0) { *d = adjacent_tetrahedron(tetrapal, *c, *b, *a); swap_int(a, c); continue; }
		if (orient3d_fast(pa, pc, pd, p) < 0) { *b = adjacent_tetrahedron(tetrapal, *d, *c, *a); swap_int(b, d); continue; }
		if (orient3d_fast(pa, pd, pb, p) < 0) { *c = adjacent_tetrahedron(tetrapal, *b, *d, *a); swap_int(c, d); continue; }
		if (orient3d_fast(pb, pd, pc, p) < 0) { *a = adjacent_tetrahedron(tetrapal, *c, *d, *b); swap_int(a, d); continue; }

		// If none of the faces return negative, then we have found the enclosing tetrahedron
		barycentric_3d(p, pa, pb, pc, pd, u, v, w, x);
		return 4;
	}

WALK_FACE:
	while (1)
	{
		// Randomly rotate indices (stochastic walk)
		rotate_triangle(a, b, c, random_range_int(&random_state, 3));
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, *c, pc);

		// Check if the point lies in an incident edge region
		if (orient2d_fast(pa, pb, p, pa, pb, pc) < 0) { goto WALK_EDGE; }
		if (orient2d_fast(pb, pc, p, pb, pc, pa) < 0) { *a = *b; *b = *c; goto WALK_EDGE; }
		if (orient2d_fast(pc, pa, p, pc, pa, pb) < 0) { *b = *a; *a = *c; goto WALK_EDGE; }

		// Point lies within this region
		barycentric_2d(p, pa, pb, pc, u, v, w);
		return 3;
	}

WALK_EDGE:
	while (1)
	{
		int j, k;
		j = adjacent_tetrahedron(tetrapal, *a, *b, VERTEX_INFINITE);
		k = adjacent_tetrahedron(tetrapal, *b, *a, VERTEX_INFINITE);

		int64_t pa[3], pb[3], pj[3], pk[3];
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, j, pj);
		get_coordinates(tetrapal, k, pk);

		// Check if point lies in an incident vertex region
		if (orient1d_fast(pa, pb, p) < 0) { goto WALK_VERTEX; }
		if (orient1d_fast(pb, pa, p) < 0) { *a = *b; goto WALK_VERTEX; }

		// Check if point lies in an incident face region
		if (orient2d_fast(pa, pb, p, pb, pa, pj) < 0) { swap_int(a, b); *c = j; goto WALK_FACE; }
		if (orient2d_fast(pb, pa, p, pa, pb, pk) < 0) { *c = k; goto WALK_FACE; }

		// Point lies within this region
		barycentric_1d(p, pa, pb, u, v);
		return 2;
	}

WALK_VERTEX:
	while (1)
	{
		// Get all the hull edges connected to this vertex and check if the orthogonal projection of [p] lies on any of these edges
		for (size_t i = 0; i < graph_list_size(tetrapal->graph, *a); i++) {

			*b = graph_get(tetrapal->graph, *a, i);

			// Only check vertices that lie on the convex hull
			if (adjacent_tetrahedron(tetrapal, *a, *b, VERTEX_INFINITE) == VERTEX_NULL)
				continue;

			int64_t pa[3], pb[3];
			get_coordinates(tetrapal, *a, pa);
			get_coordinates(tetrapal, *b, pb);

			// Check if point lies in an incident edge region
			if (orient1d_fast(pa, pb, p) > 0) { goto WALK_EDGE; }
		}

		// Point lies within this region
		*u = 1.0;
		return 1;
	}
}

static int interpolate_2d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, int* c, double* u, double* v, double* w)
{
	int64_t pa[3], pb[3], pc[3];
	uint32_t random_state = 1;

	int null;
	spatialgrid_locate(tetrapal, p, a, b, c, &null);

	//*a = tetrapal->last[0];
	//*b = tetrapal->last[1];
	//*c = tetrapal->last[2];

	// Walk as normal, but use the fast predicates (numerical stability is not that important)
	while (1)
	{
		// If we arrive at an infinite triangle, then our point lies outside the current triangulation
		if (is_infinite_triangle(*a, *b, *c))
		{
			break;
		}

		// Randomly rotate indices (stochastic walk)
		rotate_triangle(a, b, c, random_range_int(&random_state, 3));
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, *c, pc);

		// Test the query point against every face
		if (orient2d_fast(pa, pb, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *c = adjacent_triangle(tetrapal, *b, *a); swap_int(a, b); continue; }
		if (orient2d_fast(pb, pc, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *a = adjacent_triangle(tetrapal, *c, *b); swap_int(c, a); continue; }
		if (orient2d_fast(pc, pa, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *b = adjacent_triangle(tetrapal, *a, *c); swap_int(c, b); continue; }

		// If none of the faces return negative, then we have found the enclosing triangle
		barycentric_2d(p, pa, pb, pc, u, v, w);

		return 3;
	}

	// Walk the hull
	while (1)
	{
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);

		int adjacent[2] = {
			adjacent_triangle(tetrapal, *a, VERTEX_INFINITE),
			adjacent_triangle(tetrapal, VERTEX_INFINITE, *b) };

		int64_t pj[3], pk[3];
		get_coordinates(tetrapal, adjacent[0], pj);
		get_coordinates(tetrapal, adjacent[1], pk);

		if ((orient2d_fast(pj, pa, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) > 0) && (orient1d_fast(pa, pj, p) > 0)) { *b = adjacent[0]; swap_int(a, b); continue; }
		if ((orient2d_fast(pb, pk, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) > 0) && (orient1d_fast(pb, pk, p) > 0)) { *a = adjacent[1]; swap_int(a, b); continue; }

		// Check whether the point lies within a vertex region
		if (orient1d_fast(pa, pb, p) < 0) { *u = 1.0; *b = *c = 0; *v = *w = 0.0; return 1; }
		else if (orient1d_fast(pb, pa, p) < 0) { *u = 1.0; *a = *b; *b = *c = 0; *v = *w = 0.0; return 1; }

		// Get the weights for this segment
		barycentric_1d(p, pa, pb, u, v); *c = 0; *w = 0.0;

		return 2;
	}
}

static int interpolate_1d(const TetrapalData* tetrapal, const int64_t p[3], int* a, int* b, double* u, double* v)
{
	int64_t pa[3], pb[3];

	int null[2];
	spatialgrid_locate(tetrapal, p, a, b, &null[0], &null[1]);

	//*a = tetrapal->last[0];
	//*b = tetrapal->last[1];

	while (1)
	{
		// If we arrive at an infinite segment, then our point lies beyond the boundary vertex
		if (is_infinite_segment(*a, *b))
		{
			if (*a == VERTEX_INFINITE)
			{
				*a = *b;
			}

			*b = 0;
			*u = 1.0;
			*v = 0.0;

			return 1;
		}

		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);

		// Test the query point against both sides of the line
		if (orient1d_fast(pa, pb, p) < 0) { *b = adjacent_segment(tetrapal, *a, VERTEX_PREV); swap_int(a, b); continue; }
		if (orient1d_fast(pb, pa, p) < 0) { *a = adjacent_segment(tetrapal, *b, VERTEX_NEXT); swap_int(a, b); continue; }

		// If none of the faces return negative, then we have found the enclosing segment
		barycentric_1d(p, pa, pb, u, v);

		return 2;
	}
}

// ===============================================================
//		3D Triangulation
// ===============================================================	

static int triangulate_3d(TetrapalData* tetrapal)
{
	if (tetrapal->size < 4)
		return 0;

	// Determine the first valid tetrahedron
	int a = 0, b = 1, c = 2, d = 3;
	int64_t pa[3], pb[3], pc[3], pd[3];
	get_coordinates(tetrapal, a, pa);
	get_coordinates(tetrapal, b, pb);
	get_coordinates(tetrapal, c, pc);
	get_coordinates(tetrapal, d, pd);

	while (orient1d_filtered(pa, pb, pb) == 0)
		if ((size_t)b < tetrapal->size - 1) get_coordinates(tetrapal, ++b, pb);
		else return 0;

	while (orient2d_filtered(pa, pb, pc, pa, pb, pc) == 0)
		if ((size_t)c < tetrapal->size - 1) get_coordinates(tetrapal, ++c, pc);
		else return 0;

	while (orient3d_filtered(pa, pb, pc, pd) == 0)
		if ((size_t)d < tetrapal->size - 1) get_coordinates(tetrapal, ++d, pd);
		else return 0;

	// Ensure positive orientation
	if (orient3d_filtered(pa, pb, pc, pd) < 0)
		swap_int(&a, &c);

	add_tetrahedron(tetrapal, a, b, c, d);
	add_tetrahedron(tetrapal, a, c, b, VERTEX_INFINITE);
	add_tetrahedron(tetrapal, a, b, d, VERTEX_INFINITE);
	add_tetrahedron(tetrapal, a, d, c, VERTEX_INFINITE);
	add_tetrahedron(tetrapal, b, c, d, VERTEX_INFINITE);
	tetrapal->dimensions = 3;

#ifdef TETRAPAL_DEBUG
	clock_t time_start;
	TIME_LOCATE = TIME_STELLATE = 0;
#endif

	// Start incremental insertion
	for (size_t i = 0; i < tetrapal->size; i++)
	{
		a = tetrapal->last[0];
		b = tetrapal->last[1];
		c = tetrapal->last[2];
		d = tetrapal->last[3];

		int64_t pi[3];
		get_coordinates(tetrapal, i, pi);

	#ifdef TETRAPAL_DEBUG
		time_start = clock();
	#endif

		if (!locate_tetrahedron(tetrapal, pi, &a, &b, &c, &d))
			continue;

	#ifdef TETRAPAL_DEBUG
		TIME_LOCATE += clock() - time_start;
		time_start = clock();
	#endif

		delete_tetrahedron(tetrapal, a, b, c, d);
		consider_tetrahedron(tetrapal, a, b, c, i);
		consider_tetrahedron(tetrapal, a, c, d, i);
		consider_tetrahedron(tetrapal, a, d, b, i);
		consider_tetrahedron(tetrapal, b, d, c, i);

	#ifdef TETRAPAL_DEBUG
		TIME_STELLATE += clock() - time_start;
	#endif
	}

#ifdef TETRAPAL_DEBUG
	printf("LOCATE TIME: %ums\n", (unsigned int)TIME_LOCATE);
	printf("STELLATE TIME: %ums\n", (unsigned int)TIME_STELLATE);
#endif

	// Build the vertex graph by iterating through the adjacency map
	for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
	{
		const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);
		const int* vertex = hashmap_get_value(tetrapal->adjacency, i);

		graph_insert(tetrapal->graph, key->a, *vertex);
		graph_insert(tetrapal->graph, key->b, *vertex);
		graph_insert(tetrapal->graph, key->c, *vertex);

		if (is_infinite_tetrahedron(key->a, key->b, key->c, *vertex) == 0)
		{
			tetrapal->number_of_tetrahedra++;
		}

		if (is_infinite_triangle(key->a, key->b, key->c) == 0)
		{
			tetrapal->number_of_triangles++;
		}
	}

	// Every tetrahedron is recorded four times, and each triangle twice
	tetrapal->number_of_tetrahedra /= 4;
	tetrapal->number_of_triangles /= 2;

	// Build spatial grid
	spatialgrid_build(tetrapal);

	return 1;
}

static int is_infinite_tetrahedron(int a, int b, int c, int d)
{
	return (a == VERTEX_INFINITE ||
		b == VERTEX_INFINITE ||
		c == VERTEX_INFINITE ||
		d == VERTEX_INFINITE) ? 1 : 0;
}

static void rotate_tetrahedron(int* a, int* b, int* c, int* d, int rotations)
{
	int an, bn, cn, dn;

	// CCW rotation
	switch (rotations % 4) {
	case 1: an = *b; bn = *d; cn = *c; dn = *a; break;
	case 2: an = *c; bn = *d; cn = *a; dn = *b; break;
	case 3: an = *a; bn = *d; cn = *b; dn = *c; break;
	default: return;
	}

	*a = an; *b = bn; *c = cn; *d = dn;
}

static void add_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d)
{
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, d, c }, & a);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, c, d }, & b);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, d, b }, & c);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, c }, & d);

	// If this tetrahedron is infinite, find the nearest finite one
	if (a == VERTEX_INFINITE) { a = adjacent_tetrahedron(tetrapal, b, c, d); swap_int(&c, &d); }
	else if (b == VERTEX_INFINITE) { b = adjacent_tetrahedron(tetrapal, c, a, d); swap_int(&a, &d); }
	else if (c == VERTEX_INFINITE) { c = adjacent_tetrahedron(tetrapal, d, a, b); swap_int(&a, &b); }
	else if (d == VERTEX_INFINITE) { d = adjacent_tetrahedron(tetrapal, a, c, b); swap_int(&c, &b); }

	tetrapal->last[0] = a;
	tetrapal->last[1] = b;
	tetrapal->last[2] = c;
	tetrapal->last[3] = d;
}

static void delete_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d)
{
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, d, c });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, c, d });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, d, b });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, c });
}

static int adjacent_tetrahedron(const TetrapalData* tetrapal, int a, int b, int c)
{
	const int* adjacent = hashmap_find(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, c });

	if (!adjacent)
		return VERTEX_NULL;
	else
		return *adjacent;
}

static int check_conflict_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int d, int e)
{
	int64_t pa[3], pb[3], pc[3], pd[3], pe[3];
	get_coordinates(tetrapal, a, pa);
	get_coordinates(tetrapal, b, pb);
	get_coordinates(tetrapal, c, pc);
	get_coordinates(tetrapal, d, pd);
	get_coordinates(tetrapal, e, pe);

	if (a == VERTEX_INFINITE) {
		int orient = orient3d_filtered(pe, pb, pc, pd);
		return orient ? orient : incircle_filtered(pb, pc, pd, pe);
	}
	if (b == VERTEX_INFINITE) {
		int orient = orient3d_filtered(pa, pe, pc, pd);
		return orient ? orient : incircle_filtered(pd, pc, pa, pe);
	}
	if (c == VERTEX_INFINITE) {
		int orient = orient3d_filtered(pa, pb, pe, pd);
		return orient ? orient : incircle_filtered(pa, pb, pd, pe);
	}
	if (d == VERTEX_INFINITE) {
		int orient = orient3d_filtered(pa, pb, pc, pe);
		return orient ? orient : incircle_filtered(pa, pb, pc, pe);
	}

	return insphere_filtered(pa, pb, pc, pd, pe);
}

static void consider_tetrahedron(TetrapalData* tetrapal, int a, int b, int c, int e)
{
	// Find tetrahedron cbad opposite facet abc from u
	int d = adjacent_tetrahedron(tetrapal, c, b, a);

	if (d == VERTEX_NULL)
		return;

	if (check_conflict_tetrahedron(tetrapal, c, b, a, d, e) > 0) {
		delete_tetrahedron(tetrapal, c, b, a, d);
		consider_tetrahedron(tetrapal, a, b, d, e);
		consider_tetrahedron(tetrapal, d, b, c, e);
		consider_tetrahedron(tetrapal, a, d, c, e);
	}
	else
		add_tetrahedron(tetrapal, a, b, c, e);
}

static int locate_tetrahedron(TetrapalData* tetrapal, int64_t p[3], int* a, int* b, int* c, int* d)
{
	int64_t pa[3], pb[3], pc[3], pd[3];
	uint32_t random_state = 1;

	while (1)
	{
		// If we arrive at an infinite tetrahedron, then our point lies outside the current triangulation
		if (is_infinite_tetrahedron(*a, *b, *c, *d))
			return 1;

		// Randomly rotate indices (stochastic walk)
		rotate_tetrahedron(a, b, c, d, random_range_int(&random_state, 4));
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, *c, pc);
		get_coordinates(tetrapal, *d, pd);

		// Reject coincident points
		if (is_coincident(p, pa) ||
			is_coincident(p, pb) ||
			is_coincident(p, pc) ||
			is_coincident(p, pd))
			return 0;

		// Test the query point against every face
		if (orient3d_filtered(pa, pb, pc, p) < 0) { *d = adjacent_tetrahedron(tetrapal, *c, *b, *a); swap_int(a, c); continue; }
		if (orient3d_filtered(pa, pc, pd, p) < 0) { *b = adjacent_tetrahedron(tetrapal, *d, *c, *a); swap_int(b, d); continue; }
		if (orient3d_filtered(pa, pd, pb, p) < 0) { *c = adjacent_tetrahedron(tetrapal, *b, *d, *a); swap_int(c, d); continue; }
		if (orient3d_filtered(pb, pd, pc, p) < 0) { *a = adjacent_tetrahedron(tetrapal, *c, *d, *b); swap_int(a, d); continue; }

		// If none of the faces return negative, then we have found the enclosing tetrahedron
		return 1;
	}
}

// ===============================================================
//		2D Triangulation
// ===============================================================	

static int triangulate_2d(TetrapalData* tetrapal)
{
	if (tetrapal->size < 3)
		return 0;

	// Determine the first valid triangle
	int a = 0, b = 1, c = 2;
	int64_t pa[3], pb[3], pc[3];
	get_coordinates(tetrapal, a, pa);
	get_coordinates(tetrapal, b, pb);
	get_coordinates(tetrapal, c, pc);

	while (orient1d_filtered(pa, pb, pb) == 0)
		if ((size_t)b < tetrapal->size - 1) get_coordinates(tetrapal, ++b, pb);
		else return 0;

	while (orient2d_filtered(pa, pb, pc, pa, pb, pc) == 0)
		if ((size_t)c < tetrapal->size - 1) get_coordinates(tetrapal, ++c, pc);
		else return 0;

	// Record the orientation defined by the first triangle
	memcpy(tetrapal->orient_a, pa, sizeof(*tetrapal->orient_a) * 3);
	memcpy(tetrapal->orient_b, pb, sizeof(*tetrapal->orient_b) * 3);
	memcpy(tetrapal->orient_c, pc, sizeof(*tetrapal->orient_c) * 3);

	add_triangle(tetrapal, a, b, c);
	add_triangle(tetrapal, b, a, VERTEX_INFINITE);
	add_triangle(tetrapal, c, b, VERTEX_INFINITE);
	add_triangle(tetrapal, a, c, VERTEX_INFINITE);
	tetrapal->dimensions = 2;

	// Start incremental insertion
	for (size_t i = 0; i < tetrapal->size; i++)
	{
		a = tetrapal->last[0];
		b = tetrapal->last[1];
		c = tetrapal->last[2];

		int64_t pi[3];
		get_coordinates(tetrapal, i, pi);

		if (!locate_triangle(tetrapal, pi, &a, &b, &c))
			continue;

		delete_triangle(tetrapal, a, b, c);
		consider_triangle(tetrapal, a, b, i);
		consider_triangle(tetrapal, b, c, i);
		consider_triangle(tetrapal, c, a, i);
	}

	// Build the vertex graph by iterating through the adjacency map
	for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
	{
		const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);
		const int* vertex = hashmap_get_value(tetrapal->adjacency, i);

		graph_insert(tetrapal->graph, key->a, *vertex);
		graph_insert(tetrapal->graph, key->b, *vertex);

		if (is_infinite_triangle(key->a, key->b, *vertex) == 0)
		{
			tetrapal->number_of_triangles++;
		}
	}

	// Each triangle is recorded three times
	tetrapal->number_of_triangles /= 3;

	// Build spatial grid
	spatialgrid_build(tetrapal);

	return 1;
}

static int locate_triangle(TetrapalData* tetrapal, int64_t p[3], int* a, int* b, int* c)
{
	int64_t pa[3], pb[3], pc[3];
	uint32_t random_state = 1;

	while (1)
	{
		// If we arrive at an infinite triangle, then our point lies outside the current triangulation
		if (is_infinite_triangle(*a, *b, *c))
			return 1;

		// Randomly rotate indices (stochastic walk)
		rotate_triangle(a, b, c, random_range_int(&random_state, 3));
		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);
		get_coordinates(tetrapal, *c, pc);

		// Reject coincident points
		if (is_coincident(p, pa) ||
			is_coincident(p, pb) ||
			is_coincident(p, pc))
			return 0;

		// Test the query point against every face
		if (orient2d_filtered(pa, pb, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *c = adjacent_triangle(tetrapal, *b, *a); swap_int(a, b); continue; }
		if (orient2d_filtered(pb, pc, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *a = adjacent_triangle(tetrapal, *c, *b); swap_int(b, c); continue; }
		if (orient2d_filtered(pc, pa, p, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c) < 0) { *b = adjacent_triangle(tetrapal, *a, *c); swap_int(c, a); continue; }

		// If none of the faces return negative, then we have found the enclosing triangle
		return 1;
	}
}

static void consider_triangle(TetrapalData* tetrapal, int a, int b, int e)
{
	int c = adjacent_triangle(tetrapal, b, a);

	if (c == VERTEX_NULL)
		return;

	if (check_conflict_triangle(tetrapal, b, a, c, e) > 0) {
		delete_triangle(tetrapal, b, a, c);
		consider_triangle(tetrapal, a, c, e);
		consider_triangle(tetrapal, c, b, e);
	}
	else
		add_triangle(tetrapal, a, b, e);
}

static int check_conflict_triangle(TetrapalData* tetrapal, int a, int b, int c, int e)
{
	int64_t pa[3], pb[3], pc[3], pe[3];
	get_coordinates(tetrapal, a, pa);
	get_coordinates(tetrapal, b, pb);
	get_coordinates(tetrapal, c, pc);
	get_coordinates(tetrapal, e, pe);

	if (a == VERTEX_INFINITE) {
		int orient = orient2d_filtered(pb, pc, pe, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c);
		return orient ? orient : insegment_filtered(pb, pc, pe);
	}
	if (b == VERTEX_INFINITE) {
		int orient = orient2d_filtered(pc, pa, pe, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c);
		return orient ? orient : insegment_filtered(pc, pa, pe);
	}
	if (c == VERTEX_INFINITE) {
		int orient = orient2d_filtered(pa, pb, pe, tetrapal->orient_a, tetrapal->orient_b, tetrapal->orient_c);
		return orient ? orient : insegment_filtered(pa, pb, pe);
	}

	return incircle_filtered(pa, pb, pc, pe);
}

static void add_triangle(TetrapalData* tetrapal, int a, int b, int c)
{
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, VERTEX_NULL }, & c);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, c, VERTEX_NULL }, & a);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { c, a, VERTEX_NULL }, & b);

	// If this triangle is infinite, store the nearest finite one
	if (a == VERTEX_INFINITE) { a = adjacent_triangle(tetrapal, c, b); swap_int(&c, &b); }
	else if (b == VERTEX_INFINITE) { b = adjacent_triangle(tetrapal, a, c); swap_int(&a, &c); }
	else if (c == VERTEX_INFINITE) { c = adjacent_triangle(tetrapal, b, a); swap_int(&b, &a); }
	tetrapal->last[0] = a;
	tetrapal->last[1] = b;
	tetrapal->last[2] = c;
}

static void delete_triangle(TetrapalData* tetrapal, int a, int b, int c)
{
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, VERTEX_NULL });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, c, VERTEX_NULL });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { c, a, VERTEX_NULL });
}

static int adjacent_triangle(const TetrapalData* tetrapal, int a, int b)
{
	const int* adjacent = hashmap_find(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, b, VERTEX_NULL });

	if (!adjacent)
		return VERTEX_NULL;
	else
		return *adjacent;
}

static void rotate_triangle(int* a, int* b, int* c, int rotations)
{
	int an, bn, cn;

	// CCW rotation
	switch (rotations % 3) {
	case 1: an = *b; bn = *c; cn = *a; break;
	case 2: an = *c; bn = *a; cn = *b; break;
	default: return;
	}

	*a = an; *b = bn; *c = cn;
}

static int is_infinite_triangle(int a, int b, int c)
{
	return (a == VERTEX_INFINITE || b == VERTEX_INFINITE || c == VERTEX_INFINITE) ? 1 : 0;
}

// ===============================================================
//		1D Triangulation
// ===============================================================	

static int triangulate_1d(TetrapalData* tetrapal)
{
	if (tetrapal->size < 2)
		return 0;

	int a = 0, b = 1;
	int64_t pa[3], pb[3];
	get_coordinates(tetrapal, a, pa);
	get_coordinates(tetrapal, b, pb);

	while (orient1d_filtered(pa, pb, pb) == 0)
		if ((size_t)b < tetrapal->size - 1) get_coordinates(tetrapal, ++b, pb);
		else return 0;

	// Ensure positive orientation
	if (orient1d_exact(pa, pb, pb) < 0)
		swap_int(&a, &b);

	add_segment(tetrapal, a, b);
	add_segment(tetrapal, VERTEX_INFINITE, a);
	add_segment(tetrapal, b, VERTEX_INFINITE);
	tetrapal->dimensions = 1;

	// Start incremental insertion
	for (size_t i = 0; i < tetrapal->size; i++)
	{
		a = tetrapal->last[0];
		b = tetrapal->last[1];

		int64_t pi[3];
		get_coordinates(tetrapal, i, pi);

		if (!locate_segment(tetrapal, pi, &a, &b))
			continue;

		delete_segment(tetrapal, a, b);
		add_segment(tetrapal, a, i);
		add_segment(tetrapal, i, b);
	}

	// Build the vertex graph by iterating through the adjacency map
	for (size_t i = 0; i < hashmap_size(tetrapal->adjacency); i++)
	{
		const TetrapalAdjacencyKey* key = hashmap_get_key(tetrapal->adjacency, i);
		const int* vertex = hashmap_get_value(tetrapal->adjacency, i);

		graph_insert(tetrapal->graph, key->a, *vertex);
	}

	// Build spatial grid
	spatialgrid_build(tetrapal);

	return 1;
}

static int locate_segment(TetrapalData* tetrapal, int64_t p[3], int* a, int* b)
{
	int64_t pa[3], pb[3];

	while (1)
	{
		// If we arrive at an infinite segment, then our point lies beyond the boundary vertex
		if (is_infinite_segment(*a, *b))
			return 1;

		get_coordinates(tetrapal, *a, pa);
		get_coordinates(tetrapal, *b, pb);

		// Reject coincident points
		if (is_coincident(p, pa) ||
			is_coincident(p, pb))
			return 0;

		// Test the query point against both sides of the line
		if (orient1d_filtered(pa, pb, p) < 0) { *b = adjacent_segment(tetrapal, *a, VERTEX_PREV); swap_int(a, b); continue; }
		if (orient1d_filtered(pb, pa, p) < 0) { *a = adjacent_segment(tetrapal, *b, VERTEX_NEXT); swap_int(a, b); continue; }

		// If none of the faces return negative, then we have found the enclosing segment
		return 1;
	}
}

static void add_segment(TetrapalData* tetrapal, int a, int b)
{
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, VERTEX_NEXT, VERTEX_NULL }, & b);
	hashmap_insert(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, VERTEX_PREV, VERTEX_NULL }, & a);

	// If this segment is infinite, store the nearest finite one
	if (a == VERTEX_INFINITE) { a = adjacent_segment(tetrapal, b, VERTEX_NEXT); swap_int(&a, &b); }
	else if (b == VERTEX_INFINITE) { b = adjacent_segment(tetrapal, a, VERTEX_PREV); swap_int(&a, &b); }
	tetrapal->last[0] = a;
	tetrapal->last[1] = b;
}

static void delete_segment(TetrapalData* tetrapal, int a, int b)
{
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, VERTEX_NEXT, VERTEX_NULL });
	hashmap_remove(tetrapal->adjacency, &(TetrapalAdjacencyKey) { b, VERTEX_PREV, VERTEX_NULL });
}

static int adjacent_segment(const TetrapalData* tetrapal, int a, TetrapalEnum direction)
{
	// Direction must be valid
	if (direction != VERTEX_NEXT && direction != VERTEX_PREV)
		return VERTEX_NULL;

	const int* adjacent = hashmap_find(tetrapal->adjacency, &(TetrapalAdjacencyKey) { a, direction, VERTEX_NULL });

	if (!adjacent)
		return VERTEX_NULL;
	else
		return *adjacent;
}

static int is_infinite_segment(int a, int b)
{
	return (a == VERTEX_INFINITE || b == VERTEX_INFINITE) ? 1 : 0;
}

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
