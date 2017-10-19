const EXTRA_LIBS = `
declare namespace bluemath {
  declare namespace common {
declare class AABB {
    private _min;
    private _max;
    constructor(arg0: number | number[] | TypedArray, arg1?: number[] | TypedArray);
    readonly min: NDArray;
    readonly max: NDArray;
    /**
     * Update this AABB to include given coordinate
     */
    update(coord: number[] | NDArray): void;
    merge(other: AABB): void;
}

declare class Complex {
    real: number;
    imag: number;
    constructor(real?: number, imag?: number);
    clone(): Complex;
    inverse(): Complex;
    isEqual(other: Complex, tolerance?: number): boolean;
    toString(precision?: number): string;
}

declare const EPSILON = 0.000001;

declare type NumberType = 'i8' | 'ui8' | 'i16' | 'ui16' | 'i32' | 'ui32' | 'f32' | 'f64';
declare type TypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Float32Array | Float64Array;

    shape?: number[];
    datatype?: NumberType;
    fill?: number;
    idata?: number[];
}
/**
 * N-Dimensional Array
 * ===
 *
 * It can store real as well as complex numbers in n-dimensions
 * It can be used to store Vectors (1D) or Matrices (2D).
 * This class stores the data internally in flat typed arrays
 *
 * NDArray is the central class of Bluemath library.
 * It's used to input and output data to/from most of the APIs of this library.
 *
 * Construction
 * ---
 *
 * You can create an NDArray
 *
 * * With shape and/or data type
 * // 3-dimensional array with 32-bit integer storage
 * new NDArray({shape:[3,4,3],datatype:'i32'});
 *
 * * Initializing it with array data
 * // 2x3 Matrix with 64-bit floating point (double) storage
 * new NDArray([[1,1,1],[4,4,4]],{datatype:'f64'});
 *
 * * Using standard functions
 * zeros([2,2,2]); // Returns 2x2x2 NDArray of zeros
 * eye([4,4]); // Creates 4x4 Identity matrix
 *
 * Basic math operations
 * ---
 *
 * Bluemath provides functions that allow basic math operations
 * on NDArrays
 *
 * [[add]]
 *
 * [[sub]]
 *
 * [[mul]]
 *
 * [[div]]
 */
declare class NDArray {
    /**
     * Array of array dimensions. First being the outermost dimension.
     */
    private _shape;
    /**
     * Size of the data (i.e. number of real/complex numbers stored
     * in this array)
     */
    private _size;
    /**
     * Data type of each number, specified by a string code
     */
    private _datatype;
    /**
     * Real part of number elements is stored in this array
     */
    private _data;
    /**
     * If any number element of this array is Complex then its
     * imaginary part is stored in _idata sparse array object
     * indexed against its address.
     * Note that _idata is not a TypedArray as _data. This way
     * the storage is optimized for the use cases where real number
     * data is common, but in some fringe cases the number could be
     * complex.
     */
    private _idata;
    constructor(arg0: TypedArray | Array<any> | NDArrayOptions, arg1?: NDArrayOptions);
    readonly shape: number[];
    readonly size: number;
    is1D(): boolean;
    is2D(): boolean;
    /**
     * Number of elements in outermost (i.e. 0th) dimension
     */
    readonly length: number;
    readonly data: TypedArray;
    readonly datatype: NumberType;
    /**
     * Set new shape for the data stored in the array
     * The old data remains intact. If the total size with the new shape
     * is larger than the old size, then excess elements of the data are
     * fill with zero.
     * @param shape New shape
     */
    reshape(shape: number[]): this;
    /**
     * Create deep copy of the array
     */
    clone(): NDArray;
    private _calcSize();
    private _alloc(size, data?, datatype?);
    _indexToAddress(...indices: number[]): number;
    /**
     * @hidden
     */
    private static mapAddressToIndex(addr, shape);
    /**
     * @hidden
     */
    _addressToIndex(addr: number): any[];
    /**
     * Create nested array
     */
    toArray(): any;
    /**
     * Set all members of this array to given value
     */
    fill(value: number): void;
    private isSliceIndex(index);
    /**
     * Set member at given index
     * All but the last argument should specify the index.
     * The last argument is the value to set.
     */
    set(...args: (number | Complex | string | undefined | null | NDArray)[]): void;
    /**
     * Swaps matrix rows (this must be a 2D array)
     */
    swaprows(i: number, j: number): void;
    /**
     * @hidden
     */
    datacompare(otherdata: TypedArray, otheridata: number[], tolerance?: number): boolean;
    /**
     * Iterate over each element, invoke a callback with each index and value
     */
    forEach(callback: (value: number | Complex, ...index: number[]) => void): void;
    /**
     * @hidden
     */
    private static areShapesEqual(shape1, shape2);
    /**
     * Checks if the shape of this ndarray matches the shape of other
     */
    isShapeEqual(other: NDArray): boolean;
    /**
     * Does equality test for each element of the array as well as the
     * shape of the arrays
     * @param other Other NDArray to compare with
     * @param tolerance
     */
    isEqual(other: NDArray, tolerance?: number): boolean;
    /**
     * Return 1D copy of this array
     */
    flatten(): NDArray;
    /**
     * Change between Row-major and Column-major layout
     */
    swapOrder(): void;
    private createSliceRecipe(slices);
    private computeSliceShapeAndSize(slice_recipe);
    /**
     * Shorthand for get(...) method to avoid casting to <number>
     */
    getN(...slices: (string | number | undefined | null)[]): number;
    /**
     * Shorthand for get(...) method to avoid casting to <NDArray>
     */
    getA(...slices: (string | number | undefined | null)[]): NDArray;
    /**
     * Shorthand for get(...) method to avoid casting to <Complex>
     */
    getC(...slices: (string | number | undefined | null)[]): Complex;
    /**
     * Returns a specific element or a new NDArray that's a subset of
     * this array as defined by the slicing recipe.
     * Each element of the slicing recipe (i.e. any argument) can be
     * * A number specifying a specific element or slice of the array
     * in given dimension.
     * * A string of the form '<start>:<stop>', specifying the range of
     * slices in the given dimension. Both '<start>' and '<stop>' are
     * optional
     *
     * Caveats
     * ---
     * * Negative indices not supported yet
     */
    get(...slices: (string | number | undefined | null)[]): NDArray | number | Complex;
    /**
     * @hidden
     */
    take(indices: number[], axis: number): NDArray;
    /**
     * @hidden
     */
    max(axis?: number | number[]): number | NDArray;
    /**
     * @hidden
     */
    min(): void;
    /**
     * @hidden
     */
    mean(): void;
    /**
     * @hidden
     */
    all(): void;
    /**
     * @hidden
     */
    any(): void;
    /**
     * @hidden
     */
    sort(): void;
    /**
     * @hidden
     */
    argsort(): void;
    copyfrom(other: NDArray): void;
    copyto(other: NDArray): void;
    toString(precision?: number): any;
    toHTML(precision?: number): any;
}
declare class Vec2 extends NDArray {
    constructor(x: number, y: number);
}
declare class Vec3 extends NDArray {
    constructor(x: number, y: number, z: number);
}

/**
 * Convert angle to degrees
 */
declare function todeg(angleInRadians: number): number;
/**
 * Convert angle to radians
 */
declare function torad(angleInDegrees: number): number;
/**
 * Check if input equals zero within given tolerance
 */
declare function iszero(x: number, tolerance?: number): boolean;
/**
 * Check if two input numbers are equal within given tolerance
 */
declare function isequal(a: number, b: number, tolerance?: number): boolean;
/**
 * Find cube root of given number. Math.pow return NaN while taking
 * cube root of negative number, because some of the results might
 * be complex numbers. This function only return the real cubeRoot
 * of given number
 */
declare function cuberoot(x: number): number;
/**
 * Generate array of integers within given range.
 * If both a and b are specified then return [a,b)
 * if only a is specifed then return [0,a)
 */
declare function range(a: number, b?: number): NDArray;
/**
 * Creates m-by-n Identity matrix
 *
 * eye(2) // Creates 2x2 Identity matrix
 * eye([2,2]) // Creates 2x2 Identity matrix
 * eye([2,3]) // Create 2x3 Identity matrix with main diagonal set to 1
 * eye(2,'i32') // Creates 2x2 Identity matrix of 32-bit integers
 */
declare function eye(arg0: number | number[], datatype?: NumberType): NDArray;
declare function count(arr: NDArray, item: number, tolerance?: number): number;
/**
 * Creates NDArray filled with zeros
 *
 * zeros(2) // Creates array of zeros of length 2
 * zeros([2,2,2]) // Create 2x2x2 matrix of zeros
 * zeros(2,'i16') // Creates array of 2 16-bit integers filled with zeros
 */
declare function zeros(arg0: number | number[], datatype?: NumberType): NDArray;
/**
 * Creates empty NDArray of given shape or of given length if argument is
 * a number
 */
declare function empty(arg0: number | number[], datatype?: NumberType): NDArray;
/**
 * Shorthand method to create new NDArray object from Javascript Array
 */
declare function arr(arg: any[]): NDArray;
/**
 * Compute dot product of A and B, where both of them are 1D vectors of
 * same length
 */
declare function dot(A: NDArray, B: NDArray): number;
/**
 * Computes cross product of A and B
 * Only defined for A and B to 1D vectors of length at least 3
 * Only first 3 elements of A and B are used
 */
declare function cross(A: NDArray, B: NDArray): NDArray;
/**
 * Computes length or magnitude of A, where A is a 1D vector
 */
declare function length(A: NDArray): number;
/**
 * Computes direction vector of A, where A is a 1D vector
 */
declare function dir(A: NDArray): NDArray;
/**
 * Add all arguments in accordance to their types
 * The arguments could be NDArray or numbers (real/complex).
 * If some of them are NDArray's, then their shapes have to match,
 * otherwise exception is thrown
 * The order of addition starts from left to right
 */
declare function add(...args: (NDArray | number | Complex)[]): number | Complex | NDArray;
/**
 * Multiply all arguments in accordance with their data types
 * Each argument can be a number (real or complex) or NDArray.
 * If some of the arguments are NDArrays, then their shapes should
 * be compatible with the other operand of multiplication operation,
 * otherwise an exception is thrown
 * The order of multiplication starts from left to right
 */
declare function mul(...args: (NDArray | number | Complex)[]): number | Complex | NDArray;
/**
 * Subtract second argument from first
 * The arguments could be a number (real or complex) or NDArray.
 * If some of the arguments are NDArrays, then their shapes should
 * be compatible with the other operand of subtraction operation,
 * otherwise an exception is thrown
 */
declare function sub(a: number | Complex | NDArray, b: number | Complex | NDArray): number | Complex | NDArray;
/**
 * Divide first argument by second
 * The first argument can be a number (real or complex) or NDArray.
 * The second argument can be a number (real or complex)
 */
declare function div(a: number | Complex | NDArray, b: number | Complex): number | Complex | NDArray;

  }
  declare namespace linalg {

/**
 * @hidden
 */
declare function asum(X: TypedArray): number;
/**
 * @hidden
 */
declare function axpy(X: TypedArray, a: number, Y: TypedArray): void;


/**
 * @hidden
 */
declare function dot(vx: TypedArray, vy: TypedArray): any;

/**
 * @hidden
 */
declare function gemv(alpha: number, mA: TypedArray, m: number, n: number, vx: TypedArray, vy: TypedArray, beta: number): void;

/**
 * @hidden
 */
declare function gemm(mA: TypedArray, mB: TypedArray, mC: TypedArray, m: number, n: number, k: number, alpha: number, beta: number): void;

/**
 * @hidden
 */
declare const SIZE_CHAR = 1;
/**
 * @hidden
 */
declare const SIZE_INT = 4;
/**
 * @hidden
 */
declare const SIZE_DOUBLE = 8;
/**
 * @hidden
 */
declare const SIZE_SINGLE = 4;
/**
 * @hidden
 */
declare const spotrf_wrap: any;
/**
 * @hidden
 */
declare const dpotrf_wrap: any;
/**
 * @hidden
 */
declare const sgesv_wrap: any;
/**
 * @hidden
 */
declare const dgesv_wrap: any;
/**
 * @hidden
 */
declare const sgemm_wrap: any;
/**
 * @hidden
 */
declare const dgemm_wrap: any;
/**
 * @hidden
 */
declare const dgemv_wrap: any;
/**
 * @hidden
 */
declare const sgemv_wrap: any;
/**
 * @hidden
 */
declare const sdot_wrap: any;
/**
 * @hidden
 */
declare const ddot_wrap: any;
/**
 * @hidden
 */
declare const dgesdd_wrap: any;
/**
 * @hidden
 */
declare const sgesdd_wrap: any;
/**
 * @hidden
 */
declare const sgeqrf_wrap: any;
/**
 * @hidden
 */
declare const dgeqrf_wrap: any;
/**
 * @hidden
 */
declare const sorgqr_wrap: any;
/**
 * @hidden
 */
declare const dorgqr_wrap: any;
/**
 * @hidden
 */
declare const dgelsd_wrap: any;
/**
 * @hidden
 */
declare const sgetrf_wrap: any;
/**
 * @hidden
 */
declare const dgetrf_wrap: any;
/**
 * @hidden
 */
declare const sgeev_wrap: any;
/**
 * @hidden
 */
declare const dgeev_wrap: any;
/**
 * @hidden
 */
declare function defineEmVariable(type: 'i8' | 'i32' | 'f32' | 'f64', init?: number): number;
/**
 * @hidden
 */
declare function defineEmArrayVariable(type: 'i8' | 'i32' | 'f32' | 'f64', len: number, init?: TypedArray): [number, TypedArray];

/**
 * @hidden
 */
declare function geev(A: TypedArray, n: number, compleft: boolean, compright: boolean): Int8Array[];

/**
 * @hidden
 */
declare function gelsd(mA: TypedArray, m: number, n: number, nrhs: number, rcond: number, mB: TypedArray, mS: TypedArray): number;

/**
 * @hidden
 */
declare function geqrf(mA: TypedArray, m: number, n: number, mTau: TypedArray): void;

/**
 * @hidden
 */
declare function gesdd(mA: TypedArray, m: number, n: number, mU: TypedArray, mS: TypedArray, mVT: TypedArray, job: 'A' | 'N' | 'S'): void;

/**
 * @hidden
 */
declare function gesv(mA: TypedArray, mB: TypedArray, n: number, nrhs: number): TypedArray;

/**
 * @hidden
 */
declare function getrf(mA: TypedArray, m: number, n: number, mipiv: TypedArray): void;


/**
 * @hidden
 */
declare function orgqr(mA: TypedArray, m: number, n: number, k: number, mtau: TypedArray): void;

/**
 * @hidden
 */
declare function potrf(mA: TypedArray, n: number): void;

/**
 * Matrix multiplication
 *
 * At least one of the arguments has to be 2D matrix (i.e. shape mxn).
 * The other argument could be a 1D vector. It will be implicitly used
 * as 1xn matrix
 */
declare function matmul(A: NDArray, B: NDArray): NDArray;
/**
 * Computes p-norm of given Matrix or Vector
 *
 *
 * $$ \left\Vert A \right\Vert = \max_{0 \leq i < n}  \lvert a_i \rvert, p = \infty  $$
 *
 * $$ \left\Vert A \right\Vert = \min_{0 \leq i < n}  \lvert a_i \rvert, p = -\infty  $$
 *
 * $$ \left\Vert A \right\Vert = \( \lvert a_0 \rvert^p + \ldots + \lvert a_n \rvert^p \)^{1/p}, p>=1 $$
 *
 *
 * p = 'fro' will return Frobenius norm
 *
 * $$ \left\Vert A \right\Vert\_F = \sqrt { \sum\_{i=0}^m \sum\_{j=0}^n \lvert a\_{ij} \rvert ^2 } $$
 *
 */
declare function norm(A: NDArray, p?: number | 'fro'): number;
/**
 * @hidden
 * Perform LU decomposition
 *
 * $$ A = P L U $$
 */
declare function lu_custom(A: NDArray): NDArray;
/**
 * @hidden
 * Ref: Golub-Loan 3.1.1
 * System of equations that forms lower triangular system can be solved by
 * forward substitution.
 *   [ l00  0  ] [x0]  = [b0]
 *   [ l10 l11 ] [x1]    [b1]
 * Caller must ensure this matrix is Lower triangular before calling this
 * routine. Otherwise, undefined behavior
 */
/**
 * @hidden
 * System of equations that forms upper triangular system can be solved by
 * backward substitution.
 *   [ u00 u01 ] [x0]  = [b0]
 *   [ 0   u11 ] [x1]    [b1]
 * Caller must ensure this matrix is Upper triangular before calling this
 * routine. Otherwise, undefined behavior
 */
/**
 * @hidden
 * Apply permutation to vector
 * @param V Vector to undergo permutation (changed in place)
 * @param p Permutation vector
 */
/**
 * @hidden
 * Apply inverse permutation to vector
 * @param V Vector to undergo inverse permutation (changed in place)
 * @param p Permutation vector
 */
/**
 * Solves a system of linear scalar equations,
 * Ax = B
 * It computes the 'exact' solution for x. A is supposed to be well-
 * determined, i.e. full rank.
 * @param A Coefficient matrix (gets modified)
 * @param B RHS (populated with solution x upon return)
 */
declare function solve(A: NDArray, B: NDArray): void;
/**
 * Computes inner product of two 1D vectors (same as dot product).
 * Both inputs are supposed to be 1 dimensional arrays of same length.
 * If they are not same length, A.data.length must be <= B.data.length
 * Only first A.data.length elements of array B are used in case it's
 * longer than A
 * @param A 1D Vector
 * @param B 1D Vector
 */
declare function inner(A: NDArray, B: NDArray): number;
/**
 * Compute outer product of two vectors
 * @param A Vector of shape [m] or [m,1]
 * @param B Vector of shape [n] or [1,n]
 * @returns NDArray Matrix of dimension [m,n]
 */
declare function outer(A: NDArray, B: NDArray): NDArray;
/**
 * Perform Cholesky decomposition on given Matrix
 */
declare function cholesky(A: NDArray): NDArray;
/**
 * Singular Value Decomposition
 * Factors the given matrix A, into U,S,VT such that
 * A = U * diag(S) * VT
 * U and VT are Unitary matrices, S is 1D array of singular values of A
 * @param A Matrix to decompose Shape (m,n)
 * @param full_matrices If true, U and VT have shapes (m,m) and (n,n) resp.
 *  Otherwise the shapes are (m,k) and (k,n), resp. where k = min(m,n)
 * @param compute_uv Whether or not to compute U,VT in addition to S
 * @return [NDArray] [U,S,VT] if compute_uv = true, [S] otherwise
 */
declare function svd(A: NDArray, full_matrices?: boolean, compute_uv?: boolean): NDArray[];
/**
 * Rank of a matrix is defined by number of singular values of the matrix that
 * are non-zero (within given tolerance)
 * @param A Matrix to determine rank of
 * @param tol Tolerance for zero-check of singular values
 */
declare function rank(A: NDArray, tol?: number): number;
    /**
     */
    x: NDArray;
    /**
     * Sums of residuals; squared Euclidean 2-norm for each column in
     * Otherwise the shape is (k,).
     * TODO: WIP
     */
    residuals: NDArray;
    /**
     * Rank of coefficient matrix A
     */
    rank: number;
    /**
     * Singular values of coefficient matrix A
     */
    singulars: NDArray;
}
/**
 * Return the least-squares solution to a linear matrix equation.
 *
 * be under-, well-, or over- determined (i.e., the number of
 * the "exact" solution of the equation.
 *
 * @param A Coefficient matrix (m-by-n)
 * @param B Values on RHS of equation system. Could be array of length
 *          m or it could be 2D with dimensions m-by-k
 */
declare function lstsq(A: NDArray, B: NDArray, rcond?: number): lstsq_return;
/**
 * Compute sign and natural logarithm of the determinant of given Matrix
 * If an array has a very small or very large determinant, then a call to
 * issues, because it computes the logarithm of the determinant rather than
 * the determinant itself.
 * @param A Square matrix to compute sign and log-determinant of
 */
declare function slogdet(A: NDArray): number[];
/**
 * Compute determinant of a matrix
 * @param A Square matrix to compute determinant
 */
declare function det(A: NDArray): number;
/**
 * Compute (multiplicative) inverse of given matrix
 * @param A Square matrix whose inverse is to be found
 */
declare function inv(A: NDArray): NDArray;
/**
 * Create Lower triangular matrix from given matrix
 */
declare function tril(A: NDArray, k?: number): NDArray;
/**
 * Return Upper triangular matrix from given matrix
 */
declare function triu(A: NDArray, k?: number): NDArray;
/**
 * Compute QR decomposition of given Matrix
 */
declare function qr(A: NDArray): NDArray[];
/**
 * Compute Eigen values and left, right eigen vectors of given Matrix
 */
declare function eig(A: NDArray): NDArray[];

  }
  declare namespace geom {
declare class Axis {
    origin: NDArray;
    z: NDArray;
}
declare class CoordSystem {
    origin: NDArray;
    z: NDArray;
    x: NDArray;
    y: NDArray;
    constructor(origin: NDArray | number[], x: NDArray | number[], z: NDArray | number[]);
}

/**
 * Rational or polynomial bezier curve
 * If the weights are specified it's a rational Bezier curve
 */
declare class BezierCurve {
    degree: number;
    cpoints: NDArray;
    weights?: NDArray;
    constructor(degree: number, cpoints: NDArray, weights?: NDArray);
    /**
     * Dimension of the curve. Typically 2D or 3D
     */
    readonly dimension: number;
    /**
     * If the control points are defined in 2D plane, then add z=0 to each
     * of them to define them in 3D space
     */
    to3D(): void;
    /**
     * Is this Rational Bezier Curve
     */
    isRational(): boolean;
    /**
     * Evaluate the Bezier curve at given parameter value
     */
    evaluate(u: number, tess?: NDArray, tessidx?: number): null;
    /**
     * Tessellate the Bezier curve uniformly at given resolution
     */
    tessellate(resolution?: number): NDArray;
    /**
     * The curve is subdivided into two curves at the mipoint of parameter
     * range. This is done recursively until the curve becomes a straight line
     * within given tolerance.
     * The subdivision involves reparameterizing the curve, which is done using
     * blossoming or deCasteljau formula.
     */
    private static tessBezier(bezcrv, tolerance);
    /**
     * Tessellate bezier curve adaptively, within given tolerance of error
     */
    tessellateAdaptive(tolerance?: number): NDArray;
    /**
     * Checks if this Bezier curve is approximately a straight line within
     * given tolerance.
     */
    isLine(tolerance?: number): boolean;
    /**
     * Reparameterize the bezier curve within new parametric interval.
     * It uses the blossoming technique.
     */
    reparam(ua: number, ub: number): void;
    aabb(): AABB;
    clone(): BezierCurve;
    /**
     * Split into two Bezier curves at given parametric value
     */
    split(uk: number): BezierCurve[];
    toString(): string;
}
/**
 * Rational BSpline Curve
 */
declare class BSplineCurve {
    degree: number;
    cpoints: NDArray;
    knots: NDArray;
    weights?: NDArray;
    constructor(degree: number, cpoints: NDArray, knots: NDArray, weights?: NDArray);
    /**
     * Determines how many dimension the curve occupies based on shape of
     * Control points array
     */
    readonly dimension: number;
    /**
     * Convert 2D control points to 3D
     */
    to3D(): void;
    /**
     * Split the curve at given parameter value and return two bspline
     * curves. The two curves put together will exactly represent the
     * original curve.
     */
    split(uk: number): BSplineCurve[];
    /**
     * Replace the knots of this BSplineCurve with new knots
     */
    setKnots(knots: NDArray): void;
    /**
     * Set the knot at given index in the knot vector
     */
    setKnot(index: number, knot: number): void;
    /**
     * Set the weight at given index
     */
    setWeight(index: number, weight: number): void;
    /**
     * Is this Rational BSpline Curve. Determined based on whether weights
     * were specified while constructing this BSplineCurve
     */
    isRational(): boolean;
    /**
     * Evaluate basis function derivatives upto n'th
     */
    private evaluateBasisDerivatives(span, n, t);
    private evaluateBasis(span, t);
    private findSpan(t);
    private getTermDenominator(span, N);
    /**
     * Tesselate basis functions uniformly at given resolution
     */
    tessellateBasis(resolution?: number): NDArray;
    private static tessBSpline(bcrv, tolerance);
    /**
     * Tessellate this BSplineCurve adaptively within given tolerance of error
     */
    tessellateAdaptive(tolerance?: number): NDArray;
    /**
     * Checks if this Bezier curve is approximately a straight line within
     * given tolerance.
     */
    isLine(tolerance?: number): boolean;
    /**
     * Inserts knot un in the knot vector r-times
     * Algorithm A5.1 from "The NURBS Book"
     */
    insertKnot(un: number, r: number): void;
    /**
     * Inserts multiple knots into the knot vector at once
     * Algorithm A5.4 from "The NURBS Book"
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    refineKnots(ukList: number[]): void;
    /**
     * Algorithm A5.6 from "The NURBS Book"
     * The total number of bezier segments required to decompose a
     * given bspline curve
     *  = Number of internal knots + 1
     *  = Length of knot vector - 2*(p+1) + 1
     *  = (m+1) - 2*(p+1) + 1
     *  = m - 2*p
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    decompose(): BezierCurve[];
    /**
     * Evaluate the BSplineCurve at given parameter value
     * euclidean point is returned.
     */
    evaluate(t: number, tess?: NDArray, tessidx?: number): NDArray | null;
    /**
     * Evaluate the derivative of BSplineCurve at given parameter value
     * euclidean point is returned.
     */
    evaluateDerivative(t: number, d: number, tess?: NDArray, tessidx?: number): NDArray | null;
    /**
     * Tessellate the BSplineCurve uniformly at given resolution
     */
    tessellate(resolution?: number): NDArray;
    /**
     * Tessellate derivatives of BSplineCurve uniformly at given resolution
     */
    tessellateDerivatives(resolution: number | undefined, d: number): NDArray;
    clone(): BSplineCurve;
    aabb(): AABB;
    toString(): string;
}

declare class BezierSurface {
    u_degree: number;
    v_degree: number;
    cpoints: NDArray;
    weights?: NDArray;
    constructor(u_degree: number, v_degree: number, cpoints: NDArray | number[][][], weights?: NDArray | number[][]);
    readonly dimension: number;
    isRational(): boolean;
    evaluate(u: number, v: number, tess: NDArray, uidx: number, vidx: number): void;
    tessellatePoints(resolution?: number): NDArray;
    tessellate(resolution?: number): {
        vertices: TypedArray;
        faces: number[];
    };
    aabb(): AABB;
    clone(): BezierSurface;
}
declare class BSplineSurface {
    u_degree: number;
    v_degree: number;
    /**
     * cpoints is a two dimensional grid of coordinates.
     * The outermost index is along U, the next inner index is along V
     *
     *          V-->
     *      [
     *  U     [ [xa,ya,za], [xb,yb,zb], ...]
     *  |     [ [xl,yl,zl], [xm,ym,zm], ...]
     *  |     .
     *  v     .
     *      ]
     */
    cpoints: NDArray;
    u_knots: NDArray;
    v_knots: NDArray;
    weights?: NDArray;
    constructor(u_degree: number, v_degree: number, u_knots: NDArray | number[], v_knots: NDArray | number[], cpoints: NDArray | number[][][], weights?: NDArray | number[][]);
    readonly dimension: number;
    clone(): BSplineSurface;
    aabb(): AABB;
    isRational(): boolean;
    isFlat(tolerance?: number): boolean;
    setUKnots(u_knots: NDArray): void;
    setVKnots(v_knots: NDArray): void;
    evaluate(u: number, v: number, tess: NDArray, uidx: number, vidx: number): void;
    tessellatePoints(resolution?: number): NDArray;
    tessellate(resolution?: number): {
        vertices: TypedArray;
        faces: number[];
    };
    static tessellateRecursive(bsrf: BSplineSurface, tolerance?: number): {
        vertices: TypedArray;
        faces: number[];
    }[];
    tessellateAdaptive(tolerance?: number): {
        vertices: TypedArray;
        faces: number[];
    }[];
    /**
     * Split this BSplineSurface into two at uk, by refining u-knots
     */
    splitU(uk: number): BSplineSurface[];
    /**
     * Split this BSplineSurface into two at vk, by refining v-knots
     */
    splitV(vk: number): BSplineSurface[];
    /**
     * Split this BSplineSurface into four
     */
    splitUV(uk: number, vk: number): BSplineSurface[];
    /**
     * Inserts knot un in the U knot vector r-times
     * Ref: Algorithm A5.3 "The NURBS book"
     * @param un Knot to be inserted
     * @param r Number of times to insert the knot
     */
    insertKnotU(un: number, r: number): void;
    /**
     * Inserts knot vn in the V knot vector r-times
     * Ref: Algorithm A5.3 "The NURBS book"
     * @param vn Knot to be inserted
     * @param r Number of times to insert the knot
     */
    insertKnotV(vn: number, r: number): void;
    /**
     * Insert knots in U and V knot vectors
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    insertKnotUV(un: number, vn: number, ur: number, vr: number): void;
    /**
     * Inserts multiple knots into the U knot vector at once
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    refineKnotsU(uklist: number[]): void;
    /**
     * Inserts multiple knots into the V knot vector at once
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    refineKnotsV(vklist: number[]): void;
    /**
     * Inserts multiple knots into the U and V knot vectors at once
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    refineKnotsUV(uklist: number[], vklist: number[]): void;
    decomposeU(): NDArray;
    decomposeV(): NDArray;
    /**
     * Creates grid of Bezier surfaces that represent this BSpline surface.
     * The routine first computes bezier strips along u (i.e. BSpline surfaces that
     * are Bezier in one direction and BSpline in other). Subsequently
     * decompose it called on each of these strips in the v direction
     * Algorithm A5.7 from "The NURBS Book"
     * See http://www.bluemathsoftware.com/pages/nurbs/funalgo
     */
    decompose(): BezierSurface[];
    toString(): string;
}

declare class LineSegment extends BSplineCurve {
    constructor(from: number[], to: number[]);
}
declare class CircleArc extends BSplineCurve {
    constructor(coordsys: CoordSystem, radius: number, start: number, end: number);
}
declare class Circle extends CircleArc {
    constructor(coord: CoordSystem, radius: number);
}

/**
 * @hidden
 * Compute all n'th degree bernstein polynomials at given parameter value
 */
declare function bernstein(n: number, u: number): Array<number>;
/**
 * @hidden
 * @param {number} p Degree
 * @param {Array.<number>} U Knot vector
 * @param {number} u Parameter
 * @returns {number}
 */
declare function findSpan(p: number, U: Array<number> | TypedArray, u: number): number;
/**
 * @hidden
 * Evaluate basis function values
 * @param {number} p Degree
 * @param {Array.<number>} U Knot vector
 * @param {number} i Knot span index
 * @param {number} u Parameter
 * @returns {Array} Basis function values at i,u
 */
declare function getBasisFunction(p: number, U: Array<number> | TypedArray, i: number, u: number): Array<number>;
/**
 * @hidden
 * The NURBS book Algo A2.3
 * Compute non-zero basis functions and their derivatives, upto and including
 * @param {number} p Degree
 * @param {number} u Parameter
 * @param {number} i Knot span
 * @param {NDArray} knots Knot vector
 * @param {number} n nth derivative
 * @returns {NDArray} ders ders[k][j] is k'th derivative of
 *            basic function N(i-p+j,p), where 0<=k<=n and 0<=j<=p
 */
declare function getBasisFunctionDerivatives(p: number, u: number, ki: number, knots: NDArray, n: number): NDArray;
declare function blossom(cpoints: NDArray, n: number, ts: number[]): NDArray;
/**
 * Computes equation of plane passing through given 3 points
 * Eqn of plane is
 *  ax + by + cz + d = 0
 * This routine returns [a,b,c,d]
 * The direction of the normal is defined by assuming that A,B,C are in
 * counter-clockwise direction. i.e. if you curl fingers of right hand
 * in counter-clockwise direction, then the raised thumb will give the
 * direction of the plane normal
 */
declare function planeFrom3Points(A: NDArray, B: NDArray, C: NDArray): number[];
/**
 * Finds intersection between two line segments in 3D
 * First line segment extends from p1 to p2, and second extends from p3 to p4
 * Input points are assumed to be coordinates in 3D coord system
 *
 * Algorithm based on C implemention by Paul Bourke
 * http://paulbourke.net/geometry/pointlineplane/lineline.c
 *
 * The method will return a tuple with results (ua, ub),
 * where u is the parameter value on each line
 * If there is no intersection null will be returned
 */
declare function intersectLineSegLineSeg3D(p1: number[], p2: number[], p3: number[], p4: number[]): null | number[];
/**
 * The check works by constructing a line between first and last control
 * point and then finding the distance of other control points from this
 * line. Instead of actually calculating the distance from the line, we
 * do the check if the point lies on the line or not. This is done by
 * substituting the [x,y] coordinates of control point, into the equation
 * of the line. If the result is zero within the tolerance value, then
 * the control point lies on the line. If all control points lie on the line
 * then the curve can be considered a straight line.
 * @param points Array of points in 2D coord
 * @param tolerance Tolerance within which a group of points is co-linear
 */
declare function arePointsColinear(points: NDArray, tolerance: number): boolean;


declare class BilinearSurface extends BSplineSurface {
    constructor(p00: number[], p01: number[], p10: number[], p11: number[]);
}
declare class GeneralCylinder extends BSplineSurface {
    constructor(curve: BSplineCurve, direction: NDArray | number[], height: number);
}
declare class Cylinder extends GeneralCylinder {
    constructor(coordsys: CoordSystem, radius: number, height: number);
}
declare class RuledSurface extends BSplineSurface {
}
declare class RevolutionSurface extends BSplineSurface {
}
declare class Cone extends RevolutionSurface {
}
declare class Sphere extends RevolutionSurface {
}
declare class Torus extends RevolutionSurface {
}

  }
  declare namespace topo {
declare class Body {
    vertices: Vertex[];
    halfedges: HalfEdge[];
    edges: Edge[];
    faces: Face[];
    id: string;
    constructor();
    newFace(): Face;
    newVertex(coord: NDArray): Vertex;
    newHalfEdge(): HalfEdge;
    newEdge(): Edge;
    removeEdge(edge: Edge): void;
    removeVertex(vertex: Vertex): void;
    removeHalfEdge(halfEdge: HalfEdge): void;
    removeFace(face: Face): void;
    unlink(): void;
    toDOT(): string;
    toSVG(width?: number, height?: number): string;
}

declare class Edge {
    hePlus?: HalfEdge;
    heMinus?: HalfEdge;
    curve: any;
    id: string;
    constructor();
    startVertex(): Vertex;
    endVertex(): Vertex;
    unlink(): void;
}

declare type MVFS_result = {
    vertex: Vertex;
    body: Body;
    face: Face;
};
declare type MEV_result = {
    edge: Edge;
    vertex: Vertex;
};
declare type MEF_result = {
    edge: Edge;
    face: Face;
};
declare class EulerOps {
    /**
     * Make Vertex Face Solid
     * (Solid = Body in this library)
     */
    static MVFS(coord: NDArray): MVFS_result;
    /**
     * Kill Vertex Face Solid
     * (Solid = Body in this library)
     */
    static KVFS(body: Body): void;
    /**
     * Low level MEV (Make Edge Vertex)
     */
    private static LMEV(he0, he1, coord);
    /**
     * Make Edge Vertex
     */
    static MEV(face: Face, vertex: Vertex, coord: NDArray): {
        vertex: Vertex;
        edge: Edge;
    };
    /**
     * Kill Edge Vertex
     */
    static KEV(edge: Edge, vertex: Vertex): void;
    /**
     * Make Edge Face
     */
    static MEF(face: Face, fromHEV0: Vertex, fromHEV1: Vertex | undefined, toHEV0: Vertex, toHEV1: Vertex | undefined): MEF_result;
    /**
     * Kill Edge Face
     */
    static KEF(edge: Edge, face: Face): void;
    /**
     * Kill Edge Make Ring
     * (Ring = Loop in this library)
     */
    static KEMR(face: Face, v1: Vertex, v2: Vertex): Loop;
    /**
     * Make Edge Kill Ring
     * (Ring = Loop in this library)
     */
    static MEKR(faceFrom: Face, fromHEV0: Vertex, fromHEV1: Vertex, faceTo: Face, toHEV0: Vertex, toHEV1: Vertex): Edge;
}

declare class Face {
    oloop?: Loop;
    iloops: Loop[];
    body: Body;
    surface: any;
    id: string;
    constructor(body: Body);
    addLoop(loop: Loop): void;
    removeLoop(loop: Loop): void;
    setOuterloop(loop: Loop): void;
    unlink(): void;
    findHalfEdge(vtxFrom: Vertex, vtxTo?: Vertex): HalfEdge | undefined;
}

declare type heWalkHandler = (he: HalfEdge, count: number) => void;
declare class HalfEdge {
    vertex?: Vertex;
    prev?: HalfEdge;
    next?: HalfEdge;
    edge?: Edge;
    loop?: Loop;
    id: string;
    constructor(origin?: Vertex, pair?: HalfEdge, next?: HalfEdge, loop?: Loop);
    mate(): HalfEdge;
    unlink(): void;
    isSolitary(): boolean;
    prevInLoop(): HalfEdge | undefined;
    static walk(heStart: HalfEdge, callback: heWalkHandler): void;
}

declare class IDManager {
    static idmap: any;
    static init(): void;
    static genId(label: string): number;
}


declare class Loop {
    face: Face;
    halfedge?: HalfEdge;
    id: string;
    constructor(face: Face);
    insertHalfEdgeAfter(heNew: HalfEdge, heExisting: HalfEdge): void;
    removeHalfEdge(he: HalfEdge): void;
    readonly length: number;
    unlink(): void;
    toString(): string;
}


declare type walkHandler = (he: HalfEdge, count: number) => void;
declare class Vertex {
    coord?: NDArray;
    halfedge?: HalfEdge;
    id: string;
    constructor(coord?: NDArray, halfedge?: HalfEdge);
    walk(heStart: HalfEdge, callback: walkHandler): void;
    degree(): number;
    unlink(): void;
}

  }
}
`;
