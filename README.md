# linear-vect

Low-dimensional linear algebra for Haskell -- vectors, matrices, and decompositions in 2D, 3D, and 4D.

## What is this?

A Haskell library providing fixed-size vector and matrix types for 2D, 3D, and 4D, along with matrix decomposition algorithms. Forked from the [vect](http://hackage.haskell.org/package/vect) library by Balazs Komuves, with contributions from Charles Durham and capsjac.

The key difference from the original: all core operations work over the `Num` typeclass instead of being restricted to `Float` or `Double`, so the same vector/matrix types work with `Int`, custom numeric types, or any `Num` instance. Operations that inherently require `Floating` (normalization, norms) are constrained accordingly.

## What does it provide?

- **Vectors and matrices** in 2, 3, and 4 dimensions with standard algebraic operations: addition, subtraction, scalar multiplication, dot product, cross product (3D), pointwise (Hadamard) product, and outer product. Both square and rectangular matrix types are included.

- **Unit vector wrappers** (Normal types) that guarantee a vector has length 1, enforced at the type level. Useful when an algorithm requires normalized directions.

- **Orthogonal and projective matrix newtypes** that encode invariants in the type system. For example, an orthogonal matrix's inverse is its transpose -- the newtype ensures this is used correctly. Projective matrices represent affine transformations (translation, scaling, rotation) in homogeneous coordinates.

- **QR decomposition** -- two implementations that factor a matrix A into Q (orthogonal) and R (upper triangular):
  - *Householder reflections*: reflects columns onto coordinate axes using orthogonal reflection matrices. Numerically stable and the preferred method for most applications.
  - *Gram-Schmidt orthogonalization*: iteratively projects out already-found orthogonal directions to build Q column by column. Conceptually simpler but less numerically stable in floating-point arithmetic.

- **Eigenvalue decomposition** for symmetric matrices via the iterated QR algorithm with Hessenberg reduction. Returns eigenvalues (the scaling factors along principal axes) and eigenvectors (the principal axis directions). This is used downstream for fitting Bingham distributions and analyzing scatter matrices.

- **Hessenberg reduction**: transforms a matrix into upper Hessenberg form (all zeros below the first subdiagonal) using Householder transformations. This is a preprocessing step that speeds up QR iteration by reducing the work per step.

- **Random instances**: all vector, matrix, and normal types can be randomly generated. Orthogonal matrices are sampled uniformly from SO(n) using the Haar measure (Mezzadri's algorithm), which produces random rotations with no directional bias -- important for Monte Carlo simulations involving orientations.

## Example

```haskell
import Linear.Vect
import Linear.Mat
import Linear.Decomp

let a = Vec3 1.0 2.0 3.0
    b = Vec3 4.0 5.0 6.0

a &+ b           -- Vec3 5.0 7.0 9.0
a &. b           -- 32.0  (dot product)
a &^ b           -- Vec3 (-3.0) 6.0 (-3.0)  (cross product)
2.0 *& a         -- Vec3 2.0 4.0 6.0
normalize a      -- unit vector in the direction of a

let m = Mat3 (Vec3 2 1 0) (Vec3 1 3 (-1)) (Vec3 0 (-1) 6) :: Mat3 Double
let (q, r) = qrHouse m                -- QR decomposition
let (eigvecs, eigvals) = symmEigen m   -- eigenvalue decomposition
```

## Where is it used?

linear-vect is the foundational linear algebra library for this ecosystem. It is a dependency of:

- **DeUni** -- Delaunay triangulation (uses Vec2/Vec3 for point geometry, QR decomposition for circumsphere computation)
- **hammer** -- microstructure analysis (matrices, vector operations, N-dimensional arrays)
- **sledge** -- crystallographic texture analysis (quaternions are Vec4, rotation matrices are Mat3, eigendecomposition for Bingham fitting)
- **SubZero** -- subdivision surfaces (vertex positions as Vec2D/Vec3D)
- **VirMat** -- virtual microstructure generation (transitively through all of the above)

## How to build

```bash
# With Nix (recommended)
nix develop
cabal build --allow-newer

# With Cabal
cabal build

# Run tests
cabal test
```

## Credits

- **Balazs Komuves** -- original [vect](http://hackage.haskell.org/package/vect) library
- **Charles Durham** -- contributions
- **capsjac** -- `Num` generalization [fork](https://github.com/capsjac/linear-vect)

## License

BSD-3-Clause -- see [LICENSE](./LICENSE).
