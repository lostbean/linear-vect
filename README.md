# linear-vect

A low-dimensional linear algebra library for Haskell, operating on the `Num` typeclass instead of being restricted to `Float`, `Double`, or `Floating`.

Forked from the [vect](http://hackage.haskell.org/package/vect) library by Balazs Komuves, with contributions from Charles Durham, and generalized by capsjac to work over `Num`.

## Overview

`linear-vect` provides fixed-size 2D, 3D, and 4D vector and matrix types with a rich algebraic type class hierarchy. Unlike the original `vect` library, all core operations are parameterized over the `Num` typeclass, so the same vector/matrix types work with `Double`, `Float`, `Int`, or any custom numeric type. Operations that inherently require `Floating` (e.g., normalization, norms) are constrained accordingly.

The library also includes:

- Unit vector (normal) wrapper types with compile-time guarantees
- Orthogonal and projective matrix newtypes
- Non-square (rectangular) matrix types with transpose
- QR decomposition (Householder and Gram-Schmidt)
- Eigenvalue decomposition for symmetric matrices
- Hessenberg reduction
- `Unbox` instances for use with the `vector` package
- `Storable` instances for FFI interop
- `NFData` instances for deep evaluation
- `Random` instances for all types (including uniformly random orthogonal matrices via Haar measure)

## Modules

### `Linear.Class`

Defines the algebraic type class hierarchy and all operators. This is the foundation module -- it contains no data types, only class declarations and standalone utility functions.

**Algebraic classes:**

| Class | Description |
|---|---|
| `AbelianGroup` | Additive group: `(&+)`, `(&-)`, `neg`, `zero` |
| `MultSemiGroup` | Multiplicative semigroup with identity: `(.*.)`, `one` |
| `Ring` | Combines `AbelianGroup` and `MultSemiGroup` |
| `LeftModule` | Left action of a matrix on a vector: `(*.)` or `lmul` |
| `RightModule` | Right action of a vector on a matrix: `(.*)` or `rmul` |
| `LinearMap` | Scalar-vector multiplication and mapping: `(*&)`, `(&*)`, `scalarMul`, `mapVec` |

**Geometric classes:**

| Class | Description |
|---|---|
| `DotProd` | Dot product: `(&.)`, `dotprod`, `normsqr`, `lensqr` |
| `Norm` | Euclidean norm: `norm`, `vlen` (requires `Floating`) |
| `CrossProd` | Cross product: `(&^)`, `crossprod` |
| `Pointwise` | Component-wise (Hadamard) product: `(&!)`, `pointwise` |
| `Tensor` | Outer product: `outer` |

**Matrix classes:**

| Class | Description |
|---|---|
| `Transpose` | Matrix transpose: `transpose` |
| `SquareMatrix` | Identity matrix and inverse: `idmtx`, `inverse` |
| `Diagonal` | Construct diagonal matrices / extract diagonal: `diag`, `diagVec` |
| `Determinant` | Determinant: `det` |
| `Orthogonal` | Wrap/unwrap orthogonal matrix newtypes: `fromOrtho`, `toOrthoUnsafe` |
| `Projective` | Affine/projective transformations: `fromProjective`, `toProjectiveUnsafe`, `orthogonal`, `linear`, `translation`, `scaling` |
| `MatrixNorms` | Frobenius norm and matrix distance: `frobeniusNorm`, `matrixDistance` |

**Dimension conversion:**

| Class | Description |
|---|---|
| `Extend` | Convert between dimensions: `extendZero`, `extendWith`, `trim`, `extendHeadZero`, `extendHeadWith`, `trimHead` |
| `Dimension` | Query dimensionality: `dim` |

**Element access classes:**

| Class | Method | Description |
|---|---|---|
| `HasOne` | `_1` | First component |
| `HasTwo` | `_2` | Second component |
| `HasThree` | `_3` | Third component |
| `HasFour` | `_4` | Fourth component |

**Other classes:**

| Class | Description |
|---|---|
| `UnitVector` | Create/unwrap unit vectors: `mkNormal`, `fromNormal`, `toNormalUnsafe`, `fromNormalRadius` |
| `NearZero` | Epsilon-based near-zero testing: `epsilon`, `isMainlyZero` |
| `PrettyShow` | Pretty-print vectors/matrices: `showPretty` |

**Standalone utility functions:**

| Function | Description |
|---|---|
| `normalize` | Normalize a vector to unit length |
| `distance` | Euclidean distance between two vectors |
| `angle` / `angle'` | Angle between two vectors (or unit vectors) |
| `acosSafe` | Safe `acos` that tolerates small rounding past +/-1 |
| `project` / `project'` | Project vector onto hyperplane orthogonal to another |
| `flipNormal` | Negate a unit vector |
| `householder` | Householder reflection matrix from a unit vector |

### `Linear.Vect`

Concrete vector data types and their instances. Re-exports `Linear.Class`.

**Vector types:**

| Type | Constructor | Description |
|---|---|---|
| `Vec2 a` | `Vec2 !a !a` | 2D vector, strict in both components |
| `Vec3 a` | `Vec3 !a !a !a` | 3D vector, strict in all components |
| `Vec4 a` | `Vec4 !a !a !a !a` | 4D vector, strict in all components |

**Unit vector (normal) types:**

| Type | Wraps | Description |
|---|---|---|
| `Normal2 a` | `Vec2 a` | Unit-length 2D vector |
| `Normal3 a` | `Vec3 a` | Unit-length 3D vector |
| `Normal4 a` | `Vec4 a` | Unit-length 4D vector |

**Type aliases for `Double`:**

```haskell
type Vec2D    = Vec2 Double
type Vec3D    = Vec3 Double
type Vec4D    = Vec4 Double
type Normal2D = Normal2 Double
type Normal3D = Normal3 Double
type Normal4D = Normal4 Double
```

**Tuple conversion:** `mkVec2`, `mkVec3`, `mkVec4`, `unVec2`, `unVec3`, `unVec4`

**Component accessors:** `_x`, `_y`, `_z`, `_w`

### `Linear.Mat`

Matrix data types and their instances. Re-exports `Linear.Class`.

**Square matrix types** (row-major; components are row vectors):

| Type | Constructor | Description |
|---|---|---|
| `Mat2 a` | `Mat2 !(Vec2 a) !(Vec2 a)` | 2x2 matrix |
| `Mat3 a` | `Mat3 !(Vec3 a) !(Vec3 a) !(Vec3 a)` | 3x3 matrix |
| `Mat4 a` | `Mat4 !(Vec4 a) !(Vec4 a) !(Vec4 a) !(Vec4 a)` | 4x4 matrix |

**Rectangular matrix types:** `Mat2x3`, `Mat2x4`, `Mat3x2`, `Mat3x4`, `Mat4x2`, `Mat4x3`

**Orthogonal matrix newtypes:**

| Type | Wraps | Description |
|---|---|---|
| `Ortho2 a` | `Mat2 a` | Orthogonal 2x2 matrix (inverse = transpose) |
| `Ortho3 a` | `Mat3 a` | Orthogonal 3x3 matrix |
| `Ortho4 a` | `Mat4 a` | Orthogonal 4x4 matrix |

**Projective matrix newtypes:**

| Type | Wraps | Description |
|---|---|---|
| `Proj3 a` | `Mat3 a` | Projective 2D transformations |
| `Proj4 a` | `Mat4 a` | Projective 3D transformations |

### `Linear.Decomp`

Matrix decomposition algorithms. Re-exports `Linear.Class`.

| Function | Description |
|---|---|
| `qrHouse` | QR decomposition via Householder reflections. Returns `(Q, R)`. |
| `qrGram` | QR decomposition via modified Gram-Schmidt. |
| `symmEigen` | Eigenvalue decomposition for symmetric matrices via iterated QR. Returns `(eigenvectors, eigenvalues)`. |
| `hessen` | Reduce a matrix to Hessenberg form (tridiagonal for symmetric matrices). Instances for `Mat3` and `Mat4`. |

## Operator Quick Reference

| Operator | Fixity | Meaning | Example |
|---|---|---|---|
| `&+` | `infixl 6` | Vector/matrix addition | `v1 &+ v2` |
| `&-` | `infixl 6` | Vector/matrix subtraction | `v1 &- v2` |
| `&.` | `infix 7` | Dot product | `v1 &. v2` |
| `&^` | `infix 7` | Cross product | `v1 &^ v2` |
| `&!` | `infix 7` | Pointwise (Hadamard) product | `v1 &! v2` |
| `*&` | `infixr 7` | Scalar times vector (left) | `2.0 *& v` |
| `&*` | `infixl 7` | Vector times scalar (right) | `v &* 0.5` |
| `.*.` | `infixl 7` | Matrix-matrix multiplication | `m1 .*. m2` |
| `*.` | `infixr 7` | Matrix-vector multiplication | `m *. v` |
| `.*` | `infixl 7` | Vector-matrix multiplication | `v .* m` |

## Usage Examples

### Basic Vector Arithmetic

```haskell
import Linear.Vect

let a = Vec3 1.0 2.0 3.0
    b = Vec3 4.0 5.0 6.0

a &+ b           -- Vec3 5.0 7.0 9.0
a &- b           -- Vec3 (-3.0) (-3.0) (-3.0)
2.0 *& a         -- Vec3 2.0 4.0 6.0
a &. b           -- 32.0
a &^ b           -- Vec3 (-3.0) 6.0 (-3.0)
norm a           -- 3.7416...
normalize a      -- Vec3 0.2672... 0.5345... 0.8017...
distance a b     -- 5.196...
```

### Unit Vectors

```haskell
import Linear.Vect

let v = Vec3 3.0 4.0 0.0
    n = mkNormal v       -- Normal3, automatically normalized
fromNormal n             -- Vec3 0.6 0.8 0.0
```

### Matrix Operations

```haskell
import Linear.Mat

let m = Mat3 (Vec3 1 0 0) (Vec3 0 2 0) (Vec3 0 0 3) :: Mat3 Double
m .*. m          -- squares the matrix
m *. Vec3 1 1 1  -- Vec3 1.0 2.0 3.0
transpose m
inverse m
det m            -- 6.0
```

### Dimension Conversion

```haskell
import Linear.Vect

let v2 = Vec2 5.0 6.0
extendZero v2        -- Vec3 5.0 6.0 0.0
extendWith 1.0 v2    -- Vec3 5.0 6.0 1.0

let v4 = Vec4 1.0 2.0 3.0 4.0
trim v4              -- Vec3 1.0 2.0 3.0
```

### QR Decomposition

```haskell
import Linear.Decomp

let m = Mat3 (Vec3 12 (-51) 4) (Vec3 6 167 (-68)) (Vec3 (-4) 24 (-41)) :: Mat3 Double
let (q, r) = qrHouse m
-- q is orthogonal, r is upper triangular
-- m == q .*. r (up to floating point precision)
```

### Eigenvalue Decomposition (Symmetric Matrices)

```haskell
import Linear.Decomp

let m = Mat3 (Vec3 2 1 0) (Vec3 1 3 (-1)) (Vec3 0 (-1) 6) :: Mat3 Double
let (eigvecs, eigvals) = symmEigen m
```

## Notes and Limitations

- **Mat4 inverse:** `inverse` for `Mat4` is not implemented and will throw a runtime error.
- **Mat4 determinant:** `det` for `Mat4` is not implemented and will throw a runtime error.
- **Row-major storage:** Matrix constructors take row vectors. `Mat3 r1 r2 r3` means `r1` is the first row.
- **Random orthogonal matrices:** `Random` instances for `Ortho2`, `Ortho3`, `Ortho4` generate orientation-preserving (det = +1) orthogonal matrices using the Haar measure on SO(n).

## Dependencies

- `base` >= 4.7 && < 4.9
- `deepseq` >= 1.2
- `random` >= 1.0 && < 1.2
- `vector` >= 0.10
- `vector-th-unbox` >= 0.2.1

## Building

```bash
# With Stack
stack build

# With Cabal
cabal build
```

## Credits

- **Balazs Komuves** -- original [vect](http://hackage.haskell.org/package/vect) library (2008-2011)
- **Charles Durham** -- contributions (2014)
- **capsjac** -- `Num` generalization [fork](https://github.com/capsjac/linear-vect) (2014)

## Author

Edgar Gomes de Araujo (<talktoedgar@gmail.com>)

## License

BSD-3-Clause -- see [LICENSE](./LICENSE).
