{-# LANGUAGE FunctionalDependencies #-}

-- | class declarations
module Linear.Class
  ( AbelianGroup(..) , vecSum
  , MultSemiGroup(..) , Ring , semigroupProduct
  , LeftModule(..) , RightModule(..)
  , LinearMap(..) , DotProd(..) , Norm(..) , CrossProd(..)
  , normalize , distance , angle , angle'
  , UnitVector(..)
  , Pointwise(..)
  , Extend(..) , Dimension(..) , Transpose(..)
  , SquareMatrix(..) , Tensor(..) , Diagonal (..) , Determinant(..)
  , Orthogonal(..) , Projective(..) , MatrixNorms(..)
  , project , project' , projectUnsafe , flipNormal
  , householder, householderOrtho
  , HasOne   (..)
  , HasTwo   (..)
  , HasThree (..)
  , HasFour  (..)
  )
where

class AbelianGroup g where
  (&+) :: g -> g -> g
  (&-) :: g -> g -> g
  neg  :: g -> g
  zero :: g

infixl 6 &+
infixl 6 &-

vecSum :: AbelianGroup g => [g] -> g
vecSum l = foldl (&+) zero l

class MultSemiGroup r where
  (.*.) :: r -> r -> r
  one   :: r

class (AbelianGroup r, MultSemiGroup r) => Ring r

infixl 7 .*.

-- was: ringProduct :: Ring r => [r] -> r
semigroupProduct :: MultSemiGroup r => [r] -> r
semigroupProduct l = foldl (.*.) one l

class LeftModule r m | r -> m, m -> r where
  lmul :: r -> m -> m
  (*.) :: r -> m -> m
  (*.) = lmul

class RightModule m r | m -> r, r -> m where
  rmul :: m -> r -> m
  (.*) :: m -> r -> m
  (.*) = rmul

-- I'm not really sure about this.. may actually degrade the performance in some cases?
{- RULES
"matrix multiplication left"   forall m n x.  (n .*. m) *. x = n *. (m *. x)
"matrix multiplication right"  forall m n x.  x .* (m .*. n) = (x .* m) .* n
  -}

infixr 7 *.
infixl 7 .*

class AbelianGroup (v a) => LinearMap a v where
  mapVec    :: (a -> a) -> v a -> v a
  scalarMul :: a -> v a -> v a
  (*&) ::      a -> v a -> v a
  (&*) ::      v a -> a -> v a
  (*&) s v = scalarMul s v
  (&*) v s = scalarMul s v

infixr 7 *&
infixl 7 &*

{-# RULES
"scalar multiplication left"   forall (s :: Num s => s) (t :: Num t => t) x. t *& (s *& x) = (t*s) *& x
"scalar multiplication right"  forall (s :: Num s => s) (t :: Num t => t) x.  (x &* s) &* t = x &* (s*t)
  #-}

class Num a => DotProd a v where
  (&.) :: v a -> v a -> a
  dotprod :: v a -> v a -> a
  dotprod = (&.)
  normsqr :: v a -> a
  normsqr v = (v &. v)
  lensqr  :: v a -> a
  lensqr = normsqr

class (Floating a, DotProd a v) => Norm a v where
  norm    :: v a -> a
  norm = sqrt . lensqr
  vlen     :: v a -> a
  vlen = norm

infix 7 &.

{-# RULES
"vlen/square 1"   forall x.  (vlen x)*(vlen x) = lensqr x
"vlen/square 2"   forall x.  (vlen x)^2 = lensqr x
"norm/square 1"  forall x.  (norm x)*(norm x) = normsqr x
"norm/square 2"  forall x.  (norm x)^2 = normsqr x
  #-}


normalize :: (LinearMap a v, Norm a v) => v a -> v a
normalize v = scalarMul (recip (vlen v)) v

distance :: (LinearMap a v, Norm a v) => v a -> v a -> a
distance x y = norm (x &- y)

-- | the angle between two vectors
angle :: (LinearMap a v, Norm a v) => v a -> v a -> a
angle x y = acos $ (x &. y) / (norm x * norm y)

-- | the angle between two unit vectors
angle' {- ' CPP is sensitive to primes -} :: (Floating a, LinearMap a v, UnitVector a v u, DotProd a v) => u a -> u a -> a
angle' x y = acos (fromNormal x &. fromNormal y)

{-# RULES
"normalize is idempotent"  forall x. normalize (normalize x) = normalize x
  #-}

class (LinearMap a v, Norm a v) => UnitVector a v u | u -> v, v -> u where
  mkNormal         :: v a -> u a       -- ^ normalizes the input
  toNormalUnsafe   :: v a -> u a       -- ^ does not normalize the input!
  fromNormal       :: u a -> v a
  fromNormalRadius :: a -> u a -> v a
  fromNormalRadius t n = t *& fromNormal n

-- | Projects the first vector down to the hyperplane orthogonal to the second (unit) vector
project' :: (LinearMap a v, UnitVector a v u, Norm a v) => v a -> u a -> v a
project' what dir = projectUnsafe what (fromNormal dir)

-- | Direction (second argument) is assumed to be a /unit/ vector!
projectUnsafe :: (LinearMap a v, DotProd a v) => v a -> v a -> v a
projectUnsafe what dir = what &- dir &* (what &. dir)

project :: (Fractional a, LinearMap a v, DotProd a v) => v a -> v a -> v a
project what dir = what &- dir &* ((what &. dir) / (dir &. dir))

-- | Since unit vectors are not a group, we need a separate function.
flipNormal :: UnitVector a v n => n a -> n a
flipNormal = toNormalUnsafe . neg . fromNormal

-- | Cross product
class CrossProd v where
  crossprod :: v -> v -> v
  (&^)      :: v -> v -> v
  (&^) = crossprod

-- | Pointwise multiplication
class Pointwise v where
  pointwise :: v -> v -> v
  (&!)      :: v -> v -> v
  (&!) = pointwise

infix 7 &^
infix 7 &!

-- | conversion between vectors (and matrices) of different dimensions
class Extend a u v where
  extendZero     :: u a -> v a      -- ^ example: @extendZero (V2 5 6) = V4 5 6 0 0@
  extendWith     :: a -> u a -> v a -- ^ example: @extendWith 1 (V2 5 6) = V4 5 6 1 1@
  trim           :: v a -> u a      -- ^ example: @trim (V4 5 6 7 8) = V2 5 6@
  extendHeadZero :: u a -> v a      -- ^ example: @extendHeadZero (Vec2 5 6) = Vec4 0 0 5 6@
  extendHeadWith :: a -> u a -> v a -- ^ example: @extendHeadWith 1 (Vec2 5 6) = Vec4 1 1 5 6@
  trimHead       :: v a -> u a      -- ^ example: @trimHead (Vec4 5 6 7 8) = Vec2 7 8@

-- | makes a diagonal matrix from a vector
class Diagonal s t | t -> s where
  diag    :: s -> t
  diagVec :: t -> s

class Transpose m n | m -> n, n -> m where
  transpose :: m -> n

class SquareMatrix m where
  inverse :: m -> m
  idmtx :: m

{-# RULES
"transpose is an involution"  forall m. transpose (transpose m) = m
"inverse is an involution"    forall m. inverse (inverse m) = m
  #-}

class SquareMatrix (m a) => Orthogonal a m o | m -> o, o -> m where
  fromOrtho     :: o a -> m a
  toOrthoUnsafe :: m a -> o a

class (AbelianGroup m, SquareMatrix m) => MatrixNorms a m where
  frobeniusNorm  :: m -> a       -- ^ the frobenius norm (= euclidean norm in the space of matrices)
  matrixDistance :: m -> m -> a  -- ^ euclidean distance in the space of matrices
  operatorNorm   :: m -> a      -- ^ (euclidean) operator norm (not implemented yet)
  matrixDistance m n = frobeniusNorm (n &- m)
  operatorNorm = error "operatorNorm: not implemented yet"

-- | Outer product (could be unified with Diagonal?)
class Tensor t v | t -> v where
  outer :: v -> v -> t

class Determinant a m where
  det :: m -> a

class Dimension a where
  dim :: a -> Int

-- | Householder matrix, see <http://en.wikipedia.org/wiki/Householder_transformation>.
-- In plain words, it is the reflection to the hyperplane orthogonal to the input vector.
householder :: (LinearMap a v, UnitVector a v u, SquareMatrix (m a), LinearMap a m, Tensor (m a) (v a)) => u a -> m a
householder u = idmtx &- (2 *& outer v v)
  where v = fromNormal u

householderOrtho :: (LinearMap a v, UnitVector a v u, SquareMatrix (m a), LinearMap a m, Tensor (m a) (v a), Orthogonal a m o) => u a -> o a
householderOrtho = toOrthoUnsafe . householder

-- | \"Projective\" matrices have the following form: the top left corner
-- is an any matrix, the bottom right corner is 1, and the top-right
-- column is zero. These describe the affine orthogonal transformation of
-- the space one dimension less.
class (LinearMap a v, Orthogonal a n o, Diagonal (v a) (n a)) => Projective a v n o m p | m -> p, p -> m, p -> o, o -> p, p -> n, n -> p, p -> v, v -> p, n -> o, n -> v, v -> n where
  fromProjective     :: p a -> m a
  toProjectiveUnsafe :: m a -> p a
  orthogonal         :: o a -> p a
  linear             :: n a -> p a
  translation        :: v a -> p a
  scaling            :: v a -> p a


class HasOne v where
  _1 :: v a -> a

class HasTwo v where
  _2 :: v a -> a

class HasThree v where
  _3 :: v a -> a

class HasFour v where
  _4 :: v a -> a
