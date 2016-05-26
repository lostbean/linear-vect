{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}

module Linear.Mat
  ( Mat2(..) , Mat3(..) , Mat4(..)
  , Mat2x3(..) , Mat2x4(..) , Mat3x2(..)
  , Mat3x4(..) , Mat4x2(..) , Mat4x3(..)
  , Ortho2 , Ortho3 , Ortho4
  , Proj3 , Proj4
  )
where

import Control.DeepSeq
import Foreign.Storable
import Foreign.Ptr
import System.Random

import Linear.Class
import Linear.Vect

--------------------------------------------------------------------------------
-- Mat datatypes

-- | The components are /row/ vectors

data Mat2 a = Mat2 !(Vec2 a) !(Vec2 a) deriving (Read,Show)
data Mat3 a = Mat3 !(Vec3 a) !(Vec3 a) !(Vec3 a) deriving (Read,Show)
data Mat4 a = Mat4 !(Vec4 a) !(Vec4 a) !(Vec4 a) !(Vec4 a) deriving (Read,Show)
data Mat2x3 a = Mat2x3 !a !a !a !a !a !a deriving (Read,Show)
data Mat2x4 a = Mat2x4 !a !a !a !a !a !a !a !a deriving (Read,Show)
data Mat3x2 a = Mat3x2 !a !a !a !a !a !a deriving (Read,Show)
data Mat3x4 a = Mat3x4 !a !a !a !a !a !a !a !a !a !a !a !a deriving (Read,Show)
data Mat4x2 a = Mat4x2 !a !a !a !a !a !a !a !a deriving (Read,Show)
data Mat4x3 a = Mat4x3 !a !a !a !a !a !a !a !a !a !a !a !a deriving (Read,Show)

-- | Orthogonal matrices.
--
-- Note: the "Random" instances generates orthogonal matrices with determinant 1
-- (that is, orientation-preserving orthogonal transformations)!
newtype Ortho2 a = Ortho2 (Mat2 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho3 a = Ortho3 (Mat3 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho4 a = Ortho4 (Mat4 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)

-- | Projective matrices, encoding affine transformations in dimension one less.
newtype Proj3 a = Proj3 (Mat3 a) deriving (Read,Show,Storable,MultSemiGroup)
newtype Proj4 a = Proj4 (Mat4 a) deriving (Read,Show,Storable,MultSemiGroup)

--------------------------------------------------------------------------------
-- Orthogonal matrices

instance Fractional a =>  Orthogonal a Mat2 Ortho2 where
  fromOrtho (Ortho2 o) = o
  toOrthoUnsafe = Ortho2

instance Fractional a => Orthogonal a Mat3 Ortho3 where
  fromOrtho (Ortho3 o) = o
  toOrthoUnsafe = Ortho3

instance Fractional a => Orthogonal a Mat4 Ortho4 where
  fromOrtho (Ortho4 o) = o
  toOrthoUnsafe = Ortho4

------

instance Transpose (Ortho2 a) (Ortho2 a) where
  transpose (Ortho2 o) = Ortho2 (transpose o)

instance Fractional a => SquareMatrix (Ortho2 a) where
  idmtx = Ortho2 idmtx
  inverse = transpose

instance Transpose (Ortho3 a) (Ortho3 a) where
  transpose (Ortho3 o) = Ortho3 (transpose o)

instance Fractional a => SquareMatrix (Ortho3 a) where
  idmtx = Ortho3 idmtx
  inverse = transpose

instance Transpose (Ortho4 a) (Ortho4 a) where
  transpose (Ortho4 o) = Ortho4 (transpose o)

instance Fractional a => SquareMatrix (Ortho4 a) where
  idmtx = Ortho4 idmtx
  inverse = transpose

------

instance (Floating a, Ord a, Random a) => Random (Ortho2 a) where
  random g = let (o,h) = _rndOrtho2 g in (toOrthoUnsafe (_flip1stRow2 o), h)
  randomR _ = random

instance (Floating a, Ord a, Random a) => Random (Ortho3 a) where
  -- why?
  random g = let (o,h) = _rndOrtho3 g in (toOrthoUnsafe (             o), h)
  randomR _ = random

instance (Floating a, Ord a, Random a) => Random (Ortho4 a) where
  random g = let (o,h) = _rndOrtho4 g in (toOrthoUnsafe (_flip1stRow4 o), h)
  randomR _ = random

------

-- determinant will be -1
_rndOrtho2 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat2 a, g)
_rndOrtho2 g = (h2, g1) where
  h2 = householder u2
  (u2,g1) = random g

-- generates a uniformly random orthogonal 3x3 matrix
-- /with determinant +1/, with respect to the Haar measure of SO3.
--
-- see Theorem 4 in:
-- Francesco Mezzadri: How to Generate Random Matrices from the Classical Compact Groups
-- Notices of the AMS, May 2007 issue
-- <http://www.ams.org/notices/200705/fea-mezzadri-web.ps>
_rndOrtho3 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat3 a, g)
_rndOrtho3 g = ( (h3 .*. m3), g2) where
  m3 = (extendWith :: Num a => a -> Mat2 a -> Mat3 a) 1 o2
  h3 = householder u3
  (u3,g1) = random g
  (o2,g2) = _rndOrtho2 g1

-- determinant will be -1
_rndOrtho4 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat4 a, g)
_rndOrtho4 g = ( (h4 .*. m4), g2) where
  m4 = (extendWith :: Num a => a -> Mat3 a -> Mat4 a) 1 o3
  h4 = householder u4
  (u4,g1) = random g
  (o3,g2) = _rndOrtho3 g1

------

_flip1stRow2 :: Num a => Mat2 a -> Mat2 a
_flip1stRow2 (Mat2 a b) = Mat2 (neg a) b

_flip1stRow3 :: Num a => Mat3 a -> Mat3 a
_flip1stRow3 (Mat3 a b c) = Mat3 (neg a) b c

_flip1stRow4 :: Num a => Mat4 a -> Mat4 a
_flip1stRow4 (Mat4 a b c d) = Mat4 (neg a) b c d

--------------------------------------------------------------------------------
-- projective matrices

instance Fractional a => Projective a Vec2 Mat2 Ortho2 Mat3 Proj3 where
  fromProjective (Proj3 m) = m
  toProjectiveUnsafe = Proj3
  orthogonal = Proj3 . extendWith (1 :: a) . fromOrtho
  linear     = Proj3 . extendWith (1 :: a)
  translation v = Proj3 $ Mat3 (Vec3 1 0 0) (Vec3 0 1 0) (extendWith (1 :: a) v)
  scaling     v = Proj3 $ diag (extendWith (1 :: a) v)

instance Fractional a => Projective a Vec3 Mat3 Ortho3 Mat4 Proj4 where
  fromProjective (Proj4 m) = m
  toProjectiveUnsafe = Proj4
  orthogonal = Proj4 . extendWith (1 :: a) . fromOrtho
  linear     = Proj4 . extendWith (1 :: a)
  translation v = Proj4 $ Mat4 (Vec4 1 0 0 0) (Vec4 0 1 0 0) (Vec4 0 0 1 0) (extendWith (1 :: a) v)
  scaling     v = Proj4 $ diag (extendWith (1 :: a) v)

instance Transpose (Proj3 a) (Proj3 a) where
  transpose (Proj3 m) = Proj3 (transpose m)

instance Fractional a => SquareMatrix (Proj3 a) where
  idmtx = Proj3 idmtx
  inverse = _invertProj3

instance Transpose (Proj4 a) (Proj4 a) where
  transpose (Proj4 m) = Proj4 (transpose m)

instance Fractional a => SquareMatrix (Proj4 a) where
  idmtx = Proj4 idmtx
  inverse = _invertProj4

_invertProj3 :: (Fractional a, Extend a Vec2 Vec3) => Proj3 a -> Proj3 a
_invertProj3 (Proj3 mat@(Mat3 _ _ t)) =
  Proj3 $ Mat3 (extendZero a) (extendZero b) (extendWith 1 t')
  where
    t' = neg $ (trim :: Extend a Vec2 Vec3 => Vec3 a -> Vec2 a) t .* invm2
    invm2@(Mat2 a b) = inverse $ trim mat

-- Inverts a projective 4x4 matrix. But you can simply use "inverse" instead.
-- We assume that the bottom-right corner is 1.
_invertProj4 :: Fractional a => Proj4 a -> Proj4 a
_invertProj4 (Proj4 mat@(Mat4 _ _ _ t)) =
  Proj4 $ Mat4 (extendZero a) (extendZero b) (extendZero c) (extendWith 1 t')
  where
    t' = neg $ (trim :: Extend a Vec3 Vec4 => Vec4 a -> Vec3 a) t .* invm3
    invm3@(Mat3 a b c) = inverse $ trim mat


--------------------------------------------------------------------------------
-- Mat2 instances

instance Transpose (Mat2 a) (Mat2 a) where
  transpose (Mat2 (Vec2 x1 y1) (Vec2 x2 y2)) =
    Mat2 (Vec2 x1 x2)
       (Vec2 y1 y2)

instance Fractional a => SquareMatrix (Mat2 a) where
  idmtx = Mat2 (Vec2 1 0) (Vec2 0 1)
  inverse (Mat2 (Vec2 a b) (Vec2 c d)) =
    Mat2 (Vec2 (d*r) (-b*r)) (Vec2 (-c*r) (a*r))
    where r = 1.0 / (a*d - b*c)

instance Num a => AbelianGroup (Mat2 a) where
  (&+) (Mat2 r1 r2) (Mat2 s1 s2) = Mat2 (r1 &+ s1) (r2 &+ s2)
  (&-) (Mat2 r1 r2) (Mat2 s1 s2) = Mat2 (r1 &- s1) (r2 &- s2)
  neg  (Mat2 r1 r2)              = Mat2 (neg r1) (neg r2)
  zero = Mat2 zero zero

instance Num a => LinearMap a Mat2 where
  scalarMul s (Mat2 r1 r2) = Mat2 (g r1) (g r2) where g = scalarMul s
  mapVec    f (Mat2 r1 r2) = Mat2 (g r1) (g r2) where g = mapVec f

instance Fractional a => MultSemiGroup (Mat2 a) where
  (.*.) (Mat2 r1 r2) n =
    let (Mat2 c1 c2) = transpose n
    in Mat2 (Vec2 (r1 &. c1) (r1 &. c2))
            (Vec2 (r2 &. c1) (r2 &. c2))
  one = idmtx

instance Fractional a => Ring (Mat2 a)

instance Num a => LeftModule (Mat2 a) (Vec2 a) where
  lmul (Mat2 row1 row2) v = Vec2 (row1 &. v) (row2 &. v)

instance Fractional a => RightModule (Vec2 a) (Mat2 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (Vec2 a) (Mat2 a) where
  diag    (Vec2 x y) = Mat2 (Vec2 x 0) (Vec2 0 y)
  diagVec (Mat2 (Vec2 x _) (Vec2 _ y)) = (Vec2 x y)

instance Num a => Tensor (Mat2 a) (Vec2 a) where
  outer (Vec2 a b) (Vec2 x y) = Mat2
    (Vec2 (a*x) (a*y))
    (Vec2 (b*x) (b*y))

instance Num a => Determinant a (Mat2 a) where
  det (Mat2 (Vec2 a b) (Vec2 c d)) = a*d - b*c

{-
instance Show Mat2 where
  show (Mat2 r1 r2) = show r1 ++ "\n" ++ show r2
-}

instance (Num a, Storable a) => Storable (Mat2 a) where
  sizeOf    _ = 2 * sizeOf (undefined :: Vec2 a)
  alignment _ = alignment  (undefined :: Vec2 a)

  peek q = do
    let p = castPtr q :: Ptr (Vec2 a)
        k = sizeOf (undefined :: Vec2 a)
    r1 <- peek        p
    r2 <- peekByteOff p k
    return (Mat2 r1 r2)

  poke q (Mat2 r1 r2) = do
    let p = castPtr q :: Ptr (Vec2 a)
        k = sizeOf (undefined :: Vec2 a)
    poke        p   r1
    pokeByteOff p k r2

instance (Num a, Random a) => Random (Mat2 a) where
  random = randomR (Mat2 v1 v1 , Mat2 v2 v2) where
    v1 = Vec2 (-1) (-1)
    v2 = Vec2   1    1
  randomR (Mat2 a b, Mat2 c d) gen =
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (Mat2 x y, gen2)

instance Num a => Dimension (Mat2 a) where dim _ = 2

instance (Floating a) => MatrixNorms a (Mat2 a) where
  frobeniusNorm (Mat2 r1 r2) =
    sqrt $
      normsqr r1 +
      normsqr r2

instance Num a => Pointwise (Mat2 a) where
  pointwise (Mat2 x1 y1) (Mat2 x2 y2) = Mat2 (x1 &! x2) (y1 &! y2)

instance Functor Mat2 where
  fmap f (Mat2 a b) = Mat2 (fmap f a) (fmap f b)

instance Foldable Mat2 where
    foldr f acc (Mat2 a b) = foldr f (foldr f acc b) a
    foldl f acc (Mat2 a b) = foldl f (foldl f acc a) b
    foldr1 f    (Mat2 a b) = foldr1 f a `f` foldr1 f b
    foldl1 f    (Mat2 a b) = foldl1 f b `f` foldl1 f a
    null                   = const False
    length                 = const 2
    minimum = foldl1 min
    maximum = foldl1 max
    sum     (Mat2 a b) = sum     a * sum     b
    product (Mat2 a b) = product a * product b

--------------------------------------------------------------------------------
-- Mat3 instances

instance Transpose (Mat3 a) (Mat3 a) where
  transpose (Mat3 (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) (Vec3 x3 y3 z3)) =
    Mat3 (Vec3 x1 x2 x3) (Vec3 y1 y2 y3) (Vec3 z1 z2 z3)

instance Fractional a => SquareMatrix (Mat3 a) where
  idmtx = Mat3 (Vec3 1 0 0) (Vec3 0 1 0) (Vec3 0 0 1)

  inverse (Mat3 (Vec3 a b c) (Vec3 e f g) (Vec3 i j k)) =
    Mat3 (Vec3 (d11*r) (d21*r) (d31*r))
         (Vec3 (d12*r) (d22*r) (d32*r))
         (Vec3 (d13*r) (d23*r) (d33*r))
    where
      r = 1.0 / ( a*d11 + b*d12 + c*d13 )

      d11 = f*k - g*j
      d12 = g*i - e*k
      d13 = e*j - f*i

      d31 = b*g - c*f
      d32 = c*e - a*g
      d33 = a*f - b*e

      d21 = c*j - b*k
      d22 = a*k - c*i
      d23 = b*i - a*j

instance Num a => AbelianGroup (Mat3 a) where
  (&+) (Mat3 r1 r2 r3) (Mat3 s1 s2 s3) = Mat3 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3)
  (&-) (Mat3 r1 r2 r3) (Mat3 s1 s2 s3) = Mat3 (r1 &- s1) (r2 &- s2) (r3 &- s3)
  neg  (Mat3 r1 r2 r3)                 = Mat3 (neg r1) (neg r2) (neg r3)
  zero = Mat3 zero zero zero

instance Num a => LinearMap a Mat3 where
  scalarMul s (Mat3 r1 r2 r3) = Mat3 (g r1) (g r2) (g r3) where g = scalarMul s
  mapVec    f (Mat3 r1 r2 r3) = Mat3 (g r1) (g r2) (g r3) where g = mapVec f

instance Fractional a => MultSemiGroup (Mat3 a) where
  (.*.) (Mat3 r1 r2 r3) n =
    let (Mat3 c1 c2 c3) = transpose n
    in Mat3 (Vec3 (r1 &. c1) (r1 &. c2) (r1 &. c3))
            (Vec3 (r2 &. c1) (r2 &. c2) (r2 &. c3))
            (Vec3 (r3 &. c1) (r3 &. c2) (r3 &. c3))
  one = idmtx

instance Fractional a => Ring (Mat3 a)

instance Num a => LeftModule (Mat3 a) (Vec3 a) where
  lmul (Mat3 row1 row2 row3) v = Vec3 (row1 &. v) (row2 &. v) (row3 &. v)

instance Fractional a => RightModule (Vec3 a) (Mat3 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (Vec3 a) (Mat3 a) where
  diag    (Vec3 x y z) = Mat3 (Vec3 x 0 0) (Vec3 0 y 0) (Vec3 0 0 z)
  diagVec (Mat3 (Vec3 x _ _) (Vec3 _ y _) (Vec3 _ _ z)) = (Vec3 x y z)

instance Num a => Tensor (Mat3 a) (Vec3 a) where
  outer (Vec3 a b c) (Vec3 x y z) = Mat3
    (Vec3 (a*x) (a*y) (a*z))
    (Vec3 (b*x) (b*y) (b*z))
    (Vec3 (c*x) (c*y) (c*z))

instance Num a => Determinant a (Mat3 a) where
  det (Mat3 r1 r2 r3) = det (r1,r2,r3)

{-
instance Show Mat3 where
  show (Mat3 r1 r2 r3) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3
-}

instance (Num a, Storable a) => Storable (Mat3 a) where
  sizeOf    _ = 3 * sizeOf (undefined::Vec3 a)
  alignment _ = alignment  (undefined::Vec3 a)

  peek q = do
    let p = castPtr q :: Ptr (Vec3 a)
        k = sizeOf (undefined::Vec3 a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    return (Mat3 r1 r2 r3)

  poke q (Mat3 r1 r2 r3) = do
    let p = castPtr q :: Ptr (Vec3 a)
        k = sizeOf (undefined::Vec3 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3

instance (Num a, Random a) => Random (Mat3 a) where
  random = randomR (Mat3 v1 v1 v1 , Mat3 v2 v2 v2) where
    v1 = Vec3 (-1) (-1) (-1)
    v2 = Vec3   1    1    1
  randomR (Mat3 a b c, Mat3 d e f) gen =
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2
    in (Mat3 x y z, gen3)

instance Num a => Dimension (Mat3 a) where dim _ = 3

instance Floating a => MatrixNorms a (Mat3 a) where
  frobeniusNorm (Mat3 r1 r2 r3)  =
    sqrt $
      normsqr r1 +
      normsqr r2 +
      normsqr r3

instance Num a => Pointwise (Mat3 a) where
  pointwise (Mat3 x1 y1 z1) (Mat3 x2 y2 z2) = Mat3 (x1 &! x2) (y1 &! y2) (z1 &! z2)

instance Functor Mat3 where
  fmap f (Mat3 a b c) = Mat3 (fmap f a) (fmap f b) (fmap f c)

instance Foldable Mat3 where
    foldr f acc (Mat3 a b c) = foldr f (foldr f (foldr f acc c) b) a
    foldl f acc (Mat3 a b c) = foldl f (foldl f (foldl f acc a) b) c
    foldr1 f    (Mat3 a b c) = foldr1 f a `f` foldr1 f b `f` foldr1 f c
    foldl1 f    (Mat3 a b c) = foldl1 f c `f` foldl1 f b `f` foldl1 f a
    null                     = const False
    length                   = const 3
    minimum = foldl1 min
    maximum = foldl1 max
    sum     (Mat3 a b c) = sum     a * sum     b * sum     c
    product (Mat3 a b c) = product a * product b * product c

--------------------------------------------------------------------------------
-- Mat4 instances

instance Transpose (Mat4 a) (Mat4 a) where
  transpose (Mat4 (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) (Vec4 x3 y3 z3 w3) (Vec4 x4 y4 z4 w4)) =
    Mat4 (Vec4 x1 x2 x3 x4)
       (Vec4 y1 y2 y3 y4)
       (Vec4 z1 z2 z3 z4)
       (Vec4 w1 w2 w3 w4)

instance Num a => SquareMatrix (Mat4 a) where
  idmtx = Mat4 (Vec4 1 0 0 0) (Vec4 0 1 0 0) (Vec4 0 0 1 0) (Vec4 0 0 0 1)
  inverse = error "inverse/Mat4: not implemented yet"

instance Num a => AbelianGroup (Mat4 a) where
  (&+) (Mat4 r1 r2 r3 r4) (Mat4 s1 s2 s3 s4) = Mat4 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3) (r4 &+ s4)
  (&-) (Mat4 r1 r2 r3 r4) (Mat4 s1 s2 s3 s4) = Mat4 (r1 &- s1) (r2 &- s2) (r3 &- s3) (r4 &- s4)
  neg  (Mat4 r1 r2 r3 r4)                  = Mat4 (neg r1) (neg r2) (neg r3) (neg r4)
  zero = Mat4 zero zero zero zero

instance Num a => LinearMap a Mat4 where
  scalarMul s (Mat4 r1 r2 r3 r4) = Mat4 (g r1) (g r2) (g r3) (g r4) where g = scalarMul s
  mapVec    f (Mat4 r1 r2 r3 r4) = Mat4 (g r1) (g r2) (g r3) (g r4) where g = mapVec f

instance Num a => MultSemiGroup (Mat4 a) where
  (.*.) (Mat4 r1 r2 r3 r4) n =
    let (Mat4 c1 c2 c3 c4) = transpose n
    in Mat4 (Vec4 (r1 &. c1) (r1 &. c2) (r1 &. c3) (r1 &. c4))
          (Vec4 (r2 &. c1) (r2 &. c2) (r2 &. c3) (r2 &. c4))
          (Vec4 (r3 &. c1) (r3 &. c2) (r3 &. c3) (r3 &. c4))
          (Vec4 (r4 &. c1) (r4 &. c2) (r4 &. c3) (r4 &. c4))
  one = idmtx

instance Num a => Ring (Mat4 a)

instance Num a => LeftModule (Mat4 a) (Vec4 a) where
  lmul (Mat4 row1 row2 row3 row4) v = Vec4 (row1 &. v) (row2 &. v) (row3 &. v) (row4 &. v)

instance Num a => RightModule (Vec4 a) (Mat4 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (Vec4 a) (Mat4 a) where
  diag    (Vec4 x y z w) = Mat4 (Vec4 x 0 0 0) (Vec4 0 y 0 0) (Vec4 0 0 z 0) (Vec4 0 0 0 w)
  diagVec (Mat4 (Vec4 x _ _ _) (Vec4 _ y _ _) (Vec4 _ _ z _) (Vec4 _ _ _ w)) = Vec4 x y z w

instance Num a => Tensor (Mat4 a) (Vec4 a) where
  outer (Vec4 a b c d) (Vec4 x y z w) = Mat4
    (Vec4 (a*x) (a*y) (a*z) (a*w))
    (Vec4 (b*x) (b*y) (b*z) (b*w))
    (Vec4 (c*x) (c*y) (c*z) (c*w))
    (Vec4 (d*x) (d*y) (d*z) (d*w))

instance Num a => Determinant a (Mat4 a) where
  det = error "det/Mat4: not implemented yet"
  -- det (Mat4 r1 r2 r3 r4) =

{-
instance Show Mat4 where
  show (Mat4 r1 r2 r3 r4) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3 ++ "\n" ++ show r4
-}

instance (Num a, Storable a) => Storable (Mat4 a) where
  sizeOf    _ = 4 * sizeOf (undefined::Vec4 a)
  alignment _ = alignment  (undefined::Vec4 a)

  peek q = do
    let p = castPtr q :: Ptr (Vec4 a)
        k = sizeOf (undefined :: Vec4 a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    return (Mat4 r1 r2 r3 r4)

  poke q (Mat4 r1 r2 r3 r4) = do
    let p = castPtr q :: Ptr (Vec4 a)
        k = sizeOf (undefined :: Vec4 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4

instance (Num a, Random a) => Random (Mat4 a) where
  random = randomR (Mat4 v1 v1 v1 v1, Mat4 v2 v2 v2 v2) where
    v1 = Vec4 (-1) (-1) (-1) (-1)
    v2 = Vec4   1    1    1    1
  randomR (Mat4 a b c d, Mat4 e f g h) gen =
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2
        (w,gen4) = randomR (d,h) gen3
    in (Mat4 x y z w, gen4)

instance Num a => Dimension (Mat4 a) where dim _ = 4

instance Floating a => MatrixNorms a (Mat4 a) where
  frobeniusNorm (Mat4 r1 r2 r3 r4) =
    sqrt $
      normsqr r1 +
      normsqr r2 +
      normsqr r3 +
      normsqr r4

instance Num a => Pointwise (Mat4 a) where
  pointwise (Mat4 x1 y1 z1 w1) (Mat4 x2 y2 z2 w2) = Mat4 (x1 &! x2) (y1 &! y2) (z1 &! z2) (w1 &! w2)

instance Functor Mat4 where
  fmap f (Mat4 a b c d) = Mat4 (fmap f a) (fmap f b) (fmap f c) (fmap f d)

instance Foldable Mat4 where
    foldr f acc (Mat4 a b c d) = foldr f (foldr f (foldr f (foldr f acc d) c) b) a
    foldl f acc (Mat4 a b c d) = foldl f (foldl f (foldl f (foldl f acc a) b) c) d
    foldr1 f    (Mat4 a b c d) = foldr1 f a `f` foldr1 f b `f` foldr1 f c `f` foldr1 f d
    foldl1 f    (Mat4 a b c d) = foldl1 f d `f` foldl1 f c `f` foldl1 f b `f` foldl1 f a
    null                       = const False
    length                     = const 4
    minimum = foldl1 min
    maximum = foldl1 max
    sum     (Mat4 a b c d) = sum     a * sum     b * sum     c * sum     d
    product (Mat4 a b c d) = product a * product b * product c * product d

--------------------------------------------------------------------------------
-- Extend instances

instance Num a => Extend a Mat2 Mat3 where
  extendHeadZero   (Mat2 p q) = Mat3 zero (extendHeadZero p) (extendHeadZero q)
  extendHeadWith w (Mat2 p q) = Mat3 (Vec3 w 0 0) (extendHeadZero p) (extendHeadZero q)
  trimHead       (Mat3 _ p q) = Mat2 (trimHead p) (trimHead q)
  extendZero       (Mat2 p q) = Mat3 (extendZero p) (extendZero q) zero
  extendWith     w (Mat2 p q) = Mat3 (extendZero p) (extendZero q) (Vec3 0 0 w)
  trim           (Mat3 p q _) = Mat2 (trim p) (trim q)

instance Num a => Extend a Mat2 Mat4 where
  extendHeadZero   (Mat2 p q) = Mat4 zero zero (extendHeadZero p) (extendHeadZero q)
  extendHeadWith w (Mat2 p q) = Mat4 (Vec4 w 0 0 0) (Vec4 0 w 0 0) (extendHeadZero p) (extendHeadZero q)
  trimHead     (Mat4 _ _ p q) = Mat2 (trimHead p) (trimHead q)
  extendZero       (Mat2 p q) = Mat4 (extendZero p) (extendZero q) zero zero
  extendWith     w (Mat2 p q) = Mat4 (extendZero p) (extendZero q) (Vec4 0 0 w 0) (Vec4 0 0 0 w)
  trim         (Mat4 p q _ _) = Mat2 (trim p) (trim q)

instance Num a => Extend a Mat3 Mat4 where
  extendHeadZero   (Mat3 p q r) = Mat4 zero (extendHeadZero p) (extendHeadZero q) (extendHeadZero r)
  extendHeadWith w (Mat3 p q r) = Mat4 (Vec4 w 0 0 0) (extendHeadZero p) (extendHeadZero q) (extendHeadZero r)
  trimHead       (Mat4 _ p q r) = Mat3 (trimHead p) (trimHead q) (trimHead r)
  extendZero       (Mat3 p q r) = Mat4 (extendZero p) (extendZero q) (extendZero r) zero
  extendWith     w (Mat3 p q r) = Mat4 (extendZero p) (extendZero q) (extendZero r) (Vec4 0 0 0 w)
  trim           (Mat4 p q r _) = Mat3 (trim p) (trim q) (trim r)

--------------------------------------------------------------------------------
-- M_x_ instances

instance Transpose (Mat2x3 a) (Mat3x2 a) where
  transpose (Mat2x3 a b c  d e f) =
    Mat3x2 a d b e c f

instance Transpose (Mat2x4 a) (Mat4x2 a) where
  transpose (Mat2x4 a b c d  e f g h) =
    Mat4x2 a e  b f  c g  d h

instance Transpose (Mat3x2 a) (Mat2x3 a) where
  transpose (Mat3x2 a b  c d  e f) =
    Mat2x3 a c e  b d f

instance Transpose (Mat3x4 a) (Mat4x3 a) where
  transpose (Mat3x4 a b c d  e f g h  i j k l) =
    Mat4x3 a e i  b f j  c g k  d h l

instance Transpose (Mat4x2 a) (Mat2x4 a) where
  transpose (Mat4x2 a b  c d  e f  g h) =
    Mat2x4 a c e f  b d f h

instance Transpose (Mat4x3 a) (Mat3x4 a) where
  transpose (Mat4x3 a b c  d e f  g h i  j k l) =
    Mat3x4 a d g j  b e h k  c f i l

instance Storable a => Storable (Mat2x3 a) where
  sizeOf    _ = 6 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    return (Mat2x3 r1 r2 r3 r4 r5 r6)

  poke q (Mat2x3 r1 r2 r3 r4 r5 r6) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6

instance Storable a => Storable (Mat2x4 a) where
  sizeOf    _ = 8 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    r7 <- peekByteOff p (6*k)
    r8 <- peekByteOff p (7*k)
    return (Mat2x4 r1 r2 r3 r4 r5 r6 r7 r8)

  poke q (Mat2x4 r1 r2 r3 r4 r5 r6 r7 r8) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6
    pokeByteOff p (6*k) r7
    pokeByteOff p (7*k) r8

instance Storable a => Storable (Mat3x2 a) where
  sizeOf    _ = 6 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    return (Mat3x2 r1 r2 r3 r4 r5 r6)

  poke q (Mat3x2 r1 r2 r3 r4 r5 r6) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6

instance Storable a => Storable (Mat3x4 a) where
  sizeOf    _ = 12 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    r7 <- peekByteOff p (6*k)
    r8 <- peekByteOff p (7*k)
    r9 <- peekByteOff p (8*k)
    ra <- peekByteOff p (9*k)
    rb <- peekByteOff p (10*k)
    rc <- peekByteOff p (11*k)
    return (Mat3x4 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc)

  poke q (Mat3x4 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6
    pokeByteOff p (6*k) r7
    pokeByteOff p (7*k) r8
    pokeByteOff p (8*k) r9
    pokeByteOff p (9*k) ra
    pokeByteOff p (10*k) rb
    pokeByteOff p (11*k) rc

instance Storable a => Storable (Mat4x2 a) where
  sizeOf    _ = 8 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    r7 <- peekByteOff p (6*k)
    r8 <- peekByteOff p (7*k)
    return (Mat4x2 r1 r2 r3 r4 r5 r6 r7 r8)

  poke q (Mat4x2 r1 r2 r3 r4 r5 r6 r7 r8) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6
    pokeByteOff p (6*k) r7
    pokeByteOff p (7*k) r8

instance Storable a => Storable (Mat4x3 a) where
  sizeOf    _ = 12 * sizeOf (undefined :: a)
  alignment _ = alignment (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    r5 <- peekByteOff p (4*k)
    r6 <- peekByteOff p (5*k)
    r7 <- peekByteOff p (6*k)
    r8 <- peekByteOff p (7*k)
    r9 <- peekByteOff p (8*k)
    ra <- peekByteOff p (9*k)
    rb <- peekByteOff p (10*k)
    rc <- peekByteOff p (11*k)
    return (Mat4x3 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc)

  poke q (Mat4x3 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6
    pokeByteOff p (6*k) r7
    pokeByteOff p (7*k) r8
    pokeByteOff p (8*k) r9
    pokeByteOff p (9*k) ra
    pokeByteOff p (10*k) rb
    pokeByteOff p (11*k) rc

-- --------------------------------- NFData ---------------------------------------------

instance NFData a => NFData (Mat2 a) where
  rnf (Mat2 a b) = a `seq` b `seq` ()

instance NFData a => NFData (Mat3 a) where
  rnf (Mat3 a b c) = a `seq` b `seq` c `seq` ()

instance NFData a => NFData (Mat4 a) where
  rnf (Mat4 a b c d) = a `seq` b `seq` c `seq` d `seq` ()
