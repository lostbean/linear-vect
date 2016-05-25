{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}

module Linear.Mat
  ( M2(..) , M3(..) , M4(..)
  , M2x3(..) , M2x4(..) , M3x2(..)
  , M3x4(..) , M4x2(..) , M4x3(..)
  , Ortho2 , Ortho3 , Ortho4
  , Proj3 , Proj4
  )
where

import Foreign.Storable
import Foreign.Ptr
import System.Random
import Linear.Class
import Linear.Vect

--------------------------------------------------------------------------------
-- Mat datatypes

-- | The components are /row/ vectors

data M2 a = M2 !(V2 a) !(V2 a) deriving (Read,Show)
data M3 a = M3 !(V3 a) !(V3 a) !(V3 a) deriving (Read,Show)
data M4 a = M4 !(V4 a) !(V4 a) !(V4 a) !(V4 a) deriving (Read,Show)
data M2x3 a = M2x3 !a !a !a !a !a !a deriving (Read,Show)
data M2x4 a = M2x4 !a !a !a !a !a !a !a !a deriving (Read,Show)
data M3x2 a = M3x2 !a !a !a !a !a !a deriving (Read,Show)
data M3x4 a = M3x4 !a !a !a !a !a !a !a !a !a !a !a !a deriving (Read,Show)
data M4x2 a = M4x2 !a !a !a !a !a !a !a !a deriving (Read,Show)
data M4x3 a = M4x3 !a !a !a !a !a !a !a !a !a !a !a !a deriving (Read,Show)

-- | Orthogonal matrices.
--
-- Note: the "Random" instances generates orthogonal matrices with determinant 1
-- (that is, orientation-preserving orthogonal transformations)!
newtype Ortho2 a = Ortho2 (M2 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho3 a = Ortho3 (M3 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho4 a = Ortho4 (M4 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)

-- | Projective matrices, encoding affine transformations in dimension one less.
newtype Proj3 a = Proj3 (M3 a) deriving (Read,Show,Storable,MultSemiGroup)
newtype Proj4 a = Proj4 (M4 a) deriving (Read,Show,Storable,MultSemiGroup)

--------------------------------------------------------------------------------
-- Orthogonal matrices

instance Fractional a =>  Orthogonal a M2 Ortho2 where
  fromOrtho (Ortho2 o) = o
  toOrthoUnsafe = Ortho2

instance Fractional a => Orthogonal a M3 Ortho3 where
  fromOrtho (Ortho3 o) = o
  toOrthoUnsafe = Ortho3

instance Fractional a => Orthogonal a M4 Ortho4 where
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
_rndOrtho2 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (M2 a, g)
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
_rndOrtho3 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (M3 a, g)
_rndOrtho3 g = ( (h3 .*. m3), g2) where
  m3 = (extendWith :: Num a => a -> M2 a -> M3 a) 1 o2
  h3 = householder u3
  (u3,g1) = random g
  (o2,g2) = _rndOrtho2 g1

-- determinant will be -1
_rndOrtho4 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (M4 a, g)
_rndOrtho4 g = ( (h4 .*. m4), g2) where
  m4 = (extendWith :: Num a => a -> M3 a -> M4 a) 1 o3
  h4 = householder u4
  (u4,g1) = random g
  (o3,g2) = _rndOrtho3 g1

------

_flip1stRow2 :: Num a => M2 a -> M2 a
_flip1stRow2 (M2 a b) = M2 (neg a) b

_flip1stRow3 :: Num a => M3 a -> M3 a
_flip1stRow3 (M3 a b c) = M3 (neg a) b c

_flip1stRow4 :: Num a => M4 a -> M4 a
_flip1stRow4 (M4 a b c d) = M4 (neg a) b c d

--------------------------------------------------------------------------------
-- projective matrices

instance Fractional a => Projective a V2 M2 Ortho2 M3 Proj3 where
  fromProjective (Proj3 m) = m
  toProjectiveUnsafe = Proj3
  orthogonal = Proj3 . extendWith (1 :: a) . fromOrtho
  linear     = Proj3 . extendWith (1 :: a)
  translation v = Proj3 $ M3 (V3 1 0 0) (V3 0 1 0) (extendWith (1 :: a) v)
  scaling     v = Proj3 $ diag (extendWith (1 :: a) v)

instance Fractional a => Projective a V3 M3 Ortho3 M4 Proj4 where
  fromProjective (Proj4 m) = m
  toProjectiveUnsafe = Proj4
  orthogonal = Proj4 . extendWith (1 :: a) . fromOrtho
  linear     = Proj4 . extendWith (1 :: a)
  translation v = Proj4 $ M4 (V4 1 0 0 0) (V4 0 1 0 0) (V4 0 0 1 0) (extendWith (1 :: a) v)
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

_invertProj3 :: (Fractional a, Extend a V2 V3) => Proj3 a -> Proj3 a
_invertProj3 (Proj3 mat@(M3 _ _ t)) =
  Proj3 $ M3 (extendZero a) (extendZero b) (extendWith 1 t')
  where
    t' = neg $ (trim :: Extend a V2 V3 => V3 a -> V2 a) t .* invm2
    invm2@(M2 a b) = inverse $ trim mat

-- Inverts a projective 4x4 matrix. But you can simply use "inverse" instead.
-- We assume that the bottom-right corner is 1.
_invertProj4 :: Fractional a => Proj4 a -> Proj4 a
_invertProj4 (Proj4 mat@(M4 _ _ _ t)) =
  Proj4 $ M4 (extendZero a) (extendZero b) (extendZero c) (extendWith 1 t')
  where
    t' = neg $ (trim :: Extend a V3 V4 => V4 a -> V3 a) t .* invm3
    invm3@(M3 a b c) = inverse $ trim mat


--------------------------------------------------------------------------------
-- M2 instances

instance Transpose (M2 a) (M2 a) where
  transpose (M2 (V2 x1 y1) (V2 x2 y2)) =
    M2 (V2 x1 x2)
       (V2 y1 y2)

instance Fractional a => SquareMatrix (M2 a) where
  idmtx = M2 (V2 1 0) (V2 0 1)
  inverse (M2 (V2 a b) (V2 c d)) =
    M2 (V2 (d*r) (-b*r)) (V2 (-c*r) (a*r))
    where r = 1.0 / (a*d - b*c)

instance Num a => AbelianGroup (M2 a) where
  (&+) (M2 r1 r2) (M2 s1 s2) = M2 (r1 &+ s1) (r2 &+ s2)
  (&-) (M2 r1 r2) (M2 s1 s2) = M2 (r1 &- s1) (r2 &- s2)
  neg  (M2 r1 r2)              = M2 (neg r1) (neg r2)
  zero = M2 zero zero

instance Num a => Vector a M2 where
  scalarMul s (M2 r1 r2) = M2 (g r1) (g r2) where g = scalarMul s
  mapVec    f (M2 r1 r2) = M2 (g r1) (g r2) where g = mapVec f

instance Fractional a => MultSemiGroup (M2 a) where
  (.*.) (M2 r1 r2) n =
    let (M2 c1 c2) = transpose n
    in M2 (V2 (r1 &. c1) (r1 &. c2))
            (V2 (r2 &. c1) (r2 &. c2))
  one = idmtx

instance Fractional a => Ring (M2 a)

instance Num a => LeftModule (M2 a) (V2 a) where
  lmul (M2 row1 row2) v = V2 (row1 &. v) (row2 &. v)

instance Fractional a => RightModule (V2 a) (M2 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (V2 a) (M2 a) where
  diag (V2 x y) = M2 (V2 x 0) (V2 0 y)

instance Num a => Tensor (M2 a) (V2 a) where
  outer (V2 a b) (V2 x y) = M2
    (V2 (a*x) (a*y))
    (V2 (b*x) (b*y))

instance Num a => Determinant a (M2 a) where
  det (M2 (V2 a b) (V2 c d)) = a*d - b*c

{-
instance Show M2 where
  show (M2 r1 r2) = show r1 ++ "\n" ++ show r2
-}

instance (Num a, Storable a) => Storable (M2 a) where
  sizeOf    _ = 2 * sizeOf (undefined :: V2 a)
  alignment _ = alignment  (undefined :: V2 a)

  peek q = do
    let p = castPtr q :: Ptr (V2 a)
        k = sizeOf (undefined :: V2 a)
    r1 <- peek        p
    r2 <- peekByteOff p k
    return (M2 r1 r2)

  poke q (M2 r1 r2) = do
    let p = castPtr q :: Ptr (V2 a)
        k = sizeOf (undefined :: V2 a)
    poke        p   r1
    pokeByteOff p k r2

instance (Num a, Random a) => Random (M2 a) where
  random = randomR (M2 v1 v1 , M2 v2 v2) where
    v1 = V2 (-1) (-1)
    v2 = V2   1    1
  randomR (M2 a b, M2 c d) gen =
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (M2 x y, gen2)

instance Num a => Dimension (M2 a) where dim _ = 2

instance (Floating a) => MatrixNorms a (M2 a) where
  frobeniusNorm (M2 r1 r2) =
    sqrt $
      normsqr r1 +
      normsqr r2

instance Num a => Pointwise (M2 a) where
  pointwise (M2 x1 y1) (M2 x2 y2) = M2 (x1 &! x2) (y1 &! y2)


--------------------------------------------------------------------------------
-- M3 instances

instance Transpose (M3 a) (M3 a) where
  transpose (M3 (V3 x1 y1 z1) (V3 x2 y2 z2) (V3 x3 y3 z3)) =
    M3 (V3 x1 x2 x3) (V3 y1 y2 y3) (V3 z1 z2 z3)

instance Fractional a => SquareMatrix (M3 a) where
  idmtx = M3 (V3 1 0 0) (V3 0 1 0) (V3 0 0 1)

  inverse (M3 (V3 a b c) (V3 e f g) (V3 i j k)) =
    M3 (V3 (d11*r) (d21*r) (d31*r))
         (V3 (d12*r) (d22*r) (d32*r))
         (V3 (d13*r) (d23*r) (d33*r))
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

instance Num a => AbelianGroup (M3 a) where
  (&+) (M3 r1 r2 r3) (M3 s1 s2 s3) = M3 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3)
  (&-) (M3 r1 r2 r3) (M3 s1 s2 s3) = M3 (r1 &- s1) (r2 &- s2) (r3 &- s3)
  neg  (M3 r1 r2 r3)                 = M3 (neg r1) (neg r2) (neg r3)
  zero = M3 zero zero zero

instance Num a => Vector a M3 where
  scalarMul s (M3 r1 r2 r3) = M3 (g r1) (g r2) (g r3) where g = scalarMul s
  mapVec    f (M3 r1 r2 r3) = M3 (g r1) (g r2) (g r3) where g = mapVec f

instance Fractional a => MultSemiGroup (M3 a) where
  (.*.) (M3 r1 r2 r3) n =
    let (M3 c1 c2 c3) = transpose n
    in M3 (V3 (r1 &. c1) (r1 &. c2) (r1 &. c3))
            (V3 (r2 &. c1) (r2 &. c2) (r2 &. c3))
            (V3 (r3 &. c1) (r3 &. c2) (r3 &. c3))
  one = idmtx

instance Fractional a => Ring (M3 a)

instance Num a => LeftModule (M3 a) (V3 a) where
  lmul (M3 row1 row2 row3) v = V3 (row1 &. v) (row2 &. v) (row3 &. v)

instance Fractional a => RightModule (V3 a) (M3 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (V3 a) (M3 a) where
  diag (V3 x y z) = M3 (V3 x 0 0) (V3 0 y 0) (V3 0 0 z)

instance Num a => Tensor (M3 a) (V3 a) where
  outer (V3 a b c) (V3 x y z) = M3
    (V3 (a*x) (a*y) (a*z))
    (V3 (b*x) (b*y) (b*z))
    (V3 (c*x) (c*y) (c*z))

instance Num a => Determinant a (M3 a) where
  det (M3 r1 r2 r3) = det (r1,r2,r3)

{-
instance Show M3 where
  show (M3 r1 r2 r3) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3
-}

instance (Num a, Storable a) => Storable (M3 a) where
  sizeOf    _ = 3 * sizeOf (undefined::V3 a)
  alignment _ = alignment  (undefined::V3 a)

  peek q = do
    let p = castPtr q :: Ptr (V3 a)
        k = sizeOf (undefined::V3 a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    return (M3 r1 r2 r3)

  poke q (M3 r1 r2 r3) = do
    let p = castPtr q :: Ptr (V3 a)
        k = sizeOf (undefined::V3 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3

instance (Num a, Random a) => Random (M3 a) where
  random = randomR (M3 v1 v1 v1 , M3 v2 v2 v2) where
    v1 = V3 (-1) (-1) (-1)
    v2 = V3   1    1    1
  randomR (M3 a b c, M3 d e f) gen =
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2
    in (M3 x y z, gen3)

instance Num a => Dimension (M3 a) where dim _ = 3

instance Floating a => MatrixNorms a (M3 a) where
  frobeniusNorm (M3 r1 r2 r3)  =
    sqrt $
      normsqr r1 +
      normsqr r2 +
      normsqr r3

instance Num a => Pointwise (M3 a) where
  pointwise (M3 x1 y1 z1) (M3 x2 y2 z2) = M3 (x1 &! x2) (y1 &! y2) (z1 &! z2)

--------------------------------------------------------------------------------
-- M4 instances

instance Transpose (M4 a) (M4 a) where
  transpose (M4 (V4 x1 y1 z1 w1) (V4 x2 y2 z2 w2) (V4 x3 y3 z3 w3) (V4 x4 y4 z4 w4)) =
    M4 (V4 x1 x2 x3 x4)
       (V4 y1 y2 y3 y4)
       (V4 z1 z2 z3 z4)
       (V4 w1 w2 w3 w4)

instance Num a => SquareMatrix (M4 a) where
  idmtx = M4 (V4 1 0 0 0) (V4 0 1 0 0) (V4 0 0 1 0) (V4 0 0 0 1)
  inverse = error "inverse/M4: not implemented yet"

instance Num a => AbelianGroup (M4 a) where
  (&+) (M4 r1 r2 r3 r4) (M4 s1 s2 s3 s4) = M4 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3) (r4 &+ s4)
  (&-) (M4 r1 r2 r3 r4) (M4 s1 s2 s3 s4) = M4 (r1 &- s1) (r2 &- s2) (r3 &- s3) (r4 &- s4)
  neg  (M4 r1 r2 r3 r4)                  = M4 (neg r1) (neg r2) (neg r3) (neg r4)
  zero = M4 zero zero zero zero

instance Num a => Vector a M4 where
  scalarMul s (M4 r1 r2 r3 r4) = M4 (g r1) (g r2) (g r3) (g r4) where g = scalarMul s
  mapVec    f (M4 r1 r2 r3 r4) = M4 (g r1) (g r2) (g r3) (g r4) where g = mapVec f

instance Num a => MultSemiGroup (M4 a) where
  (.*.) (M4 r1 r2 r3 r4) n =
    let (M4 c1 c2 c3 c4) = transpose n
    in M4 (V4 (r1 &. c1) (r1 &. c2) (r1 &. c3) (r1 &. c4))
          (V4 (r2 &. c1) (r2 &. c2) (r2 &. c3) (r2 &. c4))
          (V4 (r3 &. c1) (r3 &. c2) (r3 &. c3) (r3 &. c4))
          (V4 (r4 &. c1) (r4 &. c2) (r4 &. c3) (r4 &. c4))
  one = idmtx

instance Num a => Ring (M4 a)

instance Num a => LeftModule (M4 a) (V4 a) where
  lmul (M4 row1 row2 row3 row4) v = V4 (row1 &. v) (row2 &. v) (row3 &. v) (row4 &. v)

instance Num a => RightModule (V4 a) (M4 a) where
  rmul v mt = lmul (transpose mt) v

instance Num a => Diagonal (V4 a) (M4 a) where
  diag (V4 x y z w) = M4 (V4 x 0 0 0) (V4 0 y 0 0) (V4 0 0 z 0) (V4 0 0 0 w)

instance Num a => Tensor (M4 a) (V4 a) where
  outer (V4 a b c d) (V4 x y z w) = M4
    (V4 (a*x) (a*y) (a*z) (a*w))
    (V4 (b*x) (b*y) (b*z) (b*w))
    (V4 (c*x) (c*y) (c*z) (c*w))
    (V4 (d*x) (d*y) (d*z) (d*w))

instance Num a => Determinant a (M4 a) where
  det = error "det/M4: not implemented yet"
  -- det (M4 r1 r2 r3 r4) =

{-
instance Show M4 where
  show (M4 r1 r2 r3 r4) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3 ++ "\n" ++ show r4
-}

instance (Num a, Storable a) => Storable (M4 a) where
  sizeOf    _ = 4 * sizeOf (undefined::V4 a)
  alignment _ = alignment  (undefined::V4 a)

  peek q = do
    let p = castPtr q :: Ptr (V4 a)
        k = sizeOf (undefined :: V4 a)
    r1 <- peek        p
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    return (M4 r1 r2 r3 r4)

  poke q (M4 r1 r2 r3 r4) = do
    let p = castPtr q :: Ptr (V4 a)
        k = sizeOf (undefined :: V4 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4

instance (Num a, Random a) => Random (M4 a) where
  random = randomR (M4 v1 v1 v1 v1, M4 v2 v2 v2 v2) where
    v1 = V4 (-1) (-1) (-1) (-1)
    v2 = V4   1    1    1    1
  randomR (M4 a b c d, M4 e f g h) gen =
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2
        (w,gen4) = randomR (d,h) gen3
    in (M4 x y z w, gen4)

instance Num a => Dimension (M4 a) where dim _ = 4

instance Floating a => MatrixNorms a (M4 a) where
  frobeniusNorm (M4 r1 r2 r3 r4) =
    sqrt $
      normsqr r1 +
      normsqr r2 +
      normsqr r3 +
      normsqr r4

instance Num a => Pointwise (M4 a) where
  pointwise (M4 x1 y1 z1 w1) (M4 x2 y2 z2 w2) = M4 (x1 &! x2) (y1 &! y2) (z1 &! z2) (w1 &! w2)

--------------------------------------------------------------------------------
-- Extend instances

instance Num a => Extend a M2 M3 where
  extendZero   (M2 p q) = M3 (extendZero p) (extendZero q) zero
  extendWith w (M2 p q) = M3 (extendZero p) (extendZero q) (V3 0 0 w)
  trim (M3 p q _) = M2 (trim p) (trim q)

instance Num a => Extend a M2 M4 where
  extendZero   (M2 p q) = M4 (extendZero p) (extendZero q) zero zero
  extendWith w (M2 p q) = M4 (extendZero p) (extendZero q) (V4 0 0 w 0) (V4 0 0 0 w)
  trim (M4 p q _ _) = M2 (trim p) (trim q)

instance Num a => Extend a M3 M4 where
  extendZero   (M3 p q r) = M4 (extendZero p) (extendZero q) (extendZero r) zero
  extendWith w (M3 p q r) = M4 (extendZero p) (extendZero q) (extendZero r) (V4 0 0 0 w)
  trim (M4 p q r _) = M3 (trim p) (trim q) (trim r)


--------------------------------------------------------------------------------
-- M_x_ instances

instance Transpose (M2x3 a) (M3x2 a) where
  transpose (M2x3 a b c  d e f) =
    M3x2 a d b e c f

instance Transpose (M2x4 a) (M4x2 a) where
  transpose (M2x4 a b c d  e f g h) =
    M4x2 a e  b f  c g  d h

instance Transpose (M3x2 a) (M2x3 a) where
  transpose (M3x2 a b  c d  e f) =
    M2x3 a c e  b d f

instance Transpose (M3x4 a) (M4x3 a) where
  transpose (M3x4 a b c d  e f g h  i j k l) =
    M4x3 a e i  b f j  c g k  d h l

instance Transpose (M4x2 a) (M2x4 a) where
  transpose (M4x2 a b  c d  e f  g h) =
    M2x4 a c e f  b d f h

instance Transpose (M4x3 a) (M3x4 a) where
  transpose (M4x3 a b c  d e f  g h i  j k l) =
    M3x4 a d g j  b e h k  c f i l

instance Storable a => Storable (M2x3 a) where
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
    return (M2x3 r1 r2 r3 r4 r5 r6)

  poke q (M2x3 r1 r2 r3 r4 r5 r6) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6

instance Storable a => Storable (M2x4 a) where
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
    return (M2x4 r1 r2 r3 r4 r5 r6 r7 r8)

  poke q (M2x4 r1 r2 r3 r4 r5 r6 r7 r8) = do
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

instance Storable a => Storable (M3x2 a) where
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
    return (M3x2 r1 r2 r3 r4 r5 r6)

  poke q (M3x2 r1 r2 r3 r4 r5 r6) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4
    pokeByteOff p (4*k) r5
    pokeByteOff p (5*k) r6

instance Storable a => Storable (M3x4 a) where
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
    return (M3x4 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc)

  poke q (M3x4 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc) = do
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

instance Storable a => Storable (M4x2 a) where
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
    return (M4x2 r1 r2 r3 r4 r5 r6 r7 r8)

  poke q (M4x2 r1 r2 r3 r4 r5 r6 r7 r8) = do
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

instance Storable a => Storable (M4x3 a) where
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
    return (M4x3 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc)

  poke q (M4x3 r1 r2 r3 r4 r5 r6 r7 r8 r9 ra rb rc) = do
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
