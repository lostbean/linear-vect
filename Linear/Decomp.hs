{-# LANGUAGE
    FlexibleContexts
  , FlexibleInstances
  , TypeSynonymInstances
  , BangPatterns
  , MultiParamTypeClasses
  , ScopedTypeVariables
  #-}
module Linear.Decomp
  ( symmEigen
  , qrGram
  , qrHouse
  , Hessenberg (..)
  , OrthoMatrix (..)
  ) where

import Data.List (foldl')

import Linear.Vect
import Linear.Mat
import Linear.Class

class Hessenberg m where
  -- | Find the tridiagonalization of symmetric matrices or transforms non-symmetric
  -- matrices to a Hessenberg form. Uses Householder transformation.
  hessen :: m -> m

-- | Class of orthogonalization function to calculate the orthonormal matrix of an
-- m-by-m square matrix A with. This can be used for QR decompositon where Q is an
-- orthogonal matrix and R is an upper triangular matrix.
class OrthoMatrix m where
  -- | Calculates the orthogonal matrix using the Householder reflection (or Householder
  -- transformation).
  orthoRowsHouse :: (Transpose m m) => m -> m
  orthoRowsHouse = transpose . orthoColsHouse . transpose

  -- | Implements the modified Gram-Schmidt algorithm to calculate an orthonormal basis of a
  -- matrix. The modified version has more stability in finite precision calculations.
  -- This particular implementation does the calculation over the matrix's rows.
  orthoRowsGram :: m -> m

  -- | Same as @orthoRowsBasis@ but applied to the matrix's columns. It has a default instance.
  orthoColsHouse :: m -> m

  -- | Same as @orthoRowsBasis@ but applied to the matrix's columns. It has a default instance.
  orthoColsGram :: (Transpose m m) => m -> m
  orthoColsGram = transpose . orthoRowsGram . transpose

-- =======================================================================================

instance (Floating a, Ord a) => OrthoMatrix (Mat2 a) where
  orthoColsHouse = transpose . getQ

  orthoRowsGram (Mat2 a1 a2) = Mat2 e1 e2
    where
      e1 = safeNormalize (Vec2 1 0) a1
      e2 = safeNormalize (Vec2 0 1) $ schimi a2 e1

-- =======================================================================================

instance (Floating a, Ord a) => Hessenberg (Mat3 a) where
  hessen m = q .*. m .*. q
    where q = getHH3 m

-- | Find the Householder transformation used for tridiagonalization of symmetric
-- matrices and for transforming non-symmetric matrices to a Hessenberg form.
getHH3 :: forall a . (Floating a, Ord a, Num a) => Mat3 a -> Mat3 a
getHH3 m = let
  x = trimHead (_1m $ transposeSqr m :: Vec3 a)
  a = let k = norm x in if _1 x > 0 then -k else k
  r = 0.5 * (a * a - (_1 x) * a)
  v = Vec2 (a / (2 * r)) 0
  u = extendHeadZero $ v &- (1 / (2 * r)) *& x
  in householder2 u

instance (Floating a, Ord a, Num a) => OrthoMatrix (Mat3 a) where
  orthoColsHouse m = let
    q1 = getQ m
    m1 = (trimHead $ q1 .*. m) :: Mat2 a
    q2 = extendHeadWith 1 $ getQ m1
    in transpose q1 .*. transpose q2

  orthoRowsGram (Mat3 a1 a2 a3) = Mat3 e1 e2 e3
    where
      e1 = safeNormalize (Vec3 1 0 0) a1
      e2 = safeNormalize (Vec3 0 1 0) $ foldl' schimi a2 [e1]
      e3 = safeNormalize (Vec3 0 0 1) $ foldl' schimi a3 [e1, e2]

-- =======================================================================================

instance (Floating a, Ord a) => Hessenberg (Mat4 a) where
  hessen m = let
    q1  = getHH4 m
    a1  = q1 .*. m .*. q1
    a1s = trimHead a1
    q2  = extendHeadWith 1 $ getHH3 a1s
    in q2 .*. a1 .*. q2

getHH4 :: forall a . (Floating a, Ord a) => Mat4 a -> Mat4 a
getHH4 m = let
  x = trimHead (_1m $ transposeSqr m :: Vec4 a)
  a = let k = norm x in if _1 x > 0 then -k else k
  r = 0.5 * (a * a - (_1 x) * a)
  v = Vec3 (a / (2 * r)) 0 0
  u = extendHeadZero $ v &- (1 / ( 2 * r)) *& x
  in householder2 u

instance (Floating a, Ord a) => OrthoMatrix (Mat4 a) where
  orthoColsHouse m = let
    q1 = getQ m
    m1 = (trimHead $ q1 .*. m)  :: Mat3 a

    q2 = getQ m1
    m2 = (trimHead $ q2 .*. m1) :: Mat2 a
    q3 = getQ m2

    q2e = extendHeadWith 1 q2
    q3e = extendHeadWith 1 q3
    in transpose q1 .*. transpose q2e .*. transpose q3e

  orthoRowsGram (Mat4 a1 a2 a3 a4) = Mat4 e1 e2 e3 e4
    where
      e1 = safeNormalize (Vec4 1 0 0 0) a1
      e2 = safeNormalize (Vec4 0 1 0 0) $ schimi a2 e1
      e3 = safeNormalize (Vec4 0 0 1 0) $ foldl' schimi a3 [e1, e2]
      e4 = safeNormalize (Vec4 0 0 0 1) $ foldl' schimi a4 [e1, e2, e3]

-- =======================================================================================

-- | Calculates the eigenpairs (eigenvalues and eigenvector) of a given
-- *symmetric* matrix. It uses only OrthoMatrix decomposition. It can obtain correct
-- eigenvalues only (the eigenvectors should be post processed) if the
-- input matrix is a tridiagonal symmetric matrix or a upper Hessenberg
-- non-symmetric matrix.
symmEigen :: ( OrthoMatrix (m a)
             , MultSemiGroup (m a)
             , Diagonal  (v a) (m a)
             , Dimension (m a)
             , Transpose (m a) (m a)
             , Functor v
             , Functor m
             , Foldable v
             , Foldable m
             , Fractional a
             , Ord a
             , NearZero a
             ) => m a -> (m a, v a)
symmEigen m = eigen q0 (r0 .*. q0) (10 * dim m)
  where
    (q0, r0) = qrHouse m
    eigen !u !a !count
      | count <= 0      = (u, diagVec a)
      | offdiag < limit = (u, diagVec a)
      | otherwise       = eigen (u .*. q) (r .*. q) (count - 1)
      where
        (diago, offdiag) = diagOffDiag a
        limit   = max (epsilon * diago) epsilon
        (q, r)  = qrHouse a

-- =======================================================================================

-- | QR decomposition using Householder method.
qrHouse :: (OrthoMatrix m, MultSemiGroup m, Transpose m m)=> m -> (m, m)
qrHouse m = (q, r)
  where q = orthoColsHouse m
        r = transpose q .*. m

-- | QR decomposition using Gram-Schmidt method. Less unstable the 'qrHouse'.
-- It also doesn't work properly when the matrix in created by an self outer product
-- e.g. "outer v v"
qrGram :: (OrthoMatrix m, MultSemiGroup m, Transpose m m)=> m -> (m, m)
qrGram m = (transpose q, r)
  where q = orthoRowsGram $ transpose m
        r = q .*. m

-- | Householder matrix, see <http://en.wikipedia.org/wiki/Householder_transformation>.
-- In plain words, it is the reflection to the hyperplane orthogonal to the input vector.
-- The input vector is normalized before the Householder transformation.
householder2 :: (Ord a, Fractional a, LinearMap a v, LinearMap a m, DotProd a v, SquareMatrix (m a), Tensor (m a) (v a)) => v a -> m a
householder2 v
  | l > 0     = idmtx &- ((2 / l) *& (outer v v))
  | otherwise = idmtx
  where l = normsqr v

-- ============================ DO NOT EXPORT OrthoMatrix algorithms =====================

-- | Find the orthogonal component of the fisrt vector in relation to
-- the second vector. Used by the Gram-Schimidt algorithm.
schimi :: (Fractional a, Ord a, LinearMap a g, DotProd a g)=> g a -> g a -> g a
schimi a b
  | lenb > 0  = a &- projAonB
  | otherwise = a
  where
    lenb     = b &. b
    projAonB = b &* ((a &. b) / lenb)

-- | Find the Householder transformation matrix. This is an orthogonal matrix.
getQ :: ( Num a
        , Ord a
        , Tensor (m a) (v a)
        , LinearMap a m
        , LinearMap a v
        , Transpose (m a) (m a)
        , SquareMatrix (m a)
        , HasSliceOne m v
        , HasOne v
        , HasE1 v
        , Norm a v
        , AbelianGroup (v a)
        ) => m a -> m a
getQ m = let
  x = _1m $ transposeSqr m
  l = norm x
  a = if _1 x > 0 then -l else l
  u = x &+ vece1 a
  in householder2 u

class HasSliceOne m v where
  _1m :: m a -> v a

instance HasSliceOne Mat2 Vec2 where
  _1m (Mat2 o _) = o

instance HasSliceOne Mat3 Vec3 where
  _1m (Mat3 o _ _) = o

instance HasSliceOne Mat4 Vec4 where
  _1m (Mat4 o _ _ _) = o

transposeSqr :: (Transpose m m) => m -> m
transposeSqr = transpose

-- | The Gershgorin circle is the sum of the absolute values of the non-diagonal entries.
diagOffDiag :: ( Num a
               , Functor m
               , Foldable m
               , Functor v
               , Foldable v
               , Diagonal (v a) (m a)
               ) => m a -> (a, a)
diagOffDiag m = let
  diago = sum $ fmap abs (diagVec m)
  total = sum $ fmap abs m
  in (diago, total - diago)

safeNormalize :: ( LinearMap a v
                 , DotProd a v
                 , Norm a v
                 , Ord a
                 ) => v a -> v a -> v a
safeNormalize fallback vec
  | l > 0     = vec &* (1 / l)
  | otherwise = fallback
  where l = norm vec

class HasE1 v where
  vece1 :: (Floating a, Num a) => a -> v a
instance HasE1 Vec2 where
  vece1 x = Vec2 x 0
instance HasE1 Vec3 where
  vece1 x = Vec3 x 0 0
instance HasE1 Vec4 where
  vece1 x = Vec4 x 0 0 0

-- ================================== Test ==================================

testData1 :: Mat3 Double -- Source <http://en.wikipedia.org/wiki/QR_decomposition>
testData1 = Mat3 (Vec3 12 (-51) 4) (Vec3 6 167 (-68)) (Vec3 (-4) 24 (-41))

testData2 :: Mat3 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData2 = Mat3 (Vec3 2 1 0) (Vec3 1 3 (-1)) (Vec3 0 (-1) 6)

testData3 :: Mat4 Double
testData3 = Mat4 (Vec4 0 10 3 9) (Vec4 10 12 6 15) (Vec4 3 6 0 7) (Vec4 9 15 7 8)

testData4 :: Mat4 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData4 = Mat4 (Vec4 4 1 (-1) 2) (Vec4 1 4 1 (-1)) (Vec4 (-1) 1 4 1) (Vec4 2 (-1) 1 4)

testData5 :: Mat2 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData5 = Mat2 (Vec2 2 1) (Vec2 1 3)

testQR :: (Transpose g g, MultSemiGroup g, AbelianGroup g, OrthoMatrix g, SquareMatrix g) => g -> g
testQR m = m &- (q .*. r)
  where (q, r) = qrHouse m

testEigen m = map (normsqr . foo) $ take n [(_1, _1), (_2, _2), (_3, _3), (_4, _4)]
  where
    n = dim m
    foo (f1, f2) = (m &- (f1 value) *& idmtx) *. (f2 $ transpose vec)
    (vec, value) = symmEigen m
