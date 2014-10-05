{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}

module Linear.Vect
  ( V2(..), V3(..), V4(..)
  , Normal2, Normal3, Normal4
  , mkV2, mkV3, mkV4
  , HasV2, HasV3, HasV4
  , _x, _y, _z, _w
  )
where

import Data.Bits
import Foreign.Storable
import Foreign.Ptr
import System.Random
import Linear.Class

--------------------------------------------------------------------------------
-- Vec datatypes

data V2 a = V2 !a !a 
  deriving (Read,Show)
data V3 a = V3 !a !a !a
  deriving (Read,Show)
data V4 a = V4 !a !a !a !a
  deriving (Read,Show)

-- | The assumption when dealing with these is always that they are of unit length.
-- Also, interpolation works differently.
newtype Normal2 a = Normal2 (V2 a) deriving (Read,Show,Storable,Dimension) 
newtype Normal3 a = Normal3 (V3 a) deriving (Read,Show,Storable,Dimension)
newtype Normal4 a = Normal4 (V4 a) deriving (Read,Show,Storable,Dimension)

deriving instance Floating a => DotProd a Normal2
deriving instance Floating a => DotProd a Normal3
deriving instance Floating a => DotProd a Normal4
deriving instance Floating a => Norm a Normal2
deriving instance Floating a => Norm a Normal3
deriving instance Floating a => Norm a Normal4

mkV2 :: (a,a) -> V2 a
mkV3 :: (a,a,a) -> V3 a
mkV4 :: (a,a,a,a) -> V4 a

mkV2 (x,y)     = V2 x y
mkV3 (x,y,z)   = V3 x y z
mkV4 (x,y,z,w) = V4 x y z w

class HasV2 v where
  getV2 :: v a -> V2 a

_x :: V2 a -> a
_x (V2 x _) = x
_y :: V2 a -> a
_y (V2 _ y) = y

class HasV3 v where
  getV3 :: v a -> V3 a

_z :: V3 a -> a
_z (V3 _ _ z) = z

class HasV4 v where
  getV4 :: v a -> V4 a

_w :: V4 a -> a
_w (V4 _ _ _ w) = w

--------------------------------------------------------------------------------
-- Unit vectors

instance Floating a => UnitVector a V2 Normal2 where
  mkNormal v = Normal2 (normalize v)
  fromNormal (Normal2 v) = v 
  toNormalUnsafe = Normal2

instance Floating a => UnitVector a V3 Normal3 where
  mkNormal v = Normal3 (normalize v)
  fromNormal (Normal3 v) = v 
  toNormalUnsafe = Normal3

instance Floating a => UnitVector a V4 Normal4 where
  mkNormal v = Normal4 (normalize v)
  fromNormal (Normal4 v) = v 
  toNormalUnsafe = Normal4

_rndUnit :: (Fractional a, Ord a, RandomGen g, Random (v a), Vector a v, DotProd a v) => g -> (v a,g)
_rndUnit g = 
  if d > 0.0001
    then ( v &* (1.0/d) , h )
    else _rndUnit h
  where
    (v,h) = random g
    d = normsqr v

instance (Floating a, Random a, Ord a) => Random (Normal2 a) where
  random g = let (v,h) = _rndUnit g in (Normal2 v, h)  
  randomR _ = random

instance (Floating a, Random a, Ord a) => Random (Normal3 a) where
  random g = let (v,h) = _rndUnit g in (Normal3 v, h)  
  randomR _ = random

instance (Floating a, Random a, Ord a) => Random (Normal4 a) where
  random g = let (v,h) = _rndUnit g in (Normal4 v, h)  
  randomR _ = random

instance Floating a => CrossProd (Normal3 a) where
  crossprod (Normal3 v) (Normal3 w) = mkNormal (crossprod v w)

--------------------------------------------------------------------------------
-- V2 instances
instance HasV2 V2 where
  getV2 = id

instance Num a => AbelianGroup (V2 a) where
  (&+) (V2 x1 y1) (V2 x2 y2) = V2 (x1+x2) (y1+y2)
  (&-) (V2 x1 y1) (V2 x2 y2) = V2 (x1-x2) (y1-y2)
  neg  (V2 x y)                = V2 (-x) (-y)
  zero = V2 0 0

instance Num a => Vector a V2 where
  scalarMul s (V2 x y) = V2 (s*x) (s*y)
  mapVec    f (V2 x y) = V2 (f x) (f y)

instance Num a => DotProd a V2 where
  (&.) (V2 x1 y1) (V2 x2 y2) = x1*x2 + y1*y2

instance Floating a => Norm a V2

instance Num a => Pointwise (V2 a) where
  pointwise (V2 x1 y1) (V2 x2 y2) = V2 (x1*x2) (y1*y2)

instance Num a => Determinant a (V2 a, V2 a) where
  det (V2 x1 y1 , V2 x2 y2) = x1*y2 - x2*y1

{-
instance Show V2 where
  show (V2 x y) = "( " ++ show x ++ " , " ++ show y ++ " )"
-}

instance (Num a, Random a) => Random (V2 a) where
  random = randomR (V2 (-1) (-1),V2 1 1)
  randomR (V2 a b, V2 c d) gen = 
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (V2 x y, gen2)

instance (Num a, Storable a) => Storable (V2 a) where
  -- 4byte aligned
  sizeOf    _ = (sizeOf (undefined :: a) * 2 + 3) .&. 252
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p k
    return (V2 x y)

  poke q (V2 x y) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p   x
    pokeByteOff p k y

instance Num a => Dimension (V2 a) where dim _ = 2

--------------------------------------------------------------------------------     
-- V3 instances
instance HasV2 V3 where
  getV2 (V3 x y _) = V2 x y

instance HasV3 V3 where
  getV3 = id

instance Num a => AbelianGroup (V3 a) where
  (&+) (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (x1+x2) (y1+y2) (z1+z2) 
  (&-) (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (x1-x2) (y1-y2) (z1-z2) 
  neg  (V3 x y z)                    = V3 (-x) (-y) (-z)
  zero = V3 0 0 0

instance Num a => Vector a V3 where
  scalarMul s (V3 x y z) = V3 (s*x) (s*y) (s*z)
  mapVec    f (V3 x y z) = V3 (f x) (f y) (f z)

instance Num a => DotProd a V3 where
  (&.) (V3 x1 y1 z1) (V3 x2 y2 z2) = x1*x2 + y1*y2 + z1*z2

instance Floating a => Norm a V3

instance Num a => Pointwise (V3 a) where
  pointwise (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (x1*x2) (y1*y2) (z1*z2)

{-
instance Show V3 where
  show (V3 x y z) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " )"
-}

instance (Num a, Random a) => Random (V3 a) where
  random = randomR (V3 (-1) (-1) (-1),V3 1 1 1)
  randomR (V3 a b c, V3 d e f) gen = 
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2  
    in (V3 x y z, gen3)

instance Num a => CrossProd (V3 a) where
  crossprod (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (y1*z2-y2*z1) (z1*x2-z2*x1) (x1*y2-x2*y1) 

instance Num a => Determinant a (V3 a, V3 a, V3 a) where
  det (u,v,w) = u &. (v &^ w)  

instance (Num a, Storable a) => Storable (V3 a) where
  sizeOf    _ = (sizeOf (undefined :: a) * 3 + 3) .&. 252
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    return (V3 x y z)

  poke q (V3 x y z) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z

instance Num a => Dimension (V3 a) where dim _ = 3

--------------------------------------------------------------------------------
-- V4 instances
instance HasV2 V4 where
  getV2 (V4 x y _ _) = V2 x y

instance HasV3 V4 where
  getV3 (V4 x y z w) = V3 x y z

instance HasV4 V4 where
  getV4 = id

instance Num a => AbelianGroup (V4 a) where
  (&+) (V4 x1 y1 z1 w1) (V4 x2 y2 z2 w2) = V4 (x1+x2) (y1+y2) (z1+z2) (w1+w2)
  (&-) (V4 x1 y1 z1 w1) (V4 x2 y2 z2 w2) = V4 (x1-x2) (y1-y2) (z1-z2) (w1-w2)
  neg  (V4 x y z w)                      = V4 (-x) (-y) (-z) (-w)
  zero = V4 0 0 0 0

instance Num a => Vector a V4 where
  scalarMul s (V4 x y z w) = V4 (s*x) (s*y) (s*z) (s*w)
  mapVec    f (V4 x y z w) = V4 (f x) (f y) (f z) (f w)

instance Num a => DotProd a V4 where
  (&.) (V4 x1 y1 z1 w1) (V4 x2 y2 z2 w2) = x1*x2 + y1*y2 + z1*z2 + w1*w2

instance Floating a => Norm a V4

instance Num a => Pointwise (V4 a) where
  pointwise (V4 x1 y1 z1 w1) (V4 x2 y2 z2 w2) = V4 (x1*x2) (y1*y2) (z1*z2) (w1*w2)

{-
instance Show V4 where
  show (V4 x y z w) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " , " ++ show w ++ " )"
-}

instance (Num a, Random a) => Random (V4 a) where
  random = randomR (V4 (-1) (-1) (-1) (-1),V4 1 1 1 1)
  randomR (V4 a b c d, V4 e f g h) gen = 
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2  
        (w,gen4) = randomR (d,h) gen3  
    in (V4 x y z w, gen4)

instance (Num a, Storable a) => Storable (V4 a) where
  sizeOf    _ = 4 * sizeOf (undefined :: a)
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    w <- peekByteOff p (3*k)
    return (V4 x y z w)

  poke q (V4 x y z w) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z
    pokeByteOff p (3*k) w

instance Num a => Dimension (V4 a) where dim _ = 4

--------------------------------------------------------------------------------
-- Extend instances

instance Num a => Extend a V2 V3 where
  extendZero   (V2 x y) = V3 x y 0
  extendWith t (V2 x y) = V3 x y t
  trim (V3 x y _)       = V2 x y

instance Num a => Extend a V2 V4 where
  extendZero   (V2 x y) = V4 x y 0 0
  extendWith t (V2 x y) = V4 x y t t
  trim (V4 x y _ _)     = V2 x y 

instance Num a => Extend a V3 V4 where
  extendZero   (V3 x y z) = V4 x y z 0
  extendWith t (V3 x y z) = V4 x y z t
  trim (V4 x y z _)       = V3 x y z

