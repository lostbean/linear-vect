{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE
    DeriveGeneric
  , FlexibleInstances
  , GeneralizedNewtypeDeriving
  , IncoherentInstances
  , MultiParamTypeClasses
  , ScopedTypeVariables
  , StandaloneDeriving
  , TemplateHaskell
  , TypeFamilies
#-}

module Linear.Vect
  ( Vec2(..), Vec3(..), Vec4(..)
  , Normal2, Normal3, Normal4
  , mkVec2, mkVec3, mkVec4
  , HasV2, HasV3, HasV4
  , _x, _y, _z, _w
  , unVec2
  , unVec3
  , unVec4

  -- * Alias
  , Vec2D
  , Vec3D
  , Vec4D
  , Normal2D
  , Normal3D
  , Normal4D

  , module Linear.Class
  ) where

import Codec.Serialise (Serialise)
import Control.DeepSeq
import Data.Bits
import Data.Foldable
import Data.Vector.Unboxed (Unbox)
import Data.Vector.Unboxed.Deriving
import Foreign.Storable
import Foreign.Ptr
import GHC.Generics
import System.Random

import Linear.Class

--------------------------------------------------------------------------------
-- Vec datatypes

data Vec2 a = Vec2 !a !a
  deriving (Read, Show, Eq, Generic)
data Vec3 a = Vec3 !a !a !a
  deriving (Read, Show, Eq, Generic)
data Vec4 a = Vec4 !a !a !a !a
  deriving (Read, Show, Eq, Generic)

-- | The assumption when dealing with these is always that they are of unit length.
-- Also, interpolation works differently.
newtype Normal2 a = Normal2 {unNormal2 :: Vec2 a} deriving (Read, Eq, Show, Generic, Storable, Dimension)
newtype Normal3 a = Normal3 {unNormal3 :: Vec3 a} deriving (Read, Eq, Show, Generic, Storable, Dimension)
newtype Normal4 a = Normal4 {unNormal4 :: Vec4 a} deriving (Read, Eq, Show, Generic, Storable, Dimension)

instance (Serialise a)=> Serialise (Vec2 a)
instance (Serialise a)=> Serialise (Vec3 a)
instance (Serialise a)=> Serialise (Vec4 a)
instance (Serialise a)=> Serialise (Normal2 a)
instance (Serialise a)=> Serialise (Normal3 a)
instance (Serialise a)=> Serialise (Normal4 a)

deriving instance Floating a => DotProd a Normal2
deriving instance Floating a => DotProd a Normal3
deriving instance Floating a => DotProd a Normal4
deriving instance Floating a => Norm a Normal2
deriving instance Floating a => Norm a Normal3
deriving instance Floating a => Norm a Normal4

mkVec2 :: (a,a) -> Vec2 a
mkVec3 :: (a,a,a) -> Vec3 a
mkVec4 :: (a,a,a,a) -> Vec4 a

mkVec2 (x,y)     = Vec2 x y
mkVec3 (x,y,z)   = Vec3 x y z
mkVec4 (x,y,z,w) = Vec4 x y z w

unVec2 :: Vec2 a -> (a, a)
unVec3 :: Vec3 a -> (a, a, a)
unVec4 :: Vec4 a -> (a, a, a, a)

unVec2 (Vec2 x y)     = (x,y)
unVec3 (Vec3 x y z)   = (x,y,z)
unVec4 (Vec4 x y z w) = (x,y,z,w)

type Vec2D = Vec2 Double
type Vec3D = Vec3 Double
type Vec4D = Vec4 Double

type Normal2D = Normal2 Double
type Normal3D = Normal3 Double
type Normal4D = Normal4 Double

class HasV2 v where
  getV2 :: v a -> Vec2 a

_x :: Vec2 a -> a
_x (Vec2 x _) = x
_y :: Vec2 a -> a
_y (Vec2 _ y) = y

class HasV3 v where
  getV3 :: v a -> Vec3 a

_z :: Vec3 a -> a
_z (Vec3 _ _ z) = z

class HasV4 v where
  getV4 :: v a -> Vec4 a

_w :: Vec4 a -> a
_w (Vec4 _ _ _ w) = w

instance HasOne Normal2 where
  _1 (Normal2 (Vec2 x _)) = x

instance HasOne Normal3 where
  _1 (Normal3 (Vec3 x _ _)) = x

instance HasOne Normal4 where
  _1 (Normal4 (Vec4 x _ _ _)) = x

instance Show a => PrettyShow (Vec2 a) where
  showPretty = wrapBars . unwords . map show . toList
instance Show a => PrettyShow (Vec3 a) where
  showPretty = wrapBars . unwords . map show . toList
instance Show a => PrettyShow (Vec4 a) where
  showPretty = wrapBars . unwords . map show . toList

instance PrettyShow Vec2D where
  showPretty = wrapBars . unwords . map showPretty . toList
instance PrettyShow Vec3D where
  showPretty = wrapBars . unwords . map showPretty . toList
instance PrettyShow Vec4D where
  showPretty = wrapBars . unwords . map showPretty . toList

instance Show a => PrettyShow (Normal2 a) where
  showPretty = showPretty . unNormal2
instance Show a => PrettyShow (Normal3 a) where
  showPretty = showPretty . unNormal3
instance Show a => PrettyShow (Normal4 a) where
  showPretty = showPretty . unNormal4

instance PrettyShow Normal2D where
  showPretty = showPretty . unNormal2
instance PrettyShow Normal3D where
  showPretty = showPretty . unNormal3
instance PrettyShow Normal4D where
  showPretty = showPretty . unNormal4

--------------------------------------------------------------------------------
-- Unit vectors

instance Floating a => UnitVector a Vec2 Normal2 where
  mkNormal v = Normal2 (normalize v)
  fromNormal (Normal2 v) = v
  toNormalUnsafe = Normal2

instance Floating a => UnitVector a Vec3 Normal3 where
  mkNormal v = Normal3 (normalize v)
  fromNormal (Normal3 v) = v
  toNormalUnsafe = Normal3

instance Floating a => UnitVector a Vec4 Normal4 where
  mkNormal v = Normal4 (normalize v)
  fromNormal (Normal4 v) = v
  toNormalUnsafe = Normal4

_rndUnit :: (Fractional a, Ord a, RandomGen g, Random (v a), LinearMap a v, DotProd a v) => g -> (v a,g)
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
-- Vec2 instances
instance HasV2 Vec2 where
  getV2 = id

instance Num a => AbelianGroup (Vec2 a) where
  (&+) (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1+x2) (y1+y2)
  (&-) (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1-x2) (y1-y2)
  neg  (Vec2 x y)                = Vec2 (-x) (-y)
  zero = Vec2 0 0

instance Num a => LinearMap a Vec2 where
  scalarMul s (Vec2 x y) = Vec2 (s*x) (s*y)
  mapVec    f (Vec2 x y) = Vec2 (f x) (f y)

instance Num a => DotProd a Vec2 where
  (&.) (Vec2 x1 y1) (Vec2 x2 y2) = x1*x2 + y1*y2

instance Floating a => Norm a Vec2

instance Num a => Pointwise (Vec2 a) where
  pointwise (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1*x2) (y1*y2)

instance Num a => Determinant a (Vec2 a, Vec2 a) where
  det (Vec2 x1 y1 , Vec2 x2 y2) = x1*y2 - x2*y1

{-
instance Show Vec2 where
  show (Vec2 x y) = "( " ++ show x ++ " , " ++ show y ++ " )"
-}

instance (Num a, Random a) => Random (Vec2 a) where
  random = randomR (Vec2 (-1) (-1),Vec2 1 1)
  randomR (Vec2 a b, Vec2 c d) gen =
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (Vec2 x y, gen2)

instance (Num a, Storable a) => Storable (Vec2 a) where
  -- 4byte aligned
  sizeOf    _ = (sizeOf (undefined :: a) * 2 + 3) .&. 252
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p
    y <- peekByteOff p k
    return (Vec2 x y)

  poke q (Vec2 x y) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p   x
    pokeByteOff p k y

instance Num a => Dimension (Vec2 a) where dim _ = 2

instance Functor Vec2 where
  fmap f (Vec2 a b) = Vec2 (f a) (f b)

instance Foldable Vec2 where
    foldr f acc (Vec2 a b) = f a (f b acc)
    foldl f acc (Vec2 a b) = acc `f` b `f` a
    foldr1 f    (Vec2 a b) = a `f` b
    foldl1 f    (Vec2 a b) = b `f` a
    null                   = const False
    length                 = const 2
    maximum     (Vec2 a b) = a `max` b
    minimum     (Vec2 a b) = a `min` b
    sum         (Vec2 a b) = a + b
    product     (Vec2 a b) = a * b

instance HasOne Vec2 where
  _1 (Vec2 x _) = x

instance HasTwo Vec2 where
  _2 (Vec2 _ y) = y

--------------------------------------------------------------------------------
-- Vec3 instances
instance HasV2 Vec3 where
  getV2 (Vec3 x y _) = Vec2 x y

instance HasV3 Vec3 where
  getV3 = id

instance Num a => AbelianGroup (Vec3 a) where
  (&+) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1+x2) (y1+y2) (z1+z2)
  (&-) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1-x2) (y1-y2) (z1-z2)
  neg  (Vec3 x y z)                    = Vec3 (-x) (-y) (-z)
  zero = Vec3 0 0 0

instance Num a => LinearMap a Vec3 where
  scalarMul s (Vec3 x y z) = Vec3 (s*x) (s*y) (s*z)
  mapVec    f (Vec3 x y z) = Vec3 (f x) (f y) (f z)

instance Num a => DotProd a Vec3 where
  (&.) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = x1*x2 + y1*y2 + z1*z2

instance Floating a => Norm a Vec3

instance Num a => Pointwise (Vec3 a) where
  pointwise (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1*x2) (y1*y2) (z1*z2)

{-
instance Show Vec3 where
  show (Vec3 x y z) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " )"
-}

instance (Num a, Random a) => Random (Vec3 a) where
  random = randomR (Vec3 (-1) (-1) (-1),Vec3 1 1 1)
  randomR (Vec3 a b c, Vec3 d e f) gen =
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2
    in (Vec3 x y z, gen3)

instance Num a => CrossProd (Vec3 a) where
  crossprod (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (y1*z2-y2*z1) (z1*x2-z2*x1) (x1*y2-x2*y1)

instance Num a => Determinant a (Vec3 a, Vec3 a, Vec3 a) where
  det (u,v,w) = u &. (v &^ w)

instance (Num a, Storable a) => Storable (Vec3 a) where
  sizeOf    _ = (sizeOf (undefined :: a) * 3 + 3) .&. 252
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    return (Vec3 x y z)

  poke q (Vec3 x y z) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z

instance Num a => Dimension (Vec3 a) where dim _ = 3

instance Functor Vec3 where
  fmap f (Vec3 a b c) = Vec3 (f a) (f b) (f c)

instance Foldable Vec3 where
    foldr f acc (Vec3 a b c) = f a (f b (f c acc))
    foldl f acc (Vec3 a b c) = acc `f` c `f` b `f` a
    foldr1 f    (Vec3 a b c) = a `f` b `f` c
    foldl1 f    (Vec3 a b c) = c `f` b `f` a
    null                     = const False
    length                   = const 3
    maximum     (Vec3 a b c) = a `max` b `max` c
    minimum     (Vec3 a b c) = a `min` b `min` c
    sum         (Vec3 a b c) = a + b + c
    product     (Vec3 a b c) = a * b * c

instance HasOne Vec3 where
  _1 (Vec3 x _ _) = x

instance HasTwo Vec3 where
  _2 (Vec3 _ y _) = y

instance HasThree Vec3 where
  _3 (Vec3 _ _ z) = z

--------------------------------------------------------------------------------
-- Vec4 instances
instance HasV2 Vec4 where
  getV2 (Vec4 x y _ _) = Vec2 x y

instance HasV3 Vec4 where
  getV3 (Vec4 x y z _) = Vec3 x y z

instance HasV4 Vec4 where
  getV4 = id

instance Num a => AbelianGroup (Vec4 a) where
  (&+) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1+x2) (y1+y2) (z1+z2) (w1+w2)
  (&-) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1-x2) (y1-y2) (z1-z2) (w1-w2)
  neg  (Vec4 x y z w)                      = Vec4 (-x) (-y) (-z) (-w)
  zero = Vec4 0 0 0 0

instance Num a => LinearMap a Vec4 where
  scalarMul s (Vec4 x y z w) = Vec4 (s*x) (s*y) (s*z) (s*w)
  mapVec    f (Vec4 x y z w) = Vec4 (f x) (f y) (f z) (f w)

instance Num a => DotProd a Vec4 where
  (&.) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = x1*x2 + y1*y2 + z1*z2 + w1*w2

instance Floating a => Norm a Vec4

instance Num a => Pointwise (Vec4 a) where
  pointwise (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1*x2) (y1*y2) (z1*z2) (w1*w2)

{-
instance Show Vec4 where
  show (Vec4 x y z w) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " , " ++ show w ++ " )"
-}

instance (Num a, Random a) => Random (Vec4 a) where
  random = randomR (Vec4 (-1) (-1) (-1) (-1),Vec4 1 1 1 1)
  randomR (Vec4 a b c d, Vec4 e f g h) gen =
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2
        (w,gen4) = randomR (d,h) gen3
    in (Vec4 x y z w, gen4)

instance (Num a, Storable a) => Storable (Vec4 a) where
  sizeOf    _ = 4 * sizeOf (undefined :: a)
  alignment _ = sizeOf (undefined :: a)

  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    w <- peekByteOff p (3*k)
    return (Vec4 x y z w)

  poke q (Vec4 x y z w) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z
    pokeByteOff p (3*k) w

instance Num a => Dimension (Vec4 a) where dim _ = 4

instance Functor Vec4 where
  fmap f (Vec4 a b c d) = Vec4 (f a) (f b) (f c) (f d)

instance Foldable Vec4 where
    foldr f acc (Vec4 a b c d) = f a (f b (f c (f d acc)))
    foldl f acc (Vec4 a b c d) = acc `f` d `f` c `f` b `f` a
    foldr1 f    (Vec4 a b c d) = a `f` b `f` c `f` d
    foldl1 f    (Vec4 a b c d) = d `f` c `f` b `f` a
    null                       = const False
    length                     = const 4
    maximum     (Vec4 a b c d) = a `max` b `max` c `max` d
    minimum     (Vec4 a b c d) = a `min` b `min` c `min` d
    sum         (Vec4 a b c d) = a + b + c + d
    product     (Vec4 a b c d) = a * b * c * d

instance HasOne Vec4 where
  _1 (Vec4 x _ _ _) = x

instance HasTwo Vec4 where
  _2 (Vec4 _ y _ _) = y

instance HasThree Vec4 where
  _3 (Vec4 _ _ z _) = z

instance HasFour Vec4 where
  _4 (Vec4 _ _ _ w) = w

--------------------------------------------------------------------------------
-- Extend instances

instance Num a => Extend a Vec2 Vec3 where
  extendHeadZero   (Vec2 x y) = Vec3 0 x y
  extendHeadWith t (Vec2 x y) = Vec3 t x y
  trimHead       (Vec3 _ x y) = Vec2 x y
  extendZero       (Vec2 x y) = Vec3 x y 0
  extendWith     t (Vec2 x y) = Vec3 x y t
  trim           (Vec3 x y _) = Vec2 x y

instance Num a => Extend a Vec2 Vec4 where
  extendHeadZero   (Vec2 x y) = Vec4 0 0 x y
  extendHeadWith t (Vec2 x y) = Vec4 t t x y
  trimHead     (Vec4 _ _ x y) = Vec2 x y
  extendZero       (Vec2 x y) = Vec4 x y 0 0
  extendWith     t (Vec2 x y) = Vec4 x y t t
  trim         (Vec4 x y _ _) = Vec2 x y

instance Num a => Extend a Vec3 Vec4 where
  extendHeadZero   (Vec3 x y z) = Vec4 0 x y z
  extendHeadWith t (Vec3 x y z) = Vec4 t x y z
  trimHead       (Vec4 _ x y z) = Vec3 x y z
  extendZero       (Vec3 x y z) = Vec4 x y z 0
  extendWith     t (Vec3 x y z) = Vec4 x y z t
  trim           (Vec4 x y z _) = Vec3 x y z

-- --------------------------------- NFData ---------------------------------------------

instance NFData a => NFData (Vec2 a) where
  rnf (Vec2 a b) = a `seq` b `seq` ()

instance NFData a => NFData (Vec3 a) where
  rnf (Vec3 a b c) = a `seq` b `seq` c `seq` ()

instance NFData a => NFData (Vec4 a) where
  rnf (Vec4 a b c d) = a `seq` b `seq` c `seq` d `seq` ()

-- ---------------------------------- Unbox ----------------------------------------------

derivingUnbox "Vec2"
    [t| forall a . (Unbox a) => Vec2 a -> (a, a) |]
    [| \ (Vec2 x y) -> (x, y) |]
    [| \ (x, y) -> (Vec2 x y) |]

derivingUnbox "Vec3"
    [t| forall a . (Unbox a) => Vec3 a -> (a, a, a) |]
    [| \ (Vec3 x y z) -> (x, y, z) |]
    [| \ (x, y, z) -> (Vec3 x y z) |]

derivingUnbox "Vec4"
    [t| forall a . (Unbox a) => Vec4 a -> (a, a, a, a) |]
    [| \ (Vec4 x y z t) -> (x, y, z, t) |]
    [| \ (x, y, z, t) -> (Vec4 x y z t) |]

derivingUnbox "Normal2"
    [t| forall a . (Unbox a) => Normal2 a -> (a, a) |]
    [| \ (Normal2 (Vec2 x y)) -> (x, y) |]
    [| \ (x, y) -> (Normal2 (Vec2 x y)) |]

derivingUnbox "Normal3"
    [t| forall a . (Unbox a) => Normal3 a -> (a, a, a) |]
    [| \ (Normal3 (Vec3 x y z)) -> (x, y, z) |]
    [| \ (x, y, z) -> (Normal3 (Vec3 x y z)) |]

derivingUnbox "Normal4"
    [t| forall a . (Unbox a) => Normal4 a -> (a, a, a, a) |]
    [| \ (Normal4 (Vec4 x y z t)) -> (x, y, z, t) |]
    [| \ (x, y, z, t) -> (Normal4 (Vec4 x y z t)) |]
