{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

module ExhaustiveSpec (spec) where

import Control.DeepSeq
import qualified Data.Foldable as F
import Data.Vector.Unboxed as U
import Foreign.Marshal.Alloc
import Foreign.Storable
import Linear.Arbitrary ()
import Linear.Class
import Linear.Decomp
import Linear.Mat
import Linear.Vect
import System.Random
import Test.Hspec
import Test.QuickCheck
import Test.QuickCheck.Monadic
import Prelude hiding (length, zipWith)
import qualified Prelude as P

-- ================================== Helper ==================================

isNearZeroV :: (F.Foldable f, NearZero a, Ord a) => f a -> Bool
isNearZeroV = F.all isMainlyZero

nearEq :: forall f a. (AbelianGroup (f a), F.Foldable f, Functor f, Norm a f, LinearMap a f, NearZero a, Ord a, Fractional a) => f a -> f a -> Bool
nearEq a b =
    let d = distance a b
        limit = 1e-8 * (1 + P.max (norm a) (norm b))
     in d < limit

nearEqMat :: forall m. (AbelianGroup (m Double), F.Foldable m, Functor m, MatrixNorms Double (m Double)) => m Double -> m Double -> Bool
nearEqMat a b =
    let d = matrixDistance a b
        na = frobeniusNorm a
        nb = frobeniusNorm b
        limit = 1e-8 * (1 + P.max na nb)
     in d < limit

-- ================================== Spec Generators ==================================

testAbelianGroup :: forall f a. (AbelianGroup (f a), F.Foldable f, Functor f, Norm a f, LinearMap a f, Show (f a), Arbitrary (f a), NearZero a, Ord a, Fractional a) => Spec
testAbelianGroup = describe "AbelianGroup" $ do
    it "commutativity: a &+ b == b &+ a" $ property $ \(a :: f a) (b :: f a) -> nearEq (a &+ b) (b &+ a)
    it "associativity: (a &+ b) &+ c == a &+ (b &+ c)" $ property $ \(a :: f a) (b :: f a) (c :: f a) -> nearEq ((a &+ b) &+ c) (a &+ (b &+ c))
    it "identity: a &+ zero == a" $ property $ \(a :: f a) -> nearEq (a &+ zero) a
    it "inverse: a &+ neg a == zero" $ property $ \(a :: f a) -> nearEq (a &+ neg a) zero

testLinearMap :: forall a f. (LinearMap a f, F.Foldable f, Functor f, Norm a f, Show (f a), Show a, Arbitrary (f a), Arbitrary a, Fractional a, Ord a, NearZero a) => Spec
testLinearMap = describe "LinearMap" $ do
    it "distributivity 1: s *& (a &+ b) == s *& a &+ s *& b" $ property $ \(s :: a) (a :: f a) (b :: f a) -> nearEq (s *& (a &+ b)) (s *& a &+ s *& b)
    it "distributivity 2: (s1 + s2) *& a == s1 *& a &+ s2 *& a" $ property $ \(s1 :: a) (s2 :: a) (a :: f a) -> nearEq ((s1 + s2) *& a) (s1 *& a &+ s2 *& a)
    it "associativity: s1 *& (s2 *& a) == (s1 * s2) *& a" $ property $ \(s1 :: a) (s2 :: a) (a :: f a) -> nearEq (s1 *& (s2 *& a)) ((s1 * s2) *& a)

testDotProd :: forall a f. (DotProd a f, Show (f a), Arbitrary (f a), Arbitrary a, Fractional a, Ord a, NearZero a) => Spec
testDotProd = describe "DotProd" $ do
    it "commutativity: a &. b == b &. a" $ property $ \(a :: f a) (b :: f a) -> abs ((a &. b) - (b &. a)) < 1e-10 * (1 + abs (a &. b))
    it "positivity: a &. a >= -1e-10" $ property $ \(a :: f a) -> (a &. a) >= -1e-10
    it "normsqr: normsqr a == a &. a" $ property $ \(a :: f a) -> isMainlyZero (normsqr a - (a &. a))

testFunctor :: forall f. (Functor f, Eq (f Int), Show (f Int), Arbitrary (f Int)) => Spec
testFunctor = describe "Functor" $ do
    it "identity: fmap id == id" $ property $ \(a :: f Int) -> fmap id a == a
    it "composition: fmap (f . g) == fmap f . fmap g" $ property $ \(a :: f Int) -> fmap ((+ 1) . (* 2)) a == (fmap (+ 1) . fmap (* 2)) a

testFoldableFull :: forall f. (F.Foldable f, Dimension (f Double), Show (f Double), Arbitrary (f Double)) => Spec
testFoldableFull = describe "Foldable Full" $ do
    it "toList length" $ property $ \(a :: f Double) -> P.length (F.toList a) == dim a
    it "sum matches" $ property $ \(a :: f Double) -> abs (F.sum a - P.sum (F.toList a)) < 1e-10 * (1 + abs (F.sum a))
    it "product matches" $ property $ \(a :: f Double) -> abs (F.product a - P.product (F.toList a)) < 1e-10 * (1 + abs (F.product a))
    it "maximum/minimum" $ property $ \(a :: f Double) -> F.maximum a == P.maximum (F.toList a) && F.minimum a == P.minimum (F.toList a)
    it "foldr/foldl" $ property $ \(a :: f Double) -> F.foldr (+) 0 a == P.foldr (+) 0 (F.toList a) && F.foldl (+) 0 a == P.foldl (+) 0 (F.toList a)

testStorable :: forall g. (Storable g, Eq g, Show g, Arbitrary g) => Spec
testStorable = it "Storable" $ property $ \(v :: g) -> monadicIO $ do v' <- run $ alloca $ \ptr -> poke ptr v >> peek ptr; assert $ v == v'

testUnboxFull :: forall g. (U.Unbox g, Eq g, Show g, Arbitrary g) => Spec
testUnboxFull = describe "Unbox Full" $ do
    it "toList/fromList" $ property $ \(l :: [g]) -> U.toList (U.fromList l) == l
    it "replicate" $ property $ \(v :: g) -> U.toList (U.replicate 3 v) == [v, v, v]

testNFData :: forall g. (NFData g, Show g, Arbitrary g) => Spec
testNFData = it "NFData" $ property $ \(v :: g) -> rnf v `seq` True

testPrettyShow :: forall g. (PrettyShow g, Show g, Arbitrary g) => Spec
testPrettyShow = it "PrettyShow" $ property $ \(v :: g) -> P.length (showPretty v) > 0

-- ================================== Spec ==================================

spec :: Spec
spec = do
    describe "Vec2 Double" $ do
        testAbelianGroup @Vec2 @Double
        testLinearMap @Double @Vec2
        testDotProd @Double @Vec2
        testFunctor @Vec2
        testFoldableFull @Vec2
        testStorable @(Vec2 Double)
        testUnboxFull @(Vec2 Double)
        testNFData @(Vec2 Double)
        testPrettyShow @(Vec2 Double)
        it "accessors" $ property $ \(x :: Double) (y :: Double) -> let v = Vec2 x y in _x v == x && _y v == y
        it "mkVec/unVec" $ property $ \(x :: Double) (y :: Double) -> unVec2 (mkVec2 (x, y)) == (x, y)
        it "Pointwise" $ property $ \(a :: Vec2 Double) (b :: Vec2 Double) -> F.toList (a &! b) == P.zipWith (*) (F.toList a) (F.toList b)
        it "HasOne/Two" $ property $ \(a :: Vec2 Double) -> _1 a == _x a && _2 a == _y a
        it "HasV2" $ property $ \(a :: Vec2 Double) -> getV2 a == a
        it "vece1" $ property $ \(x :: Double) -> let v = vece1 x :: Vec2 Double in _1 v == x && _2 v == 0

    describe "Vec3 Double" $ do
        testAbelianGroup @Vec3 @Double
        testLinearMap @Double @Vec3
        testDotProd @Double @Vec3
        testFunctor @Vec3
        testFoldableFull @Vec3
        testStorable @(Vec3 Double)
        testUnboxFull @(Vec3 Double)
        testNFData @(Vec3 Double)
        testPrettyShow @(Vec3 Double)
        it "accessors" $ property $ \(x :: Double) (y :: Double) (z :: Double) -> let v = Vec3 x y z in _1 v == x && _2 v == y && _3 v == z
        it "mkVec/unVec" $ property $ \(x :: Double) (y :: Double) (z :: Double) -> unVec3 (mkVec3 (x, y, z)) == (x, y, z)
        it "CrossProd" $ property $ \(a :: Vec3 Double) (b :: Vec3 Double) ->
            let c = a &^ b in abs (c &. a) < 1e-10 * (1 + norm c * norm a) && abs (c &. b) < 1e-10 * (1 + norm c * norm b)
        it "HasOne/Two/Three" $ property $ \(a :: Vec3 Double) -> _1 a == _1 a && _2 a == _2 a && _3 a == _3 a
        it "HasV2" $ property $ \(a :: Vec3 Double) -> let (Vec3 x y _) = a in getV2 a == Vec2 x y
        it "HasV3" $ property $ \(a :: Vec3 Double) -> getV3 a == a
        it "vece1" $ property $ \(x :: Double) -> let v = vece1 x :: Vec3 Double in _1 v == x && _2 v == 0 && _3 v == 0

    describe "Vec4 Double" $ do
        testAbelianGroup @Vec4 @Double
        testLinearMap @Double @Vec4
        testDotProd @Double @Vec4
        testFunctor @Vec4
        testFoldableFull @Vec4
        testStorable @(Vec4 Double)
        testUnboxFull @(Vec4 Double)
        testNFData @(Vec4 Double)
        testPrettyShow @(Vec4 Double)
        it "accessors" $ property $ \(x :: Double) (y :: Double) (z :: Double) (w :: Double) -> let v = Vec4 x y z w in _1 v == x && _2 v == y && _3 v == z && _4 v == w
        it "mkVec/unVec" $ property $ \(x :: Double) (y :: Double) (z :: Double) (w :: Double) -> unVec4 (mkVec4 (x, y, z, w)) == (x, y, z, w)
        it "HasOne/Two/Three/Four" $ property $ \(a :: Vec4 Double) -> _1 a == _1 a && _2 a == _2 a && _3 a == _3 a && _4 a == _4 a
        it "HasV2" $ property $ \(a :: Vec4 Double) -> let (Vec4 x y _ _) = a in getV2 a == Vec2 x y
        it "HasV3" $ property $ \(a :: Vec4 Double) -> let (Vec4 x y z _) = a in getV3 a == Vec3 x y z
        it "HasV4" $ property $ \(a :: Vec4 Double) -> getV4 a == a
        it "vece1" $ property $ \(x :: Double) -> let v = vece1 x :: Vec4 Double in _1 v == x && _2 v == 0 && _3 v == 0 && _4 v == 0

    describe "Normal Double" $ do
        it "Normal2 unit length" $ property $ \(v :: Normal2 Double) -> isMainlyZero (norm (fromNormal v) - 1)
        it "Normal3 unit length" $ property $ \(v :: Normal3 Double) -> isMainlyZero (norm (fromNormal v) - 1)
        it "Normal4 unit length" $ property $ \(v :: Normal4 Double) -> isMainlyZero (norm (fromNormal v) - 1)
        it "Normal3 CrossProd" $ property $ \(a :: Normal3 Double) (b :: Normal3 Double) ->
            let c = a &^ b in isMainlyZero (norm (fromNormal c) - 1) || isNearZeroV (fromNormal a &^ fromNormal b)
        testPrettyShow @(Normal2 Double)
        testPrettyShow @(Normal3 Double)
        testPrettyShow @(Normal4 Double)
        testUnboxFull @(Normal2 Double)
        testUnboxFull @(Normal3 Double)
        testUnboxFull @(Normal4 Double)
        testStorable @(Normal2 Double)
        testStorable @(Normal3 Double)
        testStorable @(Normal4 Double)
        it "fromNormalRadius" $ property $ \(t :: Double) (n :: Normal3 Double) -> nearEq (fromNormalRadius t n) (t *& fromNormal n)
        it "flipNormal" $ property $ \(n :: Normal3 Double) -> nearEq (fromNormal (flipNormal n)) (neg (fromNormal n))
        it "HasV2 Normal2" $ property $ \(n :: Normal2 Double) -> getV2 n == fromNormal n
        it "HasV2 Normal3" $ property $ \(n :: Normal3 Double) -> getV2 n == getV2 (fromNormal n)
        it "HasV2 Normal4" $ property $ \(n :: Normal4 Double) -> getV2 n == getV2 (fromNormal n)
        it "HasV3 Normal3" $ property $ \(n :: Normal3 Double) -> getV3 n == fromNormal n
        it "HasV3 Normal4" $ property $ \(n :: Normal4 Double) -> getV3 n == getV3 (fromNormal n)
        it "HasV4 Normal4" $ property $ \(n :: Normal4 Double) -> getV4 n == fromNormal n
        it "Random Normal2" $ property $ \(s :: Int) -> let (n, _) = random (mkStdGen s) :: (Normal2 Double, StdGen) in isMainlyZero (norm (fromNormal n) - 1)
        it "Random Normal3" $ property $ \(s :: Int) -> let (n, _) = random (mkStdGen s) :: (Normal3 Double, StdGen) in isMainlyZero (norm (fromNormal n) - 1)
        it "Random Normal4" $ property $ \(s :: Int) -> let (n, _) = random (mkStdGen s) :: (Normal4 Double, StdGen) in isMainlyZero (norm (fromNormal n) - 1)

    describe "Mat2 Double" $ do
        it "AbelianGroup" $ property $ \(a :: Mat2 Double) (b :: Mat2 Double) -> nearEqMat a b || nearEqMat (a &+ b) (b &+ a)
        it "LinearMap" $ property $ \(s :: Double) (a :: Mat2 Double) (b :: Mat2 Double) -> nearEqMat (s *& (a &+ b)) (s *& a &+ s *& b)
        it "SquareMatrix identity" $ property $ \(v :: Vec2 Double) -> nearEq (idmtx *. v) v
        it "SquareMatrix inverse" $ property $ \(m :: Mat2 Double) ->
            abs (det m :: Double) > 1e-2 ==> nearEqMat (m .*. inverse m) idmtx
        testStorable @(Mat2 Double)
        testUnboxFull @(Mat2 Double)
        testNFData @(Mat2 Double)
        testPrettyShow @(Mat2 Double)
        it "Determinant idmtx" $ (det (idmtx :: Mat2 Double) :: Double) `shouldBe` 1
        it "Transpose" $ property $ \(m :: Mat2 Double) -> transpose (transpose m) == m
        it "Pointwise" $ property $ \(a :: Mat2 Double) (b :: Mat2 Double) ->
            F.toList (a &! b) == P.zipWith (*) (F.toList a) (F.toList b)
        it "Diagonal" $ property $ \(v :: Vec2 Double) -> diagVec (diag v :: Mat2 Double) == v
        it "Tensor outer" $ property $ \(u :: Vec2 Double) (v :: Vec2 Double) ->
            let m = outer u v in nearEq (m *. v) (u &* normsqr v)
        it "HasRow" $ property $ \(m :: Mat2 Double) -> let (Mat2 r1 r2) = m in _R1 m == r1 && _R2 m == r2

    describe "Mat3 Double" $ do
        it "AbelianGroup" $ property $ \(a :: Mat3 Double) (b :: Mat3 Double) -> nearEqMat a b || nearEqMat (a &+ b) (b &+ a)
        it "LinearMap" $ property $ \(s :: Double) (a :: Mat3 Double) (b :: Mat3 Double) -> nearEqMat (s *& (a &+ b)) (s *& a &+ s *& b)
        it "SquareMatrix identity" $ property $ \(v :: Vec3 Double) -> nearEq (idmtx *. v) v
        it "SquareMatrix inverse" $ property $ \(m :: Mat3 Double) ->
            abs (det m :: Double) > 1e-2 ==> nearEqMat (m .*. inverse m) idmtx
        testStorable @(Mat3 Double)
        testUnboxFull @(Mat3 Double)
        testNFData @(Mat3 Double)
        testPrettyShow @(Mat3 Double)
        it "Determinant idmtx" $ (det (idmtx :: Mat3 Double) :: Double) `shouldBe` 1
        it "Transpose" $ property $ \(m :: Mat3 Double) -> transpose (transpose m) == m
        it "HasRow" $ property $ \(m :: Mat3 Double) -> let (Mat3 r1 r2 r3) = m in _R1 m == r1 && _R2 m == r2 && _R3 m == r3

    describe "Mat4 Double" $ do
        it "AbelianGroup" $ property $ \(a :: Mat4 Double) (b :: Mat4 Double) -> nearEqMat a b || nearEqMat (a &+ b) (b &+ a)
        it "LinearMap" $ property $ \(s :: Double) (a :: Mat4 Double) (b :: Mat4 Double) -> nearEqMat (s *& (a &+ b)) (s *& a &+ s *& b)
        it "SquareMatrix identity" $ property $ \(v :: Vec4 Double) -> nearEq (idmtx *. v) v
        it "SquareMatrix inverse" $ property $ \(m :: Mat4 Double) ->
            abs (det m :: Double) > 1e-2 ==> nearEqMat (m .*. inverse m) idmtx
        testStorable @(Mat4 Double)
        testUnboxFull @(Mat4 Double)
        testNFData @(Mat4 Double)
        testPrettyShow @(Mat4 Double)
        it "Determinant idmtx" $ (det (idmtx :: Mat4 Double) :: Double) `shouldBe` 1
        it "Transpose" $ property $ \(m :: Mat4 Double) -> transpose (transpose m) == m
        it "HasRow" $ property $ \(m :: Mat4 Double) -> let (Mat4 r1 r2 r3 r4) = m in _R1 m == r1 && _R2 m == r2 && _R3 m == r3 && _R4 m == r4

    describe "Rectangular Mat" $ do
        it "Mat2x3 transpose" $ property $ \(m :: Mat2x3 Double) -> transpose (transpose m) == m
        it "Mat2x4 transpose" $ property $ \(m :: Mat2x4 Double) -> transpose (transpose m) == m
        it "Mat3x2 transpose" $ property $ \(m :: Mat3x2 Double) -> transpose (transpose m) == m
        it "Mat3x4 transpose" $ property $ \(m :: Mat3x4 Double) -> transpose (transpose m) == m
        it "Mat4x2 transpose" $ property $ \(m :: Mat4x2 Double) -> transpose (transpose m) == m
        it "Mat4x3 transpose" $ property $ \(m :: Mat4x3 Double) -> transpose (transpose m) == m
        testStorable @(Mat2x3 Double)
        testStorable @(Mat2x4 Double)
        testStorable @(Mat3x2 Double)
        testStorable @(Mat3x4 Double)
        testStorable @(Mat4x2 Double)
        testStorable @(Mat4x3 Double)

    describe "Ortho Double" $ do
        it "Ortho2 properties" $ property $ \(o :: Ortho2 Double) ->
            let m = fromOrtho o in nearEqMat (m .*. transpose m) idmtx
        it "Ortho3 properties" $ property $ \(o :: Ortho3 Double) ->
            let m = fromOrtho o in nearEqMat (m .*. transpose m) idmtx
        it "Ortho4 properties" $ property $ \(o :: Ortho4 Double) ->
            let m = fromOrtho o in nearEqMat (m .*. transpose m) idmtx
        testPrettyShow @(Ortho2 Double)
        testPrettyShow @(Ortho3 Double)
        testPrettyShow @(Ortho4 Double)

    describe "Projective Double" $ do
        it "Proj3 properties" $ property $ \(p :: Proj3 Double) ->
            let m = fromProjective p in abs (det m :: Double) > 1e-2 ==> nearEqMat (m .*. inverse m) idmtx
        it "Proj4 properties" $ property $ \(p :: Proj4 Double) ->
            let m = fromProjective p in abs (det m :: Double) > 1e-2 ==> nearEqMat (m .*. inverse m) idmtx
        testPrettyShow @(Proj3 Double)
        testPrettyShow @(Proj4 Double)

    describe "Decomp" $ do
        it "qrHouse 3x3" $ property $ \(m :: Mat3 Double) ->
            let (q, r) = qrHouse m
             in nearEqMat (q .*. r) m
        it "qrGram 3x3" $ property $ \(m :: Mat3 Double) ->
            abs (det m :: Double) > 1e-2 ==>
                let (q, r) = qrGram m
                 in nearEqMat (q .*. r) m
        it "symmEigen 3x3" $ property $ \(m :: Mat3 Double) ->
            let mt = m .*. transpose m -- symmetric
                (q, v) = symmEigen mt
                m' = q .*. diag v .*. transpose q
             in (matrixDistance mt m' :: Double) < 1e-6 * (1 + frobeniusNorm mt)
        it "hessen 3x3" $ property $ \(m :: Mat3 Double) ->
            let h = hessen m
             in dim h == 3
        it "orthoColsHouse 3x3" $ property $ \(m :: Mat3 Double) ->
            abs (det m :: Double) > 1e-2 ==>
                let q = orthoColsHouse m
                 in nearEqMat (q .*. transpose q) idmtx
        it "orthoColsGram 2x2" $ property $ \(m :: Mat2 Double) ->
            abs (det m :: Double) > 1e-2 ==>
                let q = orthoColsGram m
                 in nearEqMat (q .*. transpose q) idmtx
        it "orthoRowsHouse 4x4" $ property $ \(m :: Mat4 Double) ->
            abs (det m :: Double) > 1e-2 ==>
                let q = orthoRowsHouse m
                 in nearEqMat (q .*. transpose q) idmtx
        it "householder property" $ property $ \(n :: Normal3 Double) ->
            let m = householder n :: Mat3 Double
                v = fromNormal n
             in nearEq (m *. v) (neg v)

    describe "Extend" $ do
        it "Vec2 -> Vec3" $ property $ \(v :: Vec2 Double) -> trim (extendZero v :: Vec3 Double) == v
        it "Vec2 -> Vec4" $ property $ \(v :: Vec2 Double) -> trim (extendZero v :: Vec4 Double) == v
        it "Vec3 -> Vec4" $ property $ \(v :: Vec3 Double) -> trim (extendZero v :: Vec4 Double) == v
        it "Mat2 -> Mat3" $ property $ \(m :: Mat2 Double) -> trim (extendZero m :: Mat3 Double) == m
        it "Mat2 -> Mat4" $ property $ \(m :: Mat2 Double) -> trim (extendZero m :: Mat4 Double) == m
        it "Mat3 -> Mat4" $ property $ \(m :: Mat3 Double) -> trim (extendZero m :: Mat4 Double) == m
        it "Vec2 -> Vec3 Head" $ property $ \(v :: Vec2 Double) -> trimHead (extendHeadZero v :: Vec3 Double) == v
        it "Vec2 -> Vec4 Head" $ property $ \(v :: Vec2 Double) -> trimHead (extendHeadZero v :: Vec4 Double) == v
        it "Vec3 -> Vec4 Head" $ property $ \(v :: Vec3 Double) -> trimHead (extendHeadZero v :: Vec4 Double) == v
        it "extendWith Vec2 -> Vec3" $ property $ \(t :: Double) (v :: Vec2 Double) ->
            let v3 = extendWith t v :: Vec3 Double in _3 v3 == t && trim v3 == v
        it "extendHeadWith Vec2 -> Vec3" $ property $ \(t :: Double) (v :: Vec2 Double) ->
            let v3 = extendHeadWith t v :: Vec3 Double in _1 v3 == t && trimHead v3 == v

    describe "Linear.Class extra" $ do
        it "angle" $ property $ \(v1 :: Vec3 Double) (v2 :: Vec3 Double) ->
            norm v1 > 1e-10 && norm v2 > 1e-10 ==>
                let a = angle v1 v2 in a >= 0 && a <= pi + 1e-10
        it "angle'" $ property $ \(n1 :: Normal3 Double) (n2 :: Normal3 Double) ->
            let a = angle' n1 n2 in a >= 0 && a <= pi + 1e-10
        it "householderOrtho" $ property $ \(n :: Normal3 Double) ->
            let o = householderOrtho n :: Ortho3 Double
                m = fromOrtho o
             in nearEqMat (m .*. transpose m) idmtx
        it "project" $ property $ \(v1 :: Vec3 Double) (v2 :: Vec3 Double) ->
            normsqr v2 > 1e-10 ==>
                let p = project v1 v2
                 in abs (p &. v2 :: Double) < 1e-8 * (1 + norm p * norm v2)
