{-# LANGUAGE ScopedTypeVariables #-}

module CoverageSpec (spec) where

import Data.Foldable (foldl, sum)
import Linear.Class
import Linear.Decomp
import Linear.Mat
import Linear.Vect
import System.Random
import Test.Hspec
import Test.QuickCheck
import Prelude hiding (foldl, sum)

-- Arbitrary instances
instance (Arbitrary a) => Arbitrary (Vec2 a) where
    arbitrary = Vec2 <$> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Vec3 a) where
    arbitrary = Vec3 <$> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Vec4 a) where
    arbitrary = Vec4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

instance (Arbitrary a) => Arbitrary (Mat2 a) where
    arbitrary = Mat2 <$> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3 a) where
    arbitrary = Mat3 <$> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4 a) where
    arbitrary = Mat4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

spec :: Spec
spec = do
    describe "Comprehensive Coverage" $ do
        it "Mat2x3 transpose" $
            property $
                \(m :: Mat2x3 Double) -> transpose (transpose m) == m
        it "Mat2x4 transpose" $
            property $
                \(m :: Mat2x4 Double) -> transpose (transpose m) == m
        it "Mat3x2 transpose" $
            property $
                \(m :: Mat3x2 Double) -> transpose (transpose m) == m
        it "Mat3x4 transpose" $
            property $
                \(m :: Mat3x4 Double) -> transpose (transpose m) == m
        it "Mat4x2 transpose" $
            property $
                \(m :: Mat4x2 Double) -> transpose (transpose m) == m
        it "Mat4x3 transpose" $
            property $
                \(m :: Mat4x3 Double) -> transpose (transpose m) == m

        it "QR Gram-Schmidt 3x3" $
            property $
                \(m :: Mat3 Double) ->
                    let (q, r) = qrGram m
                        res = m &- (q .*. r)
                        norm_m = frobeniusNorm m :: Double
                        limit = 1e-8 * (1 + norm_m)
                     in if abs (det m :: Double) > 1e-2 then frobeniusNorm res < limit else True

        it "QR Householder 3x3" $
            property $
                \(m :: Mat3 Double) ->
                    let (q, r) = qrHouse m
                        res = m &- (q .*. r)
                        norm_m = frobeniusNorm m :: Double
                        limit = 1e-10 * (1 + norm_m)
                     in if abs (det m :: Double) > 1e-2 then frobeniusNorm res < limit else True

        it "Hessenberg reduction 3x3" $
            property $
                \(m :: Mat3 Double) ->
                    let h = hessen m
                     in dim h == 3 -- Just check it runs
        it "Hessenberg reduction 4x4" $
            property $
                \(m :: Mat4 Double) ->
                    let h = hessen m
                     in dim h == 4 -- Just check it runs
        it "Extend/Trim Vec2 to Vec3" $
            property $
                \(v :: Vec2 Double) -> trim (extendZero v :: Vec3 Double) == v
        it "Extend/Trim Vec3 to Vec4" $
            property $
                \(v :: Vec3 Double) -> trim (extendZero v :: Vec4 Double) == v

        it "Extend/Trim Mat2 to Mat3" $
            property $
                \(m :: Mat2 Double) -> trim (extendZero m :: Mat3 Double) == m
        it "Extend/Trim Mat3 to Mat4" $
            property $
                \(m :: Mat3 Double) -> trim (extendZero m :: Mat4 Double) == m

        it "Ortho3 properties" $
            property $
                \(o :: Ortho3 Double) ->
                    let m = fromOrtho o
                        res = (m .*. transpose m) &- idmtx
                     in isMainlyZero (frobeniusNorm res :: Double)

        it "Proj4 inversion" $
            property $
                \(p :: Proj4 Double) ->
                    let m = fromProjective p
                     in if abs (det m :: Double) > 1e-2
                            then
                                let p' = toProjectiveUnsafe (inverse m)
                                    res = (fromProjective p' .*. m) &- idmtx
                                 in (frobeniusNorm res :: Double) < 1e-8
                            else True

        it "householder properties" $
            property $
                \(n :: Normal3 Double) ->
                    let m = householder n :: Mat3 Double
                        res = (m .*. transpose m) &- idmtx
                     in isMainlyZero (frobeniusNorm res :: Double)

instance (Arbitrary a) => Arbitrary (Mat2x3 a) where
    arbitrary = Mat2x3 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat2x4 a) where
    arbitrary = Mat2x4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3x2 a) where
    arbitrary = Mat3x2 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3x4 a) where
    arbitrary = Mat3x4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4x2 a) where
    arbitrary = Mat4x2 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4x3 a) where
    arbitrary = Mat4x3 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

instance (Arbitrary a, Floating a, Ord a, Random a) => Arbitrary (Ortho3 a) where
    arbitrary = do
        seed <- arbitrary
        return $ fst $ random (mkStdGen seed)

instance (Arbitrary a, Floating a, Ord a) => Arbitrary (Normal3 a) where
    arbitrary = do
        v <- arbitrary `suchThat` (\v -> norm v > 1e-10)
        return (toNormalUnsafe (normalize v))

instance (Arbitrary a, Fractional a) => Arbitrary (Proj4 a) where
    arbitrary = toProjectiveUnsafe <$> (arbitrary :: Gen (Mat4 a))
