{-# OPTIONS_GHC -Wno-orphans #-}

module Linear.Arbitrary () where

import Linear.Class
import Linear.Mat
import Linear.Vect
import System.Random
import Test.QuickCheck

-- Vectors
instance (Arbitrary a) => Arbitrary (Vec2 a) where arbitrary = Vec2 <$> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Vec3 a) where arbitrary = Vec3 <$> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Vec4 a) where arbitrary = Vec4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

-- Normal vectors
instance (Arbitrary a, Floating a, Ord a) => Arbitrary (Normal2 a) where
    arbitrary = toNormalUnsafe . normalize <$> (arbitrary `suchThat` (\v -> norm v > 1e-10))
instance (Arbitrary a, Floating a, Ord a) => Arbitrary (Normal3 a) where
    arbitrary = toNormalUnsafe . normalize <$> (arbitrary `suchThat` (\v -> norm v > 1e-10))
instance (Arbitrary a, Floating a, Ord a) => Arbitrary (Normal4 a) where
    arbitrary = toNormalUnsafe . normalize <$> (arbitrary `suchThat` (\v -> norm v > 1e-10))

-- Square matrices
instance (Arbitrary a) => Arbitrary (Mat2 a) where arbitrary = Mat2 <$> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3 a) where arbitrary = Mat3 <$> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4 a) where arbitrary = Mat4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

-- Rectangular matrices
instance (Arbitrary a) => Arbitrary (Mat2x3 a) where arbitrary = Mat2x3 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat2x4 a) where arbitrary = Mat2x4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3x2 a) where arbitrary = Mat3x2 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat3x4 a) where arbitrary = Mat3x4 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4x2 a) where arbitrary = Mat4x2 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
instance (Arbitrary a) => Arbitrary (Mat4x3 a) where arbitrary = Mat4x3 <$> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary

-- Orthogonal matrices
instance (Arbitrary a, Floating a, Ord a, Random a) => Arbitrary (Ortho2 a) where
    arbitrary = fst . random . mkStdGen <$> arbitrary
instance (Arbitrary a, Floating a, Ord a, Random a) => Arbitrary (Ortho3 a) where
    arbitrary = fst . random . mkStdGen <$> arbitrary
instance (Arbitrary a, Floating a, Ord a, Random a) => Arbitrary (Ortho4 a) where
    arbitrary = fst . random . mkStdGen <$> arbitrary

-- Projective matrices
instance (Arbitrary a, Fractional a) => Arbitrary (Proj3 a) where arbitrary = toProjectiveUnsafe <$> arbitrary
instance (Arbitrary a, Fractional a) => Arbitrary (Proj4 a) where arbitrary = toProjectiveUnsafe <$> arbitrary
