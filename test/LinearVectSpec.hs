{-# LANGUAGE ScopedTypeVariables #-}

module LinearVectSpec (spec) where

import Data.Foldable (all)
import Linear.Arbitrary ()
import Linear.Class
import Linear.Vect
import Test.Hspec
import Test.QuickCheck
import Prelude hiding (all)

-- Helper for near zero vectors
isNearZero :: (Foldable v, NearZero a, Ord a) => v a -> Bool
isNearZero = all isMainlyZero

spec :: Spec
spec = do
    describe "Vec2" $ do
        it "addition is commutative" $
            property $
                \(v1 :: Vec2 Double) (v2 :: Vec2 Double) -> v1 &+ v2 == v2 &+ v1
        it "addition is associative" $
            property $
                \(v1 :: Vec2 Double) (v2 :: Vec2 Double) (v3 :: Vec2 Double) ->
                    let left = (v1 &+ v2) &+ v3
                        right = v1 &+ (v2 &+ v3)
                     in isNearZero (left &- right)
        it "subtraction works" $
            property $
                \(v1 :: Vec2 Double) (v2 :: Vec2 Double) -> v1 &- v2 == v1 &+ (neg v2)
        it "dot product is commutative" $
            property $
                \(v1 :: Vec2 Double) (v2 :: Vec2 Double) -> v1 &. v2 == v2 &. v1

    describe "Vec3" $ do
        it "addition is commutative" $
            property $
                \(v1 :: Vec3 Double) (v2 :: Vec3 Double) -> v1 &+ v2 == v2 &+ v1
        it "dot product is commutative" $
            property $
                \(v1 :: Vec3 Double) (v2 :: Vec3 Double) -> v1 &. v2 == v2 &. v1
        it "cross product property: a &^ b is orthogonal to a" $
            property $
                \(v1 :: Vec3 Double) (v2 :: Vec3 Double) ->
                    let v3 = v1 &^ v2
                        dot = v3 &. v1
                        -- Use relative tolerance if values are large
                        norm_v1 = norm v1
                        norm_v3 = norm v3
                        limit = 1e-10 * (1 + norm_v1 * norm_v3)
                     in abs dot < limit

    describe "Vec4" $ do
        it "addition is commutative" $
            property $
                \(v1 :: Vec4 Double) (v2 :: Vec4 Double) -> v1 &+ v2 == v2 &+ v1
        it "dot product is commutative" $
            property $
                \(v1 :: Vec4 Double) (v2 :: Vec4 Double) -> v1 &. v2 == v2 &. v1
