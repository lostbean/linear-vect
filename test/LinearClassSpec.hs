{-# LANGUAGE ScopedTypeVariables #-}

module LinearClassSpec (spec) where

import Linear.Arbitrary ()
import Linear.Mat
import Linear.Vect
import Test.Hspec
import Test.QuickCheck

spec :: Spec
spec = do
    describe "Linear.Class" $ do
        it "normalize produces a unit vector" $
            property $
                \(v :: Vec3 Double) ->
                    let v' = normalize v
                     in if norm v > 1e-10 then abs (norm v' - 1.0) < 1e-10 else True
        it "distance is symmetric" $
            property $
                \(v1 :: Vec3 Double) (v2 :: Vec3 Double) -> distance v1 v2 =~ distance v2 v1
        it "angle is symmetric" $
            property $
                \(v1 :: Vec3 Double) (v2 :: Vec3 Double) ->
                    if norm v1 > 1e-10 && norm v2 > 1e-10
                        then abs (angle v1 v2 - angle v2 v1) < 1e-10
                        else True
        it "project' is orthogonal to the normal" $
            property $
                \(v :: Vec3 Double) (n :: Normal3 Double) ->
                    let p = project' v n
                     in isMainlyZero (p &. fromNormal n :: Double)
  where
    infix 4 =~
    (=~) :: Double -> Double -> Bool
    a =~ b = abs (a - b) < 1e-10
