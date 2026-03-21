{-# LANGUAGE ScopedTypeVariables #-}

module LinearMatSpec (spec) where

import Linear.Arbitrary ()
import Linear.Mat
import Test.Hspec
import Test.QuickCheck

spec :: Spec
spec = do
    describe "Mat2" $ do
        it "identity multiplication works" $
            property $
                \(m :: Mat2 Double) -> m .*. idmtx == m && idmtx .*. m == m
        it "matrix addition is commutative" $
            property $
                \(m1 :: Mat2 Double) (m2 :: Mat2 Double) -> m1 &+ m2 == m2 &+ m1
        it "determinant of identity is 1" $
            det (idmtx :: Mat2 Double) `shouldBe` (1 :: Double)

    describe "Mat3" $ do
        it "identity multiplication works" $
            property $
                \(m :: Mat3 Double) -> m .*. idmtx == m && idmtx .*. m == m
        it "determinant of identity is 1" $
            det (idmtx :: Mat3 Double) `shouldBe` (1 :: Double)

    describe "Mat4" $ do
        it "identity multiplication works" $
            property $
                \(m :: Mat4 Double) -> m .*. idmtx == m && idmtx .*. m == m
        it "determinant of identity is 1" $
            det (idmtx :: Mat4 Double) `shouldBe` (1 :: Double)
