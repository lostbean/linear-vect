{-# LANGUAGE ScopedTypeVariables #-}

module SerializationSpec (spec) where

import Data.Vector.Unboxed as U
import Foreign.Marshal.Alloc
import Foreign.Storable
import Linear.Arbitrary ()
import Linear.Mat
import Linear.Vect
import Test.Hspec
import Test.QuickCheck
import Test.QuickCheck.Monadic

spec :: Spec
spec = do
    describe "Storable instances" $ do
        it "Vec3 storable" $
            property $
                \(v :: Vec3 Double) -> monadicIO $ do
                    v' <- run $ alloca $ \ptr -> do
                        poke ptr v
                        peek ptr
                    assert $ v == v'
        it "Mat4 storable" $
            property $
                \(m :: Mat4 Double) -> monadicIO $ do
                    m' <- run $ alloca $ \ptr -> do
                        poke ptr m
                        peek ptr
                    assert $ m == m'

    describe "Unbox instances" $ do
        it "Vec2 unbox" $
            property $
                \(l :: [Vec2 Double]) ->
                    let v = U.fromList l
                        l' = U.toList v
                     in l == l'
        it "Vec4 unbox" $
            property $
                \(l :: [Vec4 Double]) ->
                    let v = U.fromList l
                        l' = U.toList v
                     in l == l'
