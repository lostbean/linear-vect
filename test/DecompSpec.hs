module DecompSpec (spec) where

import Linear.Class
import Linear.Decomp
import Linear.Mat
import Linear.Vect
import Test.Hspec

-- ================================== Test Data ==================================

testData1 :: Mat3 Double -- Source <http://en.wikipedia.org/wiki/QR_decomposition>
testData1 = Mat3 (Vec3 12 (-51) 4) (Vec3 6 167 (-68)) (Vec3 (-4) 24 (-41))

testData2 :: Mat3 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData2 = Mat3 (Vec3 2 1 0) (Vec3 1 3 (-1)) (Vec3 0 (-1) 6)

testData3 :: Mat4 Double
testData3 = Mat4 (Vec4 0 10 3 9) (Vec4 10 12 6 15) (Vec4 3 6 0 7) (Vec4 9 15 7 8)

testData5 :: Mat2 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData5 = Mat2 (Vec2 2 1) (Vec2 1 3)

testQR :: (Transpose g g, MultSemiGroup g, AbelianGroup g, OrthoMatrix g, SquareMatrix g) => g -> g
testQR m = m &- (q .*. r)
  where
    (q, r) = qrHouse m

testEigen3 :: Mat3 Double -> [Double]
testEigen3 m = map (normsqr . foo) [(_1, _R1), (_2, _R2), (_3, _R3)]
  where
    foo (f1, f2) = (m &- (f1 value) *& idmtx) *. (f2 $ transpose vec)
    (vec, value) = symmEigen m

testEigen4 :: Mat4 Double -> [Double]
testEigen4 m = map (normsqr . foo) [(_1, _R1), (_2, _R2), (_3, _R3), (_4, _R4)]
  where
    foo (f1, f2) = (m &- (f1 value) *& idmtx) *. (f2 $ transpose vec)
    (vec, value) = symmEigen m

spec :: Spec
spec = do
    describe "QR Decomposition" $ do
        it "3x3 QR residual is small" $ do
            let res = testQR testData1
            isMainlyZero (frobeniusNorm res :: Double) `shouldBe` True
        it "2x2 QR residual is small" $ do
            let res = testQR testData5
            isMainlyZero (frobeniusNorm res :: Double) `shouldBe` True

    describe "Eigenvalue Decomposition" $ do
        it "3x3 Eigen residuals are small" $ do
            let res = testEigen3 testData2
            all isMainlyZero res `shouldBe` True
        it "4x4 Eigen residuals are small" $ do
            let res = testEigen4 testData3
            all isMainlyZero res `shouldBe` True
