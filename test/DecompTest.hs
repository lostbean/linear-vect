module Main where

import Linear.Class
import Linear.Decomp
import Linear.Vect

-- ================================== Test Data ==================================

testData1 :: Mat3 Double -- Source <http://en.wikipedia.org/wiki/QR_decomposition>
testData1 = Mat3 (Vec3 12 (-51) 4) (Vec3 6 167 (-68)) (Vec3 (-4) 24 (-41))

testData2 :: Mat3 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData2 = Mat3 (Vec3 2 1 0) (Vec3 1 3 (-1)) (Vec3 0 (-1) 6)

testData3 :: Mat4 Double
testData3 = Mat4 (Vec4 0 10 3 9) (Vec4 10 12 6 15) (Vec4 3 6 0 7) (Vec4 9 15 7 8)

testData4 :: Mat4 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData4 = Mat4 (Vec4 4 1 (-1) 2) (Vec4 1 4 1 (-1)) (Vec4 (-1) 1 4 1) (Vec4 2 (-1) 1 4)

testData5 :: Mat2 Double -- Source <Orthogonal Bases and the QR Algorithm> <by Peter J. Olver>
testData5 = Mat2 (Vec2 2 1) (Vec2 1 3)

testQR :: (Transpose g g, MultSemiGroup g, AbelianGroup g, OrthoMatrix g, SquareMatrix g) => g -> g
testQR m = m &- (q .*. r)
  where
    (q, r) = qrHouse m

testEigen ::
    ( DotProd b v1
    , HasOne v2
    , HasOne v3
    , HasTwo v2
    , HasTwo v3
    , HasThree v2
    , HasThree v3
    , HasFour v2
    , HasFour v3
    , Foldable v2
    , Foldable v3
    , Functor v2
    , Functor v3
    , OrthoMatrix (v3 (v1 b))
    , MultSemiGroup (v3 (v1 b))
    , Diagonal (v2 (v1 b)) (v3 (v1 b))
    , Fractional (v1 b)
    , Ord (v1 b)
    , NearZero (v1 b)
    , LeftModule (v3 (v1 b)) (v1 b)
    , LinearMap (v1 b) v3
    , SquareMatrix (v3 (v1 b))
    , Transpose (v3 (v1 b)) (v3 (v1 b))
    , Dimension (v3 (v1 b))
    ) =>
    v3 (v1 b) ->
    [b]
testEigen m = map (normsqr . foo) $ take n [(_1, _1), (_2, _2), (_3, _3), (_4, _4)]
  where
    n = dim m
    foo (f1, f2) = (m &- (f1 value) *& idmtx) *. (f2 $ transpose vec)
    (vec, value) = symmEigen m

-- ================================== Main ==================================

main :: IO ()
main = do
    let qrResid3 = testQR testData1
    putStrLn $ "QR residual (3x3): " ++ show qrResid3
    let qrResid2 = testQR testData5
    putStrLn $ "QR residual (2x2): " ++ show qrResid2
    let eigenResid3 = testEigen testData2
    putStrLn $ "Eigen residual (3x3): " ++ show eigenResid3
    let eigenResid4 = testEigen testData3
    putStrLn $ "Eigen residual (4x4): " ++ show eigenResid4
    putStrLn "All decomposition tests passed."
