module FFT where

import Data.Complex
import Data.Array

-- represent Vectors as a list of coefficients from lowest degree to highest
type Vector = [Float]

extend :: Int -> Vector -> Vector
extend n p = p ++ take (n - length p) (repeat 0)

fft :: Vector -> Vector -> Vector
-- compute the convolution of two Vectors in O(n log n) time
-- using Fast Fourier Transform
fft a b
  | length a /= length b =
    let maxLen = max (length a) (length b)
    in fft (extend maxLen a) (extend maxLen b)
  | length a `mod` 2 == 1 =
    let newLen = length a + 1
    in fft (extend newLen a) (extend newLen b)
  | otherwise = map iThCoefficient [0 .. 2 * n - 1]
  where
    n = length a
    aRoots = evalRootsOfUnity a
    bRoots = evalRootsOfUnity b
    roots = zipWith (*) aRoots bRoots

    d :: Int -> Float
    d s =
      realPart $ sum $ mapInd (\rootVal j -> rootVal * rootOfUnity (2 * n) (s * j)) roots

    iThCoefficient i = 1 / (fromIntegral $ 2 * n) * d (2 * n - i)

mapInd :: (a -> Int -> b) -> [a] -> [b]
mapInd f l = zipWith f l [0..]

evalRootsOfUnity :: Vector -> [Complex Float]
-- evaluate a Vector at the 2n roots of unity
evalRootsOfUnity [] = []
evalRootsOfUnity [c] = take 2 $ repeat $ c :+ 0
evalRootsOfUnity p = map evalN [0 .. (2 * n - 1)]
  where
    n = length p

    (evens, odds) = splitEvensOdds p
    oddRoots = listArray' $ evalRootsOfUnity odds
    evenRoots = listArray' $ evalRootsOfUnity evens
    roots = rootsOfUnity (2 * n)
    evalN :: Int -> Complex Float
    -- evaluate at the nth root of unity using the fact that
    -- A(x) = A_even(x^2) + x * A_odd(x^2) and for roots of unity, squaring them
    -- returns the (n * 2)th root of unity
    evalN i = evenRoots!d + roots!i * oddRoots!d
      where d = i `mod` n


splitEvensOdds :: Vector -> (Vector, Vector)
-- splitEvensOdds [] = ([], [])
-- splitEvensOdds (even':xs) = (even':otherEvens, odds)
--   where (odds, otherEvens) = splitEvensOdds xs
splitEvensOdds = foldl (\(odds, evens) x -> (x:evens, odds)) ([], [])

rootsOfUnity :: Int -> Array Int (Complex Float)
rootsOfUnity n = listArray' $ map (rootOfUnity n) [0 .. n - 1]

rootOfUnity :: Int -> Int -> Complex Float
rootOfUnity n i = mkPolar 1 (2 * i' * pi / n')
  where
    i' = fromIntegral (i `mod` n)
    n' = fromIntegral n

listArray' :: [a] -> Array Int a
listArray' xs = listArray (0, length xs - 1) xs
