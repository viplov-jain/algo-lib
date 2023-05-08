module Common.Primes (primes, factors_pp_arr, factors_arr) where

import Common.Types

import Data.Array

union :: IntType t => [t] -> [t] -> [t]
union (x:xs) (y:ys) = case (compare x y) of
           LT -> x : union  xs  (y:ys)
           EQ -> x : union  xs     ys
           GT -> y : union (x:xs)  ys
union  xs     []    = xs
union  []     ys    = ys

infCompose :: (t -> t) -> t
infCompose g = g (infCompose g)

gapsW :: IntType a => a -> [a] -> [a] -> [a]
gapsW k (d:w) (c:cs) =
  if k < c then k : gapsW (k+d) w (c:cs) else gapsW (k+d) w cs

hitsW :: IntType a => a -> [a] -> [a] -> [[a]]
hitsW k (d:w) (p:ps) =
  if k < p
  then hitsW (k+d) w (p:ps)
  else scanl (\c d->c+p*d) (p*p) (d:w) : hitsW (k+d) w ps

joinT :: IntType t => [[t]] -> [t]
joinT ((x:xs):t) = x : (union xs . joinT . pairs) t
joinT [] = []

pairs :: IntType t => [[t]] -> [[t]]
pairs (xs:ys:t) = union xs ys : pairs t
pairs t = t

wheel :: IntType a => [a]
wheel = 2:4:2:4:6:2:6:4:2:4:6:6:2:6:4:2:6:4:6:8:4:2:4:2:
        4:8:6:4:6:2:4:6:2:6:6:4:2:4:6:2:6:4:2:4:2:10:2:10:wheel

primes :: IntType a => [a]
primes = 2:3:5:7:(infCompose ((11:) . tail . gapsW 11 wheel . joinT . hitsW 11 wheel))

upper_limit_p :: Int
upper_limit_p = 10000

factor_arr :: Array Int Int
factor_arr =
  loop primes (array (2, upper_limit_p) [(i, i) | i <- [2..upper_limit_p]])
  where
    loop (p:pt) arr = if p > 100 then arr else
      loop pt (arr//[(i, p) | i <- [p*p, p*(p+1)..upper_limit_p]])

{- A_i = (i0, p, k), i = i0*p^k-}
factors_pp_arr :: Array Int (Int, Int, Int)
factors_pp_arr =
  array (2, upper_limit_p) [(i, redi i) | i <- [2..upper_limit_p]]
  where
    redi i = reduce (i, factor_arr!i, 0)
    reduce :: (Int, Int, Int) -> (Int, Int, Int)
    reduce (n, p, acc) = let (dv, md) = divMod n p in
      if md == 0 then reduce (dv, p, acc+1) else (n, p, acc)

factors_arr :: Array Int [(Int, Int)]
factors_arr =
  array (2, upper_limit_p) [(i, factorize i []) | i <- [2..upper_limit_p]]
  where
    factorize :: Int -> [(Int, Int)] -> [(Int, Int)]
    factorize 1 l = l
    factorize i l = let (i0, p, k) = factors_pp_arr!i in factorize i0 ((p, k):l)
