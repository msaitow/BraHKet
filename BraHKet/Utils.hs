
module BraHKet.Utils
(
  makeFlat,      -- Returns flattered list of pairs
  makeFlat2,     -- Returns flattered list of pairs
  allDiff,       -- Returns true if all the elements are different from the others
  binCombi,      -- Returns a binomial combination
  permuteSign,   -- Returns a sign of permutations
  repeatN,       -- Repeat application of function f into a for N times
  makeCombi,     -- Small function to make combinations of all the object given as a list 
  uniPerm,       -- Small functions to perform permutation operation
  uniPermR,      -- Small functions to perform permutation operation
  makePerms,     -- Returns all permutations from given irreducible permutations
  makePermSFGen, -- Returns permutations for spin-free generators
  makePermSFRDM  -- Returns permutations for spin-free density-matrices
) where

import Data.List
import qualified Data.Maybe as Maybe


-- Make flattered list like
-- [(i1,e1),(i2,e2)..(in,en)] => [i1,i2..in,e1,e2..en]
makeFlat2 :: (Ord a) => [(a, a)] -> [a]
makeFlat2 pairs = (fmap fst pairs) ++ (fmap snd pairs)

-- Make flattered list of pairs
makeFlat :: (Ord a) => [(a, a)] -> [a]
makeFlat [x] = (fst x) : (snd x) : []
makeFlat (x:xs) = (fst x) : (snd x) : makeFlat xs

-- Returns True if all the elements are different
allDiff :: (Ord a) => [a] -> Bool
allDiff []     = error "allDiff: Undefined input"
allDiff [x]    = True
allDiff (x:xs) = (x `elem` xs) && allDiff xs

-- Returns binomial combination whose length is nCm
binCombi :: (Ord a) => Int -> [a] -> [[a]]
binCombi 0 _      = [[]]
binCombi _ []     = []
binCombi n (x:xs) = (map (x:) (binCombi (n-1) xs)) ++ binCombi n xs

-- Returns sign of permutations in between two lists
permuteSign :: (Ord a) => [a] -> [a] -> Int
permuteSign initList finalList = intArraySign relPerm
  where
    relPerm  = fmap posNumber initList
    posNumber kore
      | kore `elem` finalList = Maybe.fromJust $ findElem kore finalList 0
      | otherwise = error "permuteSign: Algorithmic error"
    findElem :: (Ord a) => a -> [a] -> Int -> Maybe Int
    findElem are (x:xs) len
      | length xs == 0 = if are == x then Just len else Nothing
      | otherwise      = if are == x then Just len else findElem are xs (len+1) 

-- Returns sign of the corresponding permutations
intArraySign :: [Int] -> Int
intArraySign xs = if (mod (oversetNum xs) 2) > 0 then -1 else 1
  where
    oversetNum (x:xs)
      | xs == []  = 0
      | otherwise = sum (map (\y -> if x > y then 1 else 0) xs) + oversetNum xs

-- Repeat application of function f for n times 
repeatN :: (a -> a) -> Int -> a -> a
repeatN f 1 x = f x
repeatN f n x = f $ repeatN f (n-1) x

-- Returns a list of object composed of combinations of all the object in the initial list
baseCombi :: (Ord a) => [[a]] -> Int -> [[a]]
baseCombi xs count
  | count /= 0 = [x:y    | x <- (xs !! ((length xs)-1)), y <- (baseCombi (init xs) (count-1))]
  | otherwise  = [x:y:[] | x <- (xs !! 0), y <- (xs !! 1)]                 

-- Interface of the above function
makeCombi :: (Ord a) => [[a]] -> [[a]]
makeCombi xs = baseCombi xs ((length xs)-2) 

-- Permute list of object according to the given permutation
-- Make sure if the permutation is not a contiguous numbers, everything becomes crazy.
uniPerm :: (Ord a) => [Int] -> [a] -> [a]
uniPerm (x:xs) array
  | length array < length xs = error "An invalid argument is given"
  | length xs == 0 = [array !! x]
  | otherwise = (array !! x) : uniPerm xs array

-- Small function that takes argumentes in reverse order to uniPerm
uniPermR :: (Ord a) => [a] -> [Int] -> [a]
uniPermR array xs = uniPerm xs array

-- Make all permutations from the irreducible groups of permutations
makePerms :: [[Int]] -> [[Int]]
makePerms a = genAllPerms a a
  where
    
    -- Body of the makePerms
    genAllPerms :: [[Int]] -> [[Int]] -> [[Int]]
    genAllPerms irrSource inSource
      | length outSource == length inSource = outSource
      | otherwise = genAllPerms irrSource outSource
      where
        outSource = nub $ genPerm irrSource inSource

        -- Small function which is necessary for makePerms
        genPerm :: [[Int]] -> [[Int]] -> [[Int]]
        genPerm (p:ps) inArrays
          | length ps == 0 = pArrays
          | otherwise      = pArrays ++ genPerm ps inArrays
          where pArrays = fmap (uniPerm p) inArrays
        

-- Make all permutations for the spin-free unitary group generator
makePermSFGen :: Int -> [[Int]]
makePermSFGen order =
  let base  = permutations [0..order-1]
      tails = fmap (fmap (\x -> x + order)) base
  in zipWith (++) base tails

-- Make all permutations for the spin-free reduced-density matrices
makePermSFRDM :: Int -> [[Int]]
makePermSFRDM order =
  let base  = permutations [0..order-1]
      tails = fmap (fmap (\x -> x + order)) base
  in (zipWith (++) tails base) ++ makePermSFGen order
