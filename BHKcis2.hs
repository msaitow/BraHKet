
import Data.List
import Data.Maybe
import BraHKet.Core

-- -- Implimentation1 (still doesn't work)
-- combineTerms :: [QTerm] -> [QTerm]
-- --combineTerms terms = head . foldl (combine origTerms) [] . maxTerms
-- combineTerms terms = combine origTerms maxTerms
--   where
--     maxTerms = fmap (maximum . generateAllConfs) terms
--     origTerms = nub maxTerms
-- 
--     combine :: [QTerm] -> [QTerm] -> [QTerm]
--     combine [] _ = []
--     combine (x:xs) arera = catchUpTerms xs $ catchUpTerms x arera 
-- 
--     catchUpTerms :: [QTerm] -> [QTerm] -> [QTerm]
--     catchUpTerms kore korera = if kore `elem` korera then (baseTerm newNum (tCoeff kore) (tTensor kore)) : diff else korera
--       where
--         same = filter (==kore) korera
--         diff = filter (/=kore) korera
--         newNum = sum $ fmap (tNum) same

 -- combineTerms :: [QTerm] -> [QTerm]
 -- --combineTerms terms = head . foldl (combine origTerms) [] . maxTerms
 -- combineTerms terms = fmap (changeFactors) origTerms
 --   where
 --     maxTerms = fmap (maximum . generateAllConfs) terms
 --     origTerms = nub maxTerms
 -- 
 --     changeFactors :: QTerm -> QTerm
 --     changeFactors are = baseTerm myFactor (tCoeff are) (tTensor are)
 --       where myFactor = collectFactors are maxTerms
 -- 
 --     collectFactors :: QTerm -> [QTerm] -> Double
 --     collectFactors kore korera = 
 --       let
 --         arera = [x | x <- korera, isAdditive kore x]
 --       in sum $ fmap (tNum) arera
       
main = do
  let
    -- DEclaration for indices
    i = QIndex "I" Core    False
    a = QIndex "A" Virtual False
    j = QIndex "J" Core    True
    b = QIndex "B" Virtual True      
    p = QIndex "P" Generic True      
    q = QIndex "Q" Generic True           
    r = QIndex "R" Generic True      
    s = QIndex "S" Generic True           

    -- Creaton and annihilation operators
    eL = baseSDGen [i, a]
    eH = baseSDGen [p, q, r, s]
    eR = baseSDGen [b, j]
    
    -- Declaraton for tensors
    v = baseERI [p, q, r, s]
    t = baseTensor "T1" [b, j] [[0,1],[1,0]] True
    term = baseTerm 0.5 [] [v, eL, eH, eR, t]
    --term = baseTerm 1 [] [eL, eH, eR]
    ordered = normalOrderG term
    killed = fmap (killKDeltas) ordered
    survived = fmap (fromJust) $ filter (\x -> x /= Nothing) $ killed
    maximums = fmap (maximum . generateAllConfs) survived
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
  print $ "killedKD : " ++ (show $ killed)
  print $ "survived : " ++ (show $ survived)
  print $ "maximum  : " ++ (show $ maximums)
  print $ "combined : " ++ (show $ combineTerms survived)
  print $ "num      : " ++ (show $ nub maximums)
  -- print $ "Term3    : " ++ (show $ getSummedBody g)
