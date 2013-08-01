
import Data.List
import Data.Maybe
import BraHKet.Core

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
    -- survived = filter (\x -> x /= Nothing) $ fmap (\x -> if x /= Nothing then fromJust x else id x) killed
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
  print $ "killedKD : " ++ (show $ killed)
  -- print $ "survived : " ++ (show $ survived)
  -- print $ "Term3    : " ++ (show $ getSummedBody g)
