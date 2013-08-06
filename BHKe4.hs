
import Data.List
import Data.Maybe
import BraHKet.Core
       
main = do
  let
    -- Declaration for indices
    i = QIndex "I" Active  False
    j = QIndex "J" Active  False
    a = QIndex "A" Virtual False
    b = QIndex "B" Virtual False      

    k = QIndex "K" Active  True
    l = QIndex "L" Active  True      
    c = QIndex "C" Virtual True           
    d = QIndex "D" Virtual True
    
    p = QIndex "P" Generic True      
    q = QIndex "Q" Generic True           
    r = QIndex "R" Generic True      
    s = QIndex "S" Generic True           

    -- Creaton and annihilation operators
    t  = baseTensor "@T2" [c,d,k,l] [[0,1,2,3],[1,0,3,2]] True
    h  = baseOne [p, q]
    eL = baseSFGen [i, j, a, b]
    eH = baseSFGen [p, q]
    eR = baseSFGen [c, d, k, l]
    
    term = baseTerm 1 [] [eL, h, eH, t, eR]
    ordered = normalOrderE term
    vev     = takeVEV ordered
    survived = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) vev
    maximums = fmap (maximum . generateAllConfs) survived
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
  print $ "length   : " ++ (show $ length ordered)
  print $ "vev      : " ++ (show $ vev)
  print $ "length   : " ++ (show $ length vev)
  print $ "Killed   : " ++ (show $ survived)
  print $ "maximum  : " ++ (show $ maximums)
  print $ "combined : " ++ (show $ combineTerms survived)
