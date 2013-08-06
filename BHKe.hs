
import Data.List
import Data.Maybe
import BraHKet.Core
       
main = do
  let
    -- Declaration for indices
    i = QIndex "I" Active  False
    a = QIndex "A" Virtual False
    j = QIndex "J" Active  True
    b = QIndex "B" Virtual True      
    k = QIndex "K" Active  True      
    c = QIndex "C" Virtual True           
    l = QIndex "L" Active  True      
    d = QIndex "D" Virtual True           

    -- Creaton and annihilation operators
    t  = baseTensor "@T2" [c,d,k,l] [[0,1,2,3],[1,0,3,2]] True
    eL = baseSFGen [i, a]
    eR = baseSFGen [c, d, k, l]
    
    term = baseTerm 1 [] [eL, eR, t]
    ordered = normalOrderE term
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
