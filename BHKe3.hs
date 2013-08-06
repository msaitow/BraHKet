
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
    
    p = QIndex "P" Generic True      
    q = QIndex "Q" Generic True           
    r = QIndex "R" Generic True      
    s = QIndex "S" Generic True           

    -- Creaton and annihilation operators
    -- t  = baseTensor "@T2" [c,d,k,l] [[0,1,2,3],[1,0,3,2]] True
    eL = baseSFGen [i, j, a, b]
    eH = baseSFGen [p, q, r, s]
    eR = baseSFGen [c, d, k, l]
    
    term = baseTerm 1 [] [eL, eH, eR]
    ordered = normalOrderE term
    vev     = takeVEV ordered
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
  print $ "length   : " ++ (show $ length ordered)
  print $ "vev      : " ++ (show $ vev)
  print $ "length   : " ++ (show $ length vev)
