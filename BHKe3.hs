
import Data.List
import Data.Maybe
import BraHKet.Core
       
main = do
  let
    -- Declaration for indices
    i = QIndex "i" Active  False
    a = QIndex "a" Virtual False
    j = QIndex "j" Active  False
    b = QIndex "b" Virtual False      

    k = QIndex "k" Active  True      
    c = QIndex "c" Virtual True           
    l = QIndex "l" Active  True      
    d = QIndex "d" Virtual True
    
    p = QIndex "p" Generic True      
    q = QIndex "q" Generic True           
    r = QIndex "r" Generic True      
    s = QIndex "s" Generic True           

    -- Definition of tensors
    t  = baseTensor "@T2" [c,d,k,l] [[0,1,2,3],[1,0,3,2]] True
    v  = baseERI [p, q, r, s]
    eL = baseSFGen [i, j, a, b] -- E^{ij}_{ab}
    eH = baseSFGen [p, q, r, s] -- E^{pq}_{rs}
    eR = baseSFGen [c, d, k, l] -- E^{cd}_{kl}

         --           1
         --  Define: --- E^{ij}_{ab} V^{pq}_{rs} E^{pq}_{rs} T^{cd}_{kl} E^{cd}_{kl}
         --           2
    term = baseTerm 0.5 [] [eL, v, eH, t, eR]
    ordered = normalOrderE term
    vev     = takeVEV ordered
    survived = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) vev
  print $ "inTerm   : " ++ (show term)
  print $ "outTerm  : " ++ (show $ ordered)
  print $ "length   : " ++ (show $ length ordered)
  print $ "vev      : " ++ (show $ vev)
  print $ "length   : " ++ (show $ length vev)
  print $ "combined : " ++ (show $ combineTerms survived)
