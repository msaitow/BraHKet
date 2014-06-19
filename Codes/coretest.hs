
import Data.List
import Data.Maybe
import BraHKet

-------------------------------------------------------------------------
-- Test code for contracting the core indices in the spin-free RDMs
-------------------------------------------------------------------------

main = do
  let
    -------------------------------------------------    
    -- Declaration for molecular orbital indices
    -------------------------------------------------        
    i = sfIndex "i" Active NonDummy
    j = sfIndex "j" Active NonDummy
    k = sfIndex "k" Active NonDummy
    l = sfIndex "l" Active NonDummy          
    m = sfIndex "m" Active NonDummy
    n = sfIndex "n" Active NonDummy    

    a = sfIndex "a" Core NonDummy
    b = sfIndex "b" Core NonDummy
    c = sfIndex "c" Core NonDummy
    d = sfIndex "d" Core NonDummy    
    
    v = sfIndex "v" Core NonDummy
    w = sfIndex "w" Core NonDummy
    x = sfIndex "x" Core NonDummy
    y = sfIndex "y" Core NonDummy    
    z = sfIndex "z" Core NonDummy          

    -------------------------------------------------    
    -- Definitions of tensors
    -------------------------------------------------    
    rdm1 = baseSFRDM [v, z, i, j, k,
                      m, w, l, x, n]

    thisTerm1   = baseTerm 1.0 [] [rdm1]
    contracted1 = contractCoreSF thisTerm1
    -------------------------------------------------
    rdm2 = baseSFRDM [v, w, i, j,
                      k, l, y, z]

    thisTerm2   = baseTerm 1.0 [] [rdm2]
    contracted2 = contractCoreSF thisTerm2
    -------------------------------------------------
    rdm3 = baseSFRDM [v, y, i, j,
                      k, w, l, x]

    thisTerm3   = baseTerm 1.0 [] [rdm3]
    contracted3 = contractCoreSF thisTerm3
    -------------------------------------------------
    rdm4 = baseSFRDM [c, v, i, j,
                      l, k, y, z]

    thisTerm4   = baseTerm 1.0 [] [rdm4]
    contracted4 = contractCoreSF thisTerm4
    -------------------------------------------------
    rdm5 = baseSFRDM [a, v, y, i, j,
                      b, k, w, l, x]

    thisTerm5   = baseTerm 1.0 [] [rdm5]
    contracted5 = contractCoreSF thisTerm5
    -------------------------------------------------
    rdm6 = baseSFRDM [v, y, i, j,
                      k, w, l, x]

    thisTerm6   = baseTerm 1.0 [] [rdm6]
    contracted6 = contractCoreSF thisTerm6
    -------------------------------------------------
    rdm7 = baseSFRDM [c, v, w,
                      z, y, x]

    thisTerm7   = baseTerm 1.0 [] [rdm7]
    contracted7 = contractCoreSF thisTerm7
    -------------------------------------------------

  print "---------------------------------------------------------"
  print $ "inTerm1  :" ++ (show thisTerm1)
  print $ "outTerm1 :" ++ (show contracted1)
  print "---------------------------------------------------------"  
  print $ "inTerm2  :" ++ (show thisTerm2)
  print $ "outTerm2 :" ++ (show contracted2)
  print "---------------------------------------------------------"  
  print $ "inTerm3  :" ++ (show thisTerm3)
  print $ "outTerm3 :" ++ (show contracted3)
  print "---------------------------------------------------------"  
  print $ "inTerm4  :" ++ (show thisTerm4)
  print $ "outTerm4 :" ++ (show contracted4)
  print "---------------------------------------------------------"  
  print $ "inTerm5  :" ++ (show thisTerm5)
  print $ "outTerm5 :" ++ (show contracted5)
  print "---------------------------------------------------------"  
  print $ "inTerm6  :" ++ (show thisTerm6)
  print $ "outTerm6 :" ++ (show contracted6)
  print "---------------------------------------------------------"  
  print $ "inTerm7  :" ++ (show thisTerm7)
  print $ "outTerm7 :" ++ (show contracted7)
  print "---------------------------------------------------------"  
