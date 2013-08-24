
import Data.List
import Data.Maybe
import BraHKet.Core

------------------------------------------------------------------------------------------------------------------------
-- It derives the explicit formulas for the external/external excitation block of the IC-MRCI sigma-equation
-- where the non-relativistic electronic Hamiltonian is given as
--
--     H = h^p_q E^q_p + 1/2 V^{pq}_{rs} E^{rs}_{pq}.
--
-- Also E represents the spin-free excitation operator that serves as the generator of the
-- unitary group
--
--     E^p_q       = \sum_{\sigma=\alpha,\beta} a^{\dagger}_{p\sigma} a_{q\sigma}
--
-- and
--
--     E^{pq}_{rs} = \sum_{\sigma,\tau=\alpha,\beta} a^{\dagger}_{p\sigma} a^{\dagger}_{q\tau} a_{s\tau} a_{r\sigma}. 
--
-- To compile:
--
--     shell> ghc (-O2) ./THIS_FILE.hs
--
-- For more details, see Journal of Chemical Physics 139, 044118(2013).
------------------------------------------------------------------------------------------------------------------------

main = do
  let
    -- Declaration for molecular orbital indices
    i = QIndex "i" Active  NonDummy
    a = QIndex "a" Virtual NonDummy
    j = QIndex "j" Active  NonDummy
    b = QIndex "b" Virtual NonDummy      

    k = QIndex "k" Active  Dummy      
    c = QIndex "c" Virtual Dummy           
    l = QIndex "l" Active  Dummy      
    d = QIndex "d" Virtual Dummy
    
    p = QIndex "p" Generic Dummy      
    q = QIndex "q" Generic Dummy           
    r = QIndex "r" Generic Dummy      
    s = QIndex "s" Generic Dummy           

    -- Definitions of tensors
    t   = baseTensor "@T2" [c,d,k,l] [[0,1,2,3],[1,0,3,2]] Classical -- t^{cd}_{kl}
    h   = baseOne   [p, q]       -- h^{p}_{q}
    v   = baseERI   [p, q, r, s] -- v^{pq}_{rs}
    eL  = baseSFGen [i, j, a, b] -- E^{ij}_{ab}
    eH1 = baseSFGen [p, q]       -- E^{p}_{r}
    eH2 = baseSFGen [p, q, r, s] -- E^{pq}_{rs}
    eR  = baseSFGen [c, d, k, l] -- E^{cd}_{kl}

    -----------------------------------------------------------------------------
         --           
         --  Define: E^{ij}_{ab} h^{p}_{q} E^{p}_{q} T^{cd}_{kl} E^{cd}_{kl}
         --
    -----------------------------------------------------------------------------          
    term1       = baseTerm 1.0 [] [eL, h, eH1, t, eR]
    ordered1    = normalOrderE term1
    vev1        = takeVEV ordered1
    survived1   = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) vev1
    combined1   = combineTerms survived1    
    decomposed1 = concat $ fmap contractCoreSF $ fmap masquerade $ concat $ fmap (generateInteractions [Core, Active]) $ combined1
    dsurvived1  = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) decomposed1 
    dcombined1  = combineTerms dsurvived1
                    
    -----------------------------------------------------------------------------    
         --           1
         --  Define: --- E^{ij}_{ab} V^{pq}_{rs} E^{pq}_{rs} T^{cd}_{kl} E^{cd}_{kl}
         --           2
    -----------------------------------------------------------------------------                
    term2       = baseTerm 0.5 [] [eL, v, eH2, t, eR]
    ordered2    = normalOrderE term2
    vev2        = takeVEV ordered2
    combined2   = combineTerms survived2    
    survived2   = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) vev2
    decomposed2 = concat $ fmap contractCoreSF $ fmap masquerade $ concat $ fmap (generateInteractions [Core, Active]) $ combined2
    dsurvived2  = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) decomposed2
    dcombined2  = combineTerms dsurvived2
    
  print $ "inTerm_h1    : " ++ (show term1)
  print $ "inTerm_v2    : " ++ (show term2)
  print $ "outTerm_h1   : " ++ (show ordered1)
  print $ "outTerm_v2   : " ++ (show ordered2)
  print $ "vev_h1       : " ++ (show vev1)
  print $ "vev_v2       : " ++ (show vev2)  
  print $ "combined_h1  : " ++ (show combined1)
  print $ "combined_v2  : " ++ (show combined2)  
  print $ "dsurvived1   : " ++ (show dsurvived1)
  print $ "dsurvived2   : " ++ (show dsurvived2)  
  print $ "dcombined_h1 : " ++ (show dcombined1)
  print $ "dcombined_v2 : " ++ (show dcombined2)
