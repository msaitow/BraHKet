
import Interface.SFIndex
import Interface.SFOperator
import BraHKet
import Data.Maybe

----------------------------------------------------------------------------------------------------
-- A bunch of code that derive the diagonal-preconditioner and sigma vector of the fic-MRCI
-- To compile: shell> ghc (-O2) -i:"../../" ./icciEngine.hs
----------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------
-- Print contraction
----------------------------------------------------------------------------------------------------
printContra :: (QTensor, QTerms) -> String
printContra myBlock = "  " ++ (show sigma) ++ " <-- \n" ++ sigmaPart ++ "\n    >> numTerm >> " ++ (show numTerm) ++ "\n\n\n"
  where
    numTerm = length $ snd myBlock
    sigma   = fst myBlock
    mySigma = snd myBlock

    sigmaPart = foldl (++) "" $ fmap (\x -> "    " ++ (show x) ++ "\n") mySigma

----------------------------------------------------------------------------------------------------
-- A function to generate all the icci equations
----------------------------------------------------------------------------------------------------    
sigmaEngine :: QIndices -> QIndices -> (QTensor, QTerms) -- :: (sigma, contraction terms)
sigmaEngine leftOps rightOps
  | length actives  == 6                         = (sig2, res1) 
  | length leftOps  == 0 && length rightOps == 0 = (sig0, res2) 
  | length leftOps  == 0                         = (sig0, res3) 
  | length rightOps == 0                         = (sig2, res4) 
  | otherwise                                    = (sig2, res5) 
  where

    ----------------------------------------------------------------------------------
    res1 = symmSigmaComm leftOps rightOps  -- Commutator block
    res2 = [baseTerm 1.0 (["E0"]) []]                               -- Ref/Ref block
    res3 = refICblock rightOps                                      -- Ref/ICB block
    res4 = icRefblock leftOps                                       -- ICB/Ref block
    res5 = constructSigma leftOps rightOps -- ICB/ICB block
    ----------------------------------------------------------------------------------
    sig0 = baseTensor "@S0" leftOps []     Classical
    sig2 = baseTensor "@S2" leftOps mySymm Classical
    ----------------------------------------------------------------------------------

    -- Only for the ccaa/ccaa, aavv/aavv and ccvv/ccvv elements, the two-fold symmetry of the sigma-vector is allowed
    mySymm
      | length leftOps < 4 || length rightOps < 4 = [[0,1,2,3]]
      | otherwise                                 = if lSymm && rSymm then [[0,1,2,3],[1,0,3,2]] else [[0,1,2,3]]
      where
        lSymm = iSpace ( leftOps !! 0) == iSpace ( leftOps !! 1) && iSpace ( leftOps !! 2) == iSpace ( leftOps !! 3)
        rSymm = iSpace (rightOps !! 0) == iSpace (rightOps !! 1) && iSpace (rightOps !! 2) == iSpace (rightOps !! 3)        

    tOps     = transposeIndices rightOps 
    myAmp    = baseTensor "@T2" tOps [[0,1,2,3],[1,0,3,2]] Classical
    actives  = filter (\x -> x == Active) $ fmap iSpace (leftOps ++ rightOps)

    h  = baseOne   [g0, g1]         -- h^{g0}_{g1}
    v  = baseERI   [g0, g1, g2, g3] -- v^{g0g1}_{g2g3}
    eh = baseSFGen [g0, g1]         -- E^{g0}_{g1}
    ev = baseSFGen [g0, g1, g2, g3] -- E^{g0g1}_{g2g3}

    -- s0 <-- <0| H Tr Er |0> 
    refICblock :: QIndices -> QTerms
    refICblock myInds = constructCoreFock $ processTermsE [Core, Active, Virtual] [baseTerm 1.0 [] [h, eh, myAmp, eR], baseTerm 0.5 [] [v, ev, myAmp, eR]]
      where
        eR = baseSFGen myInds

    -- s2 <-- <0| El H T0 |0>
    icRefblock myInds = constructCoreFock $ processTermsE [Core, Active, Virtual] [baseTerm 1.0 ["T0"] [eL, h, eh], baseTerm 0.5 ["T0"] [eL, v, ev]]
      where
        eL = baseSFGen myInds

    -- s2 <-- Tr 1/2(1+P_{LR}) <0| [ El^t, [ H, Er ] ] |0> + E0 <0| El^t Er |0>
    symmSigmaComm :: QIndices -> QIndices -> QTerms
    symmSigmaComm lOps rOps = screenTerms $ combineTerms $ fmap (scaleTerm 0.5) $ fockTerms1 ++ fockTerms2
      where
        tlOps = transposeIndices rOps 
        trOps = transposeIndices lOps 
        fockTerms1 = constructSigmaCommBase  lOps  rOps tOps
        fockTerms2 = constructSigmaCommBase tlOps trOps tOps

------------------------------------------------------------
-- A function to generate all the icci preconditioner
------------------------------------------------------------
precondEngine :: QIndices -> (QTensor, QTerms)
precondEngine leftOps
  | length actives == 3 = (h2, res1)
  | length leftOps == 0 = (h0, res2)
  | otherwise           = (h2, res3)
  where
    actives  = filter (\x -> x == Active) $ fmap iSpace leftOps
    
    --------------------------------------------------------------------
    res1 = constructHamiltonianComm leftOps (transposeIndices leftOps)
    res2 = [baseTerm 1.0 (["1/E0"]) []]
    res3 = constructHamiltonian     leftOps (transposeIndices leftOps)
    --------------------------------------------------------------------
    h0   = baseTensor "@H0" leftOps []                    Classical
    h2   = baseTensor "@H2" leftOps [[0,1,2,3],[1,0,3,2]] Classical
    --------------------------------------------------------------------

--------------------------------------
-- Left-side excitation operators
--------------------------------------
eL = [
  [],           -- unity 
  [p, q, a, b], -- (a,a) -> (v,v)
  [p, q, r, a], -- (a,a) -> (a,v)
  [i, j, a, b], -- (c,c) -> (v,v) 
  [i, p, a, b], -- (c,a) -> (v,v)
  [i, j, p, q], -- (c,c) -> (a,a)
  [i, j, p, a], -- (c,c) -> (a,v)
  [i, p, q, r], -- (c,a) -> (a,a)
  [i, p, q, a], -- (c,a) -> (a,v)
  [p, i, q, a]  -- (a,c) -> (a,v)
  ]

--------------------------------------
-- Right-side excitation operators
--------------------------------------
eR = [
  [],           -- unity        
  [v0, v1, a0, a1], -- (a,a) -> (v,v)
  [a0, v0, a1, a2], -- (a,a) -> (a,v)
  [v0, v1, c0, c1], -- (c,c) -> (v,v)
  [v0, v1, c0, a0], -- (c,a) -> (v,v)
  [a0, a1, c0, c1], -- (c,c) -> (a,a)
  [a0, v0, c0, c1], -- (c,c) -> (a,v)
  [a0, a1, c0, a2], -- (c,a) -> (a,a)
  [a0, v0, c0, a1], -- (c,a) -> (a,v)
  [a0, v0, a1, c0]  -- (a,c) -> (a,v)
  ]
--------------------------------------

---------------------------------------------------
-- Write the contractions into the tensor file
---------------------------------------------------
writeContractions :: (FilePath, String) -> IO ()
writeContractions (fileName, contractions) = do
  let dirName = "tensors/"
  writeFile (dirName ++ fileName) contractions
  

main = do
  -- Construct sigma-equation 
  let sigmaElements = [ ( (a,b), sigmaEngine a b ) | a <- eL, b <- eR]
  -- Construct preconditioner
  let precondElements = fmap precondEngine eL

  let printedData1 = foldl (++) "" $ fmap printContra $ fmap snd sigmaElements
  let printedData2 = foldl (++ ) "" $ fmap printContra $ precondElements
    
  putStr "\n"
  putStr $ ">>>>>>> Sigma-equation <<<<<<<\n\n"
  putStr $ printedData1 ++ "\n\n"
  putStr $ ">>>>>>> Preconditioner <<<<<<<\n\n"
  putStr $ printedData2 ++ "\n\n"
  putStr " << Total number of terms >> "
  putStr $ (show $ (foldl (+) 0 $ fmap (\x -> length $ snd $ snd x) sigmaElements) + (foldl (+) 0 $ fmap (\x -> length $ snd x) precondElements)) ++ "\n\n"
  --putStr $ (show $ (foldl (+) 0 $ fmap (\x -> length $ snd $ snd x) sigmaElements)) ++ "\n\n"  
   
