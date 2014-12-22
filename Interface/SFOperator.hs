-- ///////////////////////////////////////////////////////////////////////////
--
--  Copyright (c) 2013 - 2014 Masaaki Saitow. All rights reserved.
--  The source code in this file is a part of
--
--  Fafnir :: A Symbolic Manipulation Module for Many-Fermionic Operators.
--
--  Any modifications or use of this software should cite the following
--  paper on which the algorithm implemented in this package is based and
--  femto itself:
--
--  [1] M. Saitow, Y.Kurashige and T. Yanai, J. Chem. Phys. 139, 044118 (2013)
--  [2] M. Saitow, Femto :: An Integrated Toolset for the Automated Tensor 
--  Generation
--
--  Currently, the written license for obtaining, storing and using this code
--  has never been prepared. Note that the author is the only holder of the 
--  rights to distribute both precompiled form and the source code of this 
--  software, i.e. no others can re-distribute a copy of this software without 
--  any permission of the author.
--
--  For questions concerning the software or copyright write to the author:
--
--                             Masaaki Saitow
--                         (msaitow514@gmail.com)
--
-- ///////////////////////////////////////////////////////////////////////////

--------------------------------------
-- Handy interface -------------------
--------------------------------------
module Interface.SFOperator(

printList,
strList,
constructSigmaCommBase,
constructSigmaComm,
constructSigmaBase,
constructSigma,
constructHamiltonian,
constructHamiltonianComm,
scaleTerm,
screenTerms,
processTermBaseSF,
processTermsE,
processTermsCommE,
transposeIndices,

) where

import Data.List
import Data.Maybe
import Interface.SFIndex
--import Fafnir
import BraHKet

----------------------------------------------------------
-- Dummy indices used here
----------------------------------------------------------
_p = sfIndex "_p" Generic Dummy      
_q = sfIndex "_q" Generic Dummy           
_r = sfIndex "_r" Generic Dummy      
_s = sfIndex "_s" Generic Dummy
----------------------------------------------------------

-- Print each member into the screen
printList :: (Show a) => (Ord a) => [a] -> IO ()
printList myList = do
  putStr $ strList myList

-- Convert each member in the given list into a formatted string
strList :: (Show a) => (Ord a) => [a] -> String
strList myList = if length myList == 0 then "" else foldl (++) "" $ fmap (\x -> x ++ "\n") $ zipWith (++) myNumbers $ fmap show myList
  where myNumbers = fmap (\x -> "  [" ++ (show x) ++ "] " ) [0..(length myList)-1]

------------------------------------------------------------------------------------------------------------------------
-- It derives the explicit formulas of the IC-MRCI sigma-equation (Hamiltonian-Vector product)
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

-- Returns the commutator-based Hamiltonian-Vector product
constructSigmaCommBase :: QIndices -> QIndices -> QIndices -> QTerms
constructSigmaCommBase eLinds eRinds eTinds = fockTerms
  where

    -- Definitions of tensors
    t   = sfTensor "@T2" eTinds [[0,1,2,3],[1,0,3,2]] Classical 
    eL  = baseSFGen eLinds
    eR  = baseSFGen eRinds 
    h   = baseOne   [_p, _q]         -- h^{p}_{q}
    v   = baseERI   [_p, _q, _r, _s] -- v^{pq}_{rs}
    eH1 = baseSFGen [_p, _q]         -- E^{p}_{r}
    eH2 = baseSFGen [_p, _q, _r, _s] -- E^{pq}_{rs}

    ----------------------------------------------------------------------------------------------

         --           
         --  Define: E^{left} h^{p}_{q} E^{p}_{q} T^{right} E_{right}
         --

         --                       1
         --                    + --- E^{left} V^{pq}_{rs} E^{pq}_{rs} T^{right} E_{right}
         --                       2

  ------------------------------------------------------------------------------------------------          
    terms      = [baseTerm 1.0 [] [eL, h, eH1, t, eR], baseTerm 0.5 [] [eL, v, eH2, t, eR]]

    ordered    = combineTerms $ screenTerms $ takeVEV $ (concat $ fmap normalOrderCommE terms) ++ (normalOrderE $ baseTerm 1.0 ["E0"] [eL, t, eR])
    decomposed = combineTerms $ screenTerms $ concat $ fmap contractCoreSF $ fmap masquerade $ contractVirtSF $ concat $ fmap (generateInteractions [Core, Active, Virtual]) ordered
    fockTerms  = constructCoreFock $ decomposed

constructSigmaComm :: QIndices -> QIndices -> QTerms
constructSigmaComm eLinds eRinds = constructSigmaCommBase eLinds eRinds eT
  where
    eT = [(eRinds !! 2), (eRinds !! 3), (eRinds !! 0), (eRinds !! 1)]

-- Returns the Hamiltonian-Vector product
constructSigmaBase :: QIndices -> QIndices -> QIndices -> QTerms
constructSigmaBase eLinds eRinds eTinds = fockTerms
  where

    -- Definitions of tensors
    t   = sfTensor "@T2" eTinds [[0,1,2,3],[1,0,3,2]] Classical 
    eL  = baseSFGen eLinds
    eR  = baseSFGen eRinds 
    h   = baseOne   [_p, _q]         -- h^{p}_{q}
    v   = baseERI   [_p, _q, _r, _s] -- v^{pq}_{rs}
    eH1 = baseSFGen [_p, _q]         -- E^{p}_{r}
    eH2 = baseSFGen [_p, _q, _r, _s] -- E^{pq}_{rs}

    ----------------------------------------------------------------------------------------------

         --           
         --  Define: E^{left} h^{p}_{q} E^{p}_{q} T^{right} E_{right}
         --

         --                       1
         --                    + --- E^{left} V^{pq}_{rs} E^{pq}_{rs} T^{right} E_{right}
         --                       2

  ------------------------------------------------------------------------------------------------          
    terms      = [baseTerm 1.0 [] [eL, h, eH1, t, eR], baseTerm 0.5 [] [eL, v, eH2, t, eR]]

    ordered    = combineTerms $ screenTerms $ takeVEV $ concat $ fmap normalOrderE terms
    decomposed = combineTerms $ screenTerms $ concat $ fmap contractCoreSF $ fmap masquerade $ contractVirtSF $ concat $ fmap (generateInteractions [Core, Active, Virtual]) ordered
    fockTerms  = constructCoreFock $ decomposed

constructSigma :: QIndices -> QIndices -> QTerms
constructSigma eLinds eRinds = constructSigmaBase eLinds eRinds eT
  where
    eT = [(eRinds !! 2), (eRinds !! 3), (eRinds !! 0), (eRinds !! 1)]

-- Returns the explicit Hamiltonian element
constructHamiltonian :: QIndices -> QIndices -> QTerms
constructHamiltonian eLinds eRinds = fockTerms
  where

    -- Definitions of tensors
    eL  = baseSFGen eLinds
    eR  = baseSFGen eRinds 
    h   = baseOne   [_p, _q]         -- h^{p}_{q}
    v   = baseERI   [_p, _q, _r, _s] -- v^{pq}_{rs}
    eH1 = baseSFGen [_p, _q]         -- E^{p}_{r}
    eH2 = baseSFGen [_p, _q, _r, _s] -- E^{pq}_{rs}

    ----------------------------------------------------------------------------------------------

         --           
         --  Define: E^{left} h^{p}_{q} E^{p}_{q} T^{right} E_{right}
         --

         --                       1
         --                    + --- E^{left} V^{pq}_{rs} E^{pq}_{rs} T^{right} E_{right}
         --                       2

  ------------------------------------------------------------------------------------------------          
    terms      = [baseTerm 1.0 [] [eL, h, eH1, eR], baseTerm 0.5 [] [eL, v, eH2, eR]]

    ordered    = combineTerms $ screenTerms $ takeVEV $ concat $ fmap normalOrderE terms
    decomposed = combineTerms $ screenTerms $ concat $ fmap contractCoreSF $ fmap masquerade $ contractVirtSF $ concat $ fmap (generateInteractions [Core, Active, Virtual]) ordered
    fockTerms  = constructCoreFock $ decomposed

-- Returns the explicit Hamiltonian element
constructHamiltonianComm :: QIndices -> QIndices -> QTerms
constructHamiltonianComm eLinds eRinds = fockTerms
  where

    -- Definitions of tensors
    eL  = baseSFGen eLinds
    eR  = baseSFGen eRinds 
    h   = baseOne   [_p, _q]         -- h^{p}_{q}
    v   = baseERI   [_p, _q, _r, _s] -- v^{pq}_{rs}
    eH1 = baseSFGen [_p, _q]         -- E^{p}_{r}
    eH2 = baseSFGen [_p, _q, _r, _s] -- E^{pq}_{rs}

    ----------------------------------------------------------------------------------------------

         --           
         --  Define: E^{left} h^{p}_{q} E^{p}_{q} T^{right} E_{right}
         --

         --                       1
         --                    + --- E^{left} V^{pq}_{rs} E^{pq}_{rs} T^{right} E_{right}
         --                       2

  ------------------------------------------------------------------------------------------------          
    terms      = [baseTerm 1.0 [] [eL, h, eH1, eR], baseTerm 0.5 [] [eL, v, eH2, eR] ]
    
    ordered    = combineTerms $ screenTerms $ takeVEV $ (concat $ fmap normalOrderCommE terms) ++ (normalOrderE $ baseTerm 1.0 ["E0"] [eL, eR]) 
    decomposed = combineTerms $ screenTerms $ concat $ fmap contractCoreSF $ fmap masquerade $ contractVirtSF $ concat $ fmap (generateInteractions [Core, Active, Virtual]) ordered
    fockTerms  = constructCoreFock $ decomposed


-- Returns the terms divided by a factor
scaleTerm :: Double -> QTerm -> QTerm
scaleTerm num inTerm = baseTerm ((tNum inTerm) * num) (tCoeff inTerm) (tTensor inTerm) 

-- Screen the NULL terms
screenTerms :: QTerms -> QTerms
screenTerms inTerm = fmap (fromJust) $ filter (\x -> x /= Nothing) $ fmap (killKDeltas) inTerm

--------------------------------------------------------------------------------------------------------
-- Transform the given many-body ansatz into a stream of the programmable tensor contraction terms
--------------------------------------------------------------------------------------------------------
-- [Input] 
--   givenNormalOrderSF :: A function for the normal ordering (QTerm->QTerms)
--   targetSpace        :: Group of the interaction to which the generic indices are decomposed
--   inTerms            :: The many-body ansatz
--------------------------------------------------------------------------------------------------------
processTermBaseSF :: (QTerm -> QTerms) -> [QSpace] -> QTerms -> QTerms
processTermBaseSF givenNormalOrderSF targetSpace inTerms = combineTerms $ screenTerms $ contractVirtSF $ concat $ fmap contractCoreSF $ concat $ fmap (generateInteractions targetSpace) $ contractVirtSF $ takeVEV $ concat $ fmap givenNormalOrderSF inTerms

-- Transform the normal ordered terms with the spin-free generators into the programable form
processTermsE :: [QSpace] -> QTerms -> QTerms
processTermsE targetSpace inTerms = processTermBaseSF normalOrderE targetSpace inTerms

-- Transform the normal ordered terms with the spin-free generators into the programable form through the normal ordering of the multiple-commutator, i.e. <|[A,[B,[C,[..,X]..]]]|>
processTermsCommE :: [QSpace] -> QTerms -> QTerms
processTermsCommE targetSpace inTerms = processTermBaseSF normalOrderCommE targetSpace inTerms

-- Returns a transposed list of indices, e.g. (a^i a^j a_s a_r)^t -> a^r a^s a_j a_i
transposeIndices :: QIndices -> QIndices
transposeIndices myInds = if odd $ length myInds then error $ "transposeIndices: The length of given indices is odd. >> " ++ (show myInds) ++ " << \n"
                          else (reverse $ take order $ reverse myInds) ++ take order myInds
  where
    order = floor $ (fromIntegral $ length myInds) / (fromIntegral 2)
    
