-- Copyright (C) 2013-2014 by Masaaki Saitow (msaitow514@gmail.com)
--
-- This program is free software; you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2 of the License, or
-- (at your option) any later version
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program; if not, write to the Free Software
-- Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

-----------------------------------------------------
-- Wick expansion of the many-fermionic ansatz
-----------------------------------------------------
----module Fafnir.Wick
module BraHKet.Wick
(

  normalOrderA,         -- NO for Cre/Des operators
  normalOrderG,         -- NO for spin-dependent generator
  normalOrderE,         -- NO for spin-free generator
  normalOrderCommE,     -- NO for the multiple commutators of the spin-free generator
  contractGen,          -- Contract two spin-free generators
  takeVEV,              -- Take vacuum expectation value
  contractCoreSF,       -- Contarct core operator in the spin-free density matrix
  contractVirtSF        -- Screen terms with the  virtual operator in the density
  
) where

import qualified BraHKet.Utils as Utils 
import qualified Data.Map as Map
import qualified Data.Maybe as Maybe
import qualified Data.List as List
import BraHKet.Core

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Normal ordering function for creation and annihilation operators
normalOrderA :: QTerm -> QTerms
normalOrderA term =
  if length cres == length des then zipWith (makeTerm) signs kDeltas
  else error "normalOrderA: Numbers of creation and annihilation operators should be equal for the current implementation"
  where
    tensors   = tTensor term
    operators = [x | x <- tensors, tLabel x == creName_ || tLabel x == desName_]
    others    = [x | x <- tensors, tLabel x /= creName_ && tLabel x /= desName_]
    cres      = [x | x <- operators, tLabel x == creName_]
    des       = [x | x <- operators, tLabel x == desName_]
    contPairs = [(x, y) | x <- cres, y <- des]

    -- Actually survived contraction pairs
    survivedPairs = fmap (fmap reorderOps) $ filter (makeAllOps) $ Utils.binCombi (length cres) contPairs
    flatPairs     = fmap Utils.makeFlat survivedPairs
    signs         = fmap (Utils.permuteSign operators) flatPairs

    kDeltas = fmap (fmap makekDeltas) survivedPairs

    -- Returns proper contraction pairs (assuming all index pairs are composed of one creation and one annihilation operators)
    reorderOps :: (QTensor, QTensor) -> (QTensor, QTensor)
    reorderOps (op1, op2)
      | iSpace ((tIndices op1) !! 0) == Virtual && iSpace ((tIndices op2) !! 0) == Virtual = (op2, op1)
      | iSpace ((tIndices op1) !! 0) == Generic && iSpace ((tIndices op2) !! 0) == Virtual = (op2, op1)
      | iSpace ((tIndices op1) !! 0) == Virtual && iSpace ((tIndices op2) !! 0) == Generic = (op2, op1)
      | iSpace ((tIndices op1) !! 0) == Active  || iSpace ((tIndices op2) !! 0) == Active  = error "normalOrderA: Can't handle index for the active MOs"
      | otherwise = (op1, op2)
                                                                                         
    makeTerm :: Int -> QTensors -> QTerm
    makeTerm sign kDs = QTerm ((fromIntegral sign)*(tNum term)) (tCoeff term) (others ++ kDs)

    makekDeltas :: (QTensor, QTensor) -> QTensor
    makekDeltas (ind1, ind2) = baseKD $ (tIndices ind1) ++ (tIndices ind2)
    
    makeAllOps :: [(QTensor, QTensor)] -> Bool
    makeAllOps contras = (makeCre contras) && (makeDes contras)
      where
        makeCre kore = (List.sort $ fmap (fst) kore) == List.sort cres
        makeDes are  = (List.sort $ fmap (snd) are)  == List.sort des


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Normal ordering function for spin-dependent excitation operators
normalOrderG :: QTerm -> QTerms
normalOrderG term = zipWith (makeTerm) signs kDeltas
  where
    tensors   = tTensor term
    operators = [x | x <- tensors, take (length sdGenName_) (tLabel x) == sdGenName_]
    others    = [x | x <- tensors, take (length sdGenName_) (tLabel x) /= sdGenName_]

    -- Extract groups of creation and annihilation operators of generator
    creGroup = fmap (transformCre) operators
    desGroup = fmap (transformDes) operators
    opGroup  = fmap (transformBoth) operators

    cres = concat creGroup
    des  = concat desGroup
    ops  = concat opGroup
    conPairs = [(x, y) | x <- cres, y <- des]

    -- List of fully-contracted operator pairs 
    survivedPairs = fmap (fmap reorderOps) $ filter (makeAllOps) $ Utils.binCombi (length cres) conPairs
    deadPairs = concat $ fmap (makeDeadPairs) opGroup
    allPairs = filter (isAlive) survivedPairs
    
    flatPairs = fmap Utils.makeFlat allPairs
    signs     = fmap (Utils.permuteSign ops) flatPairs

    kDeltas = fmap (fmap makekDeltas) allPairs

    -- Returns proper contraction pairs (assuming all index pairs are composed of one creation and one annihilation operators)
    reorderOps :: (QIndex, QIndex) -> (QIndex, QIndex)
    reorderOps (op1, op2)
      | iSpace op1 == Virtual && iSpace op2 == Virtual = (op2, op1)
      | iSpace op1 == Generic && iSpace op2 == Virtual = (op2, op1)
      | iSpace op1 == Virtual && iSpace op2 == Generic = (op2, op1)
      | iSpace op1 == Active  || iSpace op2 == Active  = error "normalOrderG: Can't handle index for the active MOs"
      | otherwise = (op1, op2)

    -- If contras contains either f the dear pairs, contras vanish
    isAlive :: [(QIndex, QIndex)] -> Bool
    isAlive contras = not $ allNothing deadPairs
      where
        allNothing [x]    = x `elem` contras
        allNothing (x:xs) = x `elem` contras || allNothing xs

    -- Make list of intra-group contractions
    makeDeadPairs :: QIndices -> [(QIndex, QIndex)]
    makeDeadPairs opList =
      let
        order = floor $ (fromIntegral . length $ opList)/2
        creOps = take order opList
        desOps = take order $ reverse opList
      in [(x,y) | x <- creOps, y <- desOps]
    
    transformCre :: QTensor -> QIndices
    transformCre sdGen = 
      let
        order = floor $ (fromIntegral $ length $ tIndices sdGen) / (fromIntegral 2)
      in take order (tIndices sdGen)

    transformDes :: QTensor -> QIndices
    transformDes sdGen = 
      let
        order = floor $ (fromIntegral $ length $ tIndices sdGen) / (fromIntegral 2)
      in take order $ reverse (tIndices sdGen)

    transformBoth :: QTensor -> QIndices
    transformBoth sdGen = (transformCre sdGen) ++ (transformDes sdGen)

    makeAllOps :: [(QIndex, QIndex)] -> Bool
    makeAllOps contras = (makeCre contras) && (makeDes contras)
      where
        makeCre kore = (List.sort $ fmap (fst) kore) == List.sort cres
        makeDes are  = (List.sort $ fmap (snd) are)  == List.sort des

    makeTerm :: Int -> QTensors -> QTerm
    makeTerm sign kDs = QTerm ((fromIntegral sign)*(tNum term)) (tCoeff term) (others ++ kDs)

    makekDeltas :: (QIndex, QIndex) -> QTensor
    makekDeltas (ind1, ind2) = baseKD $ ind1 : ind2 : []


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Normal ordering function for spin-free excitaiton operators
normalOrderE :: QTerm -> QTerms
normalOrderE term = if length otherOps /= 0
                    then error "normalOrderE: Normal ordering among other types of excitation operators is not yet implemented."
                    else iterateNormalOrder restEops [initTerm]
  where
    tensors   = tTensor term
    operators = [x | x <- tensors, take (length sfGenName_) (tLabel x) == sfGenName_]
    others    = [x | x <- tensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
    otherOps  = [x | x <- others, tComm x == Operator]        
    initTerm  = makeTerm' $ others ++ [operators !! (length operators-1)]
    restEops  = take (length operators-1) operators
    numSteps  = (length operators) - 1

    iterateNormalOrder :: QTensors -> QTerms -> QTerms
    iterateNormalOrder restEs currentTerms
      | length restEs == 0 = currentTerms
      | otherwise          = iterateNormalOrder (init restEs) orderedTerms
      where
        newTensors   = zipWith (:) ( replicate (length currentTerms) (last restEs) ) $ fmap (tTensor) currentTerms
        newTerms     = fmap (makeTerm') newTensors -- :: QTerms
        orderedTerms = concat $ fmap (contractGen) newTerms
        
    -- Returns terms with correct prefactors and tensors
    makeTerm' :: QTensors -> QTerm
    makeTerm' tenten = baseTerm (tNum term) (tCoeff term) (tenten)


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Normal ordering function for the multiple commutators of the spin-free excitaiton operators, [E1, [E2, [E3, .. En]] .. ]
normalOrderCommE :: QTerm -> QTerms
normalOrderCommE term = if length otherOps /= 0
                    then error "normalOrderCommE: Normal ordering among other types of excitation operators is not yet implemented."
                    else iterateNormalOrder restEops [initTerm]
  where
    tensors   = tTensor term
    operators = [x | x <- tensors, take (length sfGenName_) (tLabel x) == sfGenName_]
    others    = [x | x <- tensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
    otherOps  = [x | x <- others, tComm x == Operator]    
    initTerm  = makeTerm' (tNum term) $ others ++ [operators !! (length operators-1)]
    restEops  = take (length operators-1) operators
    numSteps  = (length operators) - 1

    iterateNormalOrder :: QTensors -> QTerms -> QTerms
    iterateNormalOrder restEs currentTerms
      | length restEs == 0 = currentTerms
      | otherwise          = iterateNormalOrder (init restEs) orderedTerms
      where
        newTensors   = zipWith (:) ( replicate (length currentTerms) (last restEs) ) $ fmap (tTensor) currentTerms
        factors      = fmap (tNum) currentTerms
        newTerms     = zipWith (makeTerm') factors newTensors -- :: QTerms
        orderedTerms = concat $ fmap (contractComm) newTerms
        
    -- Returns terms with correct prefactors and tensors
    makeTerm' :: Double -> QTensors -> QTerm
    makeTerm' myFactor tenten = baseTerm myFactor (tCoeff term) (tenten)

    -- Returns normal ordered results of [E1, E2]
    contractComm :: QTerm -> QTerms
    contractComm thisTerm
      | length es > 2  = error "contractComm: Length of operators given should be  shorter than or equal to 2."
      | length es == 2 = withoutMaximum_e1e2 ++ withoutMaximum_e2e1
      | otherwise      = [thisTerm]
      where
        thisTensors = tTensor thisTerm 
        commutables = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
        es = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) == sfGenName_]
        negativeTerm = baseTerm (negate $ tNum thisTerm) (tCoeff thisTerm) (commutables ++ [es !! 1, es !! 0])
        withoutMaximum_e1e2 = filter (\x -> not $ isMaxE x) $ contractGen thisTerm
        withoutMaximum_e2e1 = filter (\x -> not $ isMaxE x) $ contractGen negativeTerm

        -- Returns true if the term has the generator of maximum rank that can be formed from E1 and E2
        isMaxE :: QTerm -> Bool
        isMaxE kore
          | length myEs /= 1 = False
          | otherwise        = length (tIndices $ myEs !! 0) == length (tIndices $ es !! 0) + length (tIndices $ es !! 1)
          where
            myTensors = tTensor kore
            myEs = [x | x <- myTensors, take (length sfGenName_) (tLabel x) == sfGenName_]


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Contract given two spin-free excitation operators
contractGen :: QTerm -> QTerms
contractGen thisTerm
  | allSF /= True               = error "contractGen: The spin-dependent indices detected."
  | (length otherOps) /= 0      = error "contractGen: Normal ordering with other type of operator is not yet implemented."
  | (length incommutables)  > 2 = error "contractGen: Length of operators given should be shorter than or qeual to 2."
  | (length incommutables) == 2 = concat . fmap (fmap makeTerm) $ fmap (makeContractions) [0..(min order1 order2)]
  | otherwise                                             = [thisTerm]
  where
    allSF = foldl (&&) True $ fmap (\x -> iSpin x == SpinFree) $ concat $ fmap tIndices $ tTensor thisTerm                
    thisTensors = tTensor thisTerm
    commutables   = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
    incommutables = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) == sfGenName_]
    otherOps      = [x | x <- commutables, tComm x == Operator]
    e1 = incommutables !! (length incommutables-2)
    e2 = incommutables !! (length incommutables-1)
    order1     = floor $ (fromIntegral $ length $ tIndices e1) / (fromIntegral 2)
    order2     = floor $ (fromIntegral $ length $ tIndices e2) / (fromIntegral 2)
    upperInds1 = take (order1) (tIndices e1)
    upperInds2 = take (order2) (tIndices e2)
    lowerInds1 = reverse $ take (order1) $ reverse (tIndices e1)
    lowerInds2 = reverse $ take (order2) $ reverse (tIndices e2) 

    -- Lower indices of e1 and upper indices of e2 are used to make contraction
    contPairs = [(x, y) | x <- lowerInds1, y <- upperInds2]
    e1PairInds = Map.fromList $ zip lowerInds1 upperInds1  -- -> Lower indices are keys -> Upper indices
    e2PairInds = Map.fromList $ zip upperInds2 lowerInds2  -- -> Upper indices are keys -> Lower indices

    e1Pairs = zip upperInds1 lowerInds1  
    e2Pairs = zip upperInds2 lowerInds2  

    e1Maps = Map.fromList $ zip upperInds1 e1Pairs  -- -> Upper indices are keys -> (Upper, Lower)
    e2Maps = Map.fromList $ zip lowerInds2 e2Pairs  -- -> Lower indices are keys -> (Upper, Lower)

    -- Returns terms with correct prefactors and tensors
    makeTerm :: QTensors -> QTerm
    makeTerm tenten = baseTerm (tNum thisTerm) (tCoeff thisTerm) (commutables ++ tenten)

    -- Returns contracted operators [input; contraction order]
    makeContractions :: Int -> [QTensors]
    makeContractions myOrder = filter (isAppropriate) $ zipWith (:) (fmap (baseSFGen) $ fmap (Utils.makeFlat2) totalInds) kDeltas
      where
        myContractions = Utils.binCombi myOrder contPairs
        kDeltas      = fmap (fmap makeKD) myContractions       -- Kronecker's deltas [[kDs]]
        newPairs     = fmap (fmap makeNewPairs) myContractions -- Newly formed pairs [[(uind_e1, lind_e2)]]
        -----------------------------
        deadE1Pairs  = fmap (fmap Maybe.fromJust) $ fmap (filter (\x -> x /= Nothing)) $ fmap (fmap findDeadE1) $ fmap (fmap fst) newPairs
          where
            findDeadE1 :: QIndex -> Maybe (QIndex, QIndex)
            findDeadE1 konoIndex = if Map.lookup konoIndex e1Maps /= Nothing then Map.lookup konoIndex e1Maps else Nothing
        
        deadE2Pairs  = fmap (fmap Maybe.fromJust) $ fmap (filter (\x -> x /= Nothing)) $ fmap (fmap findDeadE2) $ fmap (fmap snd) newPairs
          where
            findDeadE2 :: QIndex -> Maybe (QIndex, QIndex)
            findDeadE2 konoIndex = if Map.lookup konoIndex e2Maps /= Nothing then Map.lookup konoIndex e2Maps else Nothing
        -----------------------------
        aliveE1Pairs = fmap (isAliveE1) deadE1Pairs
        aliveE2Pairs = fmap (isAliveE2) deadE2Pairs            
        totalInds    = zipWith (++) aliveE1Pairs $ zipWith (++) aliveE2Pairs newPairs

        makeKD :: (QIndex, QIndex) -> QTensor
        makeKD (lind1, uind2) = baseKD [lind1, uind2]

        makeNewPairs :: (QIndex, QIndex) -> (QIndex, QIndex)
        makeNewPairs (lind1, uind2) = (Maybe.fromJust $ Map.lookup lind1 e1PairInds, Maybe.fromJust $ Map.lookup uind2 e2PairInds) -- -> Return (Upper, Lower) pairs

        isAliveE1 :: [(QIndex, QIndex)] -> [(QIndex, QIndex)]
        isAliveE1 koreraE1 = [x| x <- e1Pairs, not $ x `elem` koreraE1]

        isAliveE2 :: [(QIndex, QIndex)] -> [(QIndex, QIndex)]
        isAliveE2 koreraE2 = [x| x <- e2Pairs, not $ x `elem` koreraE2]

        isAppropriate :: QTensors -> Bool
        isAppropriate tenten = List.sort (myIndices) == List.sort (Utils.makeFlat $ e1Pairs ++ e2Pairs)
          where
            myIndices = concat $ fmap (tIndices) tenten
        

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Construct spin-free density matrix from the normal ordered spin-free generators by taking the vacuum expectaiton values
takeVEV :: QTerms -> QTerms
takeVEV terms = fmap Maybe.fromJust $ filter (\x -> x /= Nothing) $ fmap (uniVEV) terms
  where
    -- VEV for one term
    uniVEV :: QTerm -> Maybe QTerm
    uniVEV term
      | length sfOps > 1 = Just term
      | otherwise        = if Virtual `elem` allSpaces then Nothing else Just newTerm
      where
        tensors = tTensor term
        sfOps = [x | x <- tensors, take (length sfGenName_) (tLabel x) == sfGenName_]
        allSpaces = fmap (iSpace) $ concat $ fmap (tIndices) sfOps
        newTerm = baseTerm (tNum term) (tCoeff term) (fmap (makeNewTensor) tensors)　

        makeNewTensor :: QTensor -> QTensor
        makeNewTensor konoTensor
          | take (length sfGenName_) (tLabel konoTensor) == sfGenName_ = baseSFRDM (tIndices konoTensor)
          | otherwise                                                  = konoTensor


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Contract core operator in the spin-free density matrix
contractCoreSF :: QTerm -> QTerms
contractCoreSF thisTerm
  | length rdms /= 1                       = [thisTerm]
  | allSF /= True                          = error   "contractCoreSF: The spin-dependent indices detected inthe input term"
  | length genericIndices /= 0             = error $ "contractCoreSF: There are at least one Generic indices which cannot be treated in this function. \n" ++ (show thisTerm) ++ "\n"
  | length ucorePairs /= length lcorePairs = []
  | otherwise                              = zipWith3 baseTerm newFactors newCoeffs newTensors
  where
    allSF = foldl (&&) True $ fmap (\x -> iSpin x == SpinFree) $ concat $ fmap tIndices $ tTensor thisTerm            
    factor  = tNum thisTerm
    tensors = tTensor thisTerm
    rdms    = [x | x <- tensors, take (length sfRDMName_) (tLabel x) == sfRDMName_]
    others  = [x | x <- tensors, take (length sfRDMName_) (tLabel x) /= sfRDMName_]
    theRDM  = rdms !! 0
    order   = floor $ (fromIntegral $ length $ tIndices theRDM) / (fromIntegral 2)
    genericIndices = [x | x <- tIndices theRDM, (iSpace x) == Generic]

    -- Index-Number pairs [(QIndex, Int)]
    indPairs   = zip (tIndices theRDM) [0..((length $ tIndices theRDM)-1)]
    upperPairs = take (order) indPairs
    lowerPairs = fmap (\(x,y)->(x,y-order)) $ reverse $ take (order)$ reverse indPairs 

    -- Number-Index pairs [(Int, QIndex)]
    numPairs  = fmap (\(x,y)->(y,x)) indPairs
    unumPairs = fmap (\(x,y)->(y,x)) upperPairs
    lnumPairs = fmap (\(x,y)->(y,x)) lowerPairs

    ------------------
    unumMap   = Map.fromList unumPairs -- Position :: Int -> QIndex
    lnumMap   = Map.fromList lnumPairs -- Position :: Int -> QIndex
    ------------------
    
    ucorePairs = [x | x <- upperPairs, (iSpace $ fst x) == Core]
    lcorePairs = [x | x <- lowerPairs, (iSpace $ fst x) == Core]
    
    ucPos = fmap (snd) ucorePairs -- Position of the upper core index
    lcPos = fmap (snd) lcorePairs -- Position of the lower core index

    -- Form the core contraction pairs to constract the Kronecker's delta
    lcTemp = List.permutations lcPos
    ucTemp = take (length lcTemp) $ repeat ucPos
    corePos = zip ucTemp lcTemp

    ------------------------------------
    -- Body of this function -----------
    ------------------------------------
    kDeltas   = fmap formKDeltas corePos
    newGenPos = fmap formNewGen corePos
    newGens   = fmap formSFGen newGenPos
    --totalPos  = zipWith (combineListP) corePos newGenPos
    totalPos  = zipWith (\(x1,y1) (x2,y2) -> (x1++x2, y1++y2)) corePos newGenPos    
    signs     = fmap determineSign totalPos
    factors   = fmap determineFactor corePos

    newFactors = fmap (factor*) $ zipWith (*) signs factors
    newTensors = fmap (++others) $ zipWith (combineTensorList) kDeltas newGens 
    newCoeffs  = take (length newTensors) $ repeat (tCoeff thisTerm)
    -------------------------------------
    -------------------------------------
    
    -- Returns the index list for forming the new spin-free excitation operator
    formNewGen :: ([Int], [Int]) -> ([Int], [Int])
    formNewGen myCoreInds = if List.sort lncInds == List.sort newlncInds then (uncInds, newlncInds) else error "fromNewGen: Algorithmic error occured."
      where
        --ucMap = Map.fromList $ zip (fst myCoreInds) (snd myCoreInds) -- -> Upper index is a key -> lower index
        lcMap = Map.fromList $ zip (snd myCoreInds) (fst myCoreInds) -- -> Lower index is a key -> upper index

        uncInds    = filter (\x -> not $ x `elem` (fst myCoreInds)) [0..(order-1)]
        lncInds    = filter (\x -> not $ x `elem` (snd myCoreInds)) [0..(order-1)]
        newlncInds = fmap searchIndex uncInds
        
        -- Seeks for the lower index of the generator
        searchIndex :: Int -> Int
        searchIndex thisPos
          | nextUc == Nothing                                     = thisPos
          | not ((Maybe.fromJust nextUc) `elem` (snd myCoreInds)) = Maybe.fromJust nextUc
          | otherwise                                             = searchIndex $ Maybe.fromJust nextUc
          where nextUc = Map.lookup thisPos lcMap

    -- Determine the factors of each term
    determineFactor :: ([Int], [Int]) -> Double
    determineFactor myCoreInds = fromIntegral $ 2 ^ (length $ Utils.decompPerm $ Utils.eliminatePerm myCoreInds)
    
    -- Determine the signs of the term
    determineSign :: ([Int], [Int]) -> Double
    determineSign myAllInds = fromIntegral $ Utils.permuteSign (fst myAllInds) (snd myAllInds)

    -- Returns upper index
    returnUInds :: Int -> QIndex
    returnUInds myPos = if newPos == Nothing then error "returnUInds: Algorithmic error occured." else Maybe.fromJust newPos
      where newPos = Map.lookup myPos unumMap

    -- Returns lower index
    returnLInds :: Int -> QIndex
    returnLInds myPos = if newPos == Nothing then error "returnLInds: Algorithmic error occured." else Maybe.fromJust newPos
      where newPos = Map.lookup myPos lnumMap

    -- Form kdeltas
    formKDeltas :: ([Int], [Int]) -> QTensors
    formKDeltas myCoreInds = fmap (\(x,y) -> baseKD [x,y]) kDInds
      where kDInds = zip (fmap returnUInds $ fst myCoreInds) (fmap returnLInds $ snd myCoreInds)

    -- Form SFGen
    formSFGen :: ([Int], [Int]) -> Maybe QTensor
    formSFGen myGenInds
      | length (fst myGenInds) /= 0 = Just $ baseSFRDM $ (fmap returnUInds (fst myGenInds)) ++ (fmap returnLInds (snd myGenInds))
      | otherwise                   = Nothing

    -- Combine the Maybe SFGen and Kronecker's deltas to form newTensors :: QTensors
    combineTensorList :: QTensors -> Maybe QTensor -> QTensors
    combineTensorList kds gen
      | gen == Nothing = kds
      | otherwise      = (Maybe.fromJust gen) : kds


---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Screen terms with virtual operator in the spin-free reduced density matrix
contractVirtSF :: QTerms -> QTerms
contractVirtSF koreraTerms = filter (allClear) koreraTerms
  where
    allClear :: QTerm -> Bool
    allClear kore
      | length rdms   == 0 = True
      | length virOps == 0 = True
      | otherwise          = False
      where
        tensors = tTensor kore
        rdms    = [x | x <- tensors, take (length sfRDMName_) (tLabel x) == sfRDMName_]
        virOps  = filter (\x -> iSpace x == Virtual) $ List.nub $ concat $fmap (tIndices) rdms
