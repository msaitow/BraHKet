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

-----------------------------------------------
-- Formation of the various Fock matrix
-----------------------------------------------
---module Fafnir.Fock
module BraHKet.Fock
(
  makeEffectiveInts, -- Construct the effective integrals
  constructCoreFock  -- Construct core Fock matrix by combining the core integrals
  
) where

import qualified BraHKet.Utils as Utils 
import qualified Data.Map as Map
import qualified Data.Maybe as Maybe
import qualified Data.List as List
import BraHKet.Core

-----------------------------------------------------------------
-- @@ Definition of the effective integrals
-----------------------------------------------------------------
h1Eff_  = "@h1Int_"
h2Eff_  = "@h2Int_"
h3Eff_  = "@h3Int_"
h4Eff_  = "@h4Int_"
h5Eff_  = "@h5Int_"
fcEff_  = "@Fcore"    
ecEff_  = "Ecore"
casEff_ = "Ecas"    
-----------------------------------------------------------------

-------------------------------------------------------------------------------------------------------------
-- << Construct the effective integrals >> ------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------
-- I. Replace the core integrals by the effective integrals -------------------------------------------------
-------------------------------------------------------------------------------------------------------------
-- Make sure this algorithm is valid only for the CI, CC and CT ansatz. In general, it cannot be used for
-- the PT approaches since term of the PT equations may contain more than one integrals.
-------------------------------------------------------------------------------------------------------------
makeEffectiveInts :: QTerm -> QTerm
makeEffectiveInts givenTerm
  | hasH1     = replaceH1 givenTerm
  | hasH2     = replaceH2 givenTerm
  | hasCFock  = replaceCFock givenTerm
  | otherwise = givenTerm
  where
    allNames = fmap tLabel $ tTensor givenTerm
    hasH1    = h1Name_      `elem` allNames
    hasH2    = h2Name_      `elem` allNames
    hasCFock = casFockName_ `elem` allNames

    -- Replace one-body integral by the h1 effective int
    replaceH1 myTerm    = replaceOneBody h1Name_ h1Eff_ myTerm

    -- Replace the trace of the Fock matrix by the effective int
    replaceCFock myTerm = moveCasEffToCoeff $ replaceOneBody casFockName_ casEff_ myTerm
      where
        moveCasEffToCoeff :: QTerm -> QTerm
        moveCasEffToCoeff thisTerm
          | hasCasFock = baseTerm (tNum thisTerm) ((take (length focks) $ repeat casEff_) ++ tCoeff thisTerm) others
          | otherwise  = thisTerm
          where
            hasCasFock      = casEff_ `elem` (fmap tLabel $ tTensor thisTerm)
            (focks, others) = List.partition (\x -> tLabel x == casEff_) $ tTensor thisTerm
        
    -- Replace one-body integral by the h1 effective int
    replaceOneBody :: String -> String -> QTerm -> QTerm
    replaceOneBody myName newName termh1 = if (inds !! 0) == (inds !! 1) && iSpace (inds !! 0) == Core && isDummy (inds !! 0) == Dummy then repTerm else termh1 
      where
        (h1Int, theRest) = List.partition (\x -> tLabel x == myName) $ tTensor termh1
        inds    = tIndices (h1Int !! 0)
        repTerm = baseTerm (tNum termh1) (tCoeff termh1) ((QTensor newName [] [] Classical) : theRest) 

    -- Replace two-body integral by either h2, h3, h4, or h5 effective int
    replaceH2 :: QTerm -> QTerm
    replaceH2 termh2
      | (caseType2 allInds && closedInds (h2Int !! 0) == 2) = baseTerm (tNum termh2) (tCoeff termh2) ((baseTensor h2Eff_ [] [] Classical) : theRest)
      | (caseType3 allInds && closedInds (h2Int !! 0) == 2) = baseTerm (tNum termh2) (tCoeff termh2) ((baseTensor h3Eff_ [] [] Classical) : theRest)
      | (caseType4 allInds && closedInds (h2Int !! 0) == 1) = baseTerm (tNum termh2) (tCoeff termh2) ((baseTensor h4Eff_ restInds [[0,1], [1,0]] Classical) : theRest)
      | (caseType5 allInds && closedInds (h2Int !! 0) == 1) = baseTerm (tNum termh2) (tCoeff termh2) ((baseTensor h5Eff_ restInds [[0,1], [1,0]] Classical) : theRest)
      | otherwise = termh2
      where
        (h2Int, theRest) = List.partition (\x -> tLabel x == h2Name_) $ tTensor termh2            
        allInds  = fmap tIndices $ tgetConfs (h2Int !! 0)
        restInds = openIndices $ h2Int !! 0
            
        ---------------------------------------------------------------------------------------------------------
        -- Now, analyze what kind of index dependence h2Int has to determine the type of the effective integrals
        ---------------------------------------------------------------------------------------------------------
        
        -------------------------------------------------------------
        -- Count the number of the closed indices
        -------------------------------------------------------------            
        closedInds :: QTensor -> Integer
        closedInds tempT = sum $ filter (\x -> x == 1) $ fmap sum $ fmap (haveInd tensorh2) myInds
          where
            tensorh2 = tTensor termh2
            myInds   = List.nub $ filter (\x -> isDummy x == Dummy) $ tIndices tempT

            -- If the given tensors have the index, returns 1, otherwise returns 0 
            haveInd :: QTensors -> QIndex -> [Integer]
            haveInd tens ind = fmap (hasInd ind) tens
              where
                hasInd :: QIndex -> QTensor -> Integer
                hasInd ind2 ten  = if ind2 `elem` (tIndices ten) then 1 else 0

        -------------------------------------------------------------
        -- Returns the open indices
        -------------------------------------------------------------            
        openIndices :: QTensor -> QIndices
        openIndices tempT = noDinds ++ (fmap snd $ filter (\x -> fst x /= 1) $ zip (fmap sum $ fmap (haveInd tensorh2) myInds) myInds)
          where
            tensorh2 = tTensor termh2
            myInds   = filter (\x -> isDummy x == Dummy)    $ tIndices tempT
            noDinds  = filter (\x -> isDummy x == NonDummy) $ tIndices tempT                

            -- If the given tensors have the index, returns 1, otherwise returns 0 
            haveInd :: QTensors -> QIndex -> [Integer]
            haveInd tens ind = fmap (hasInd ind) tens
              where
                hasInd :: QIndex -> QTensor -> Integer
                hasInd ind2 ten  = if ind2 `elem` (tIndices ten) then 1 else 0

        -------------------------------------------------------------
        -- In case of type2 integral
        -- >> Type2 integral :: h2() <- <c1c2|c1c2>
        -------------------------------------------------------------
        caseType2 :: [QIndices] -> Bool
        caseType2 tempI = foldl (||) False $ fmap evalType2 tempI
          where
            evalType2 :: QIndices -> Bool
            evalType2 inds2 = (i0 == i2) && (i1 == i3) && (iSpace i0 == Core) && (iSpace i1 == Core) && (isDummy i0 == Dummy) && (isDummy i1 == Dummy)
              where
                i0 = inds2 !! 0 
                i1 = inds2 !! 1    
                i2 = inds2 !! 2 
                i3 = inds2 !! 3 

        -------------------------------------------------------------
        -- In case of type3 integral
        -- >> Type3 integral :: h3() <- <c1c1|c2c2>
        -------------------------------------------------------------
        caseType3 :: [QIndices] -> Bool
        caseType3 tempI = foldl (||) False $ fmap evalType3 tempI
          where
            evalType3 :: QIndices -> Bool
            evalType3 inds2 = (i0 == i1) && (i2 == i3) && (iSpace i0 == Core) && (iSpace i2 == Core) && (isDummy i0 == Dummy) && (isDummy i2 == Dummy)
              where
                i0 = inds2 !! 0
                i1 = inds2 !! 1
                i2 = inds2 !! 2
                i3 = inds2 !! 3

        -------------------------------------------------------------
        -- In case of type4 integral
        -- >> Type4 integral :: h4(*,*) <- <c1*|c1*>
        -------------------------------------------------------------                    
        caseType4 :: [QIndices] -> Bool
        caseType4 tempI = foldl (||) False $ fmap evalType4 tempI
          where
            evalType4 :: QIndices -> Bool
            evalType4 inds2 = (i0 == i2) && (iSpace i0 == Core) && (isDummy i0 == Dummy) -- && ((isDummy i1 /= Dummy) || (isDummy i3 /= Dummy))
              where
                i0 = inds2 !! 0
                i1 = inds2 !! 1
                i2 = inds2 !! 2
                i3 = inds2 !! 3
            
        -------------------------------------------------------------
        -- In case of type5 integral
        -- >> Type5 integral :: h5(*,*) <- <c1c1|**>
        -------------------------------------------------------------                    
        caseType5 :: [QIndices] -> Bool
        caseType5 tempI = foldl (||) False $ fmap evalType5 tempI
          where
            evalType5 :: QIndices -> Bool
            evalType5 inds2 = (i0 == i1) && (iSpace i0 == Core) && (isDummy i0 == Dummy) -- && ((isDummy i2 /= Dummy) || (isDummy i3 /= Dummy))
              where
                i0 = inds2 !! 0
                i1 = inds2 !! 1
                i2 = inds2 !! 2
                i3 = inds2 !! 3

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Make core Fock matrix by combinimg the core operators
constructCoreFock :: QTerms -> QTerms
constructCoreFock inTerms
  | numReplaced      /= 0 = replaceEffectiveInts replacedTerms
  | numReplacedCFock /= 0 = replacedTerms 
  | otherwise             = inTerms
  where

    ----------------------------------------------------------------------------------
    -- If the number of the replaced integral terms is zero, return inTerms untached
    ----------------------------------------------------------------------------------
    replacedTerms    = fmap makeEffectiveInts inTerms
    numReplaced      = length $ filter (\x -> x == h1Eff_ || x == h2Eff_ || x == h3Eff_ || x == h4Eff_) $ fmap tLabel $ concat $ fmap tTensor replacedTerms
    numReplacedCFock = length $ filter (\x -> x == casEff_) $ concat $ fmap tCoeff replacedTerms  

    -------------------------------------------------------------------------------------------------------------
    -- II. Combine the effective inetgrals to construct the Fock matrix -----------------------------------------
    -------------------------------------------------------------------------------------------------------------
    -- Make sure this algorithm is valid only for the CI, CC and CT ansatz. In general, it cannot be used for
    -- the PT approaches since term of the PT equations may contain more than one integrals.
    -------------------------------------------------------------------------------------------------------------
    replaceEffectiveInts :: QTerms -> QTerms
    replaceEffectiveInts givenTerms = uniReplaceFock $ fmap canonicalizeTerm givenTerms
      where
        --------------------------------------------
        -- Body of the function
        --------------------------------------------
        uniReplaceFock :: QTerms -> QTerms
        uniReplaceFock myTerms = fockTerms ++ filteredTerms
          where
            --------------------------------------------------------------------------
            -- Strategy to construct the core integrals
            --------------------------------------------------------------------------
            -- [1] Extract the terms with H1 and H4 integrals from the other terms
            -- [2] Construct the core integrals fromt the extracted terms
            -- [3] Erase the rest effective integrals from the rest terms
            --------------------------------------------------------------------------
            (targets, extractedTerms) = List.partition (\x -> hasH1 x ||  hasH4 x) myTerms -- H1 and H4 terms are filtered
            fockTerms = fmap (\x -> if hasH1 x then makeEc x else makeFc x) targets        -- Replaced by the core integrals
            filteredTerms = filterInts targets extractedTerms

            -- Takes fock terms to screen the effective integrals
            filterInts :: QTerms -> QTerms -> QTerms
            filterInts [] _ = []
            filterInts _ [] = []
            filterInts (f:fs) toBeFiltered
              | length fs == 0 = if hasH1 f then screenInts "H1" f toBeFiltered else screenInts "H4" f toBeFiltered
              | hasH1 f        = filterInts fs $ screenInts "H1" f toBeFiltered
              | hasH4 f        = filterInts fs $ screenInts "H4" f toBeFiltered
              | otherwise = error "filterInts: Algorithmic error"

            -----------------------------------------------------------------------------------------
            -- Term rewriting functions
            -----------------------------------------------------------------------------------------
            screenInts :: String -> QTerm -> QTerms -> QTerms
            screenInts mode thisTerm terms
              | mode == "H1" = screenh1 thisTerm terms
              | mode == "H4" = screenh4 thisTerm terms
              | otherwise    = error $ "verifyInts: Invalid option specified. >> " ++ mode ++ " << "
              where
                
                --------------------------------------------------------------
                -- Make core energy:
                -- >> Ecore <- 2 h1Int_() + 2 h2Int_() - h3Int_()
                --------------------------------------------------------------
                --------------
                screenh1 :: QTerm -> QTerms -> QTerms
                screenh1 tempH1 tempTs = if extractedNames == List.sort [h2Eff_, h3Eff_]
                                         then alives
                                         else error $ "screenh1: Core energy construction failed.  for " ++ (show tempH1) ++ " >> myJ is: " ++ (show myJ) ++ " myK: " ++ (show myK) ++ " the rests are >> " ++ (show tempTs) ++  " >> Found candidates are: " ++ (show killed)
                  where
                    (h1Int, theRest) = List.partition (\x -> tLabel x == h1Eff_) $ tTensor tempH1                    
                    myJ              = canonicalizeTerm $
                                       baseTerm (tNum tempH1)            (tCoeff tempH1) (baseTensor h2Eff_ (tIndices $ h1Int !! 0) (tSymm $ h1Int !! 0) (tComm $ h1Int !! 0) : theRest)
                    myK              = canonicalizeTerm $
                                       baseTerm (-1.0*(tNum tempH1)/2.0) (tCoeff tempH1) (baseTensor h3Eff_ (tIndices $ h1Int !! 0) (tSymm $ h1Int !! 0) (tComm $ h1Int !! 0) : theRest)

                    (killed, alives) = List.partition (\x -> x == myJ || x == myK) tempTs
                    extractedNames   = List.sort $ filter (\x -> x == h2Eff_ || x == h3Eff_) $ fmap tLabel $ concat $ fmap tTensor killed
                --------------

                --------------------------------------------------------------
                -- Make core Fock matrix:
                -- >> Fcore(p,q) <- h1(p,q) + 2 h4Int_(p,q) - h5Int_(p,q)
                --------------------------------------------------------------
                --------------
                screenh4 :: QTerm -> QTerms -> QTerms
                screenh4 tempH4 tempTs = if extractedNames == List.sort [h1Name_, h5Eff_]
                                         then alives
                                         else error $ "screenh4: Core Fock matrix construction failed. >> " ++ show extractedNames ++ " for " ++ show tempH4
                                              ++ " myH: " ++ show myH ++ " myK: " ++ show myK ++ " the rest terms are: " ++ show tempTs ++ " << "
                  where
                    (h4Int, theRest) = List.partition (\x -> tLabel x == h4Eff_) $ tTensor tempH4                    
                    myH              = canonicalizeTerm $
                                       baseTerm ( 1.0*(tNum tempH4)/2.0) (tCoeff tempH4) (baseOne           (tIndices $ h4Int !! 0)                                           : theRest)
                    myK              = canonicalizeTerm $
                                       baseTerm (-1.0*(tNum tempH4)/2.0) (tCoeff tempH4) (baseTensor h5Eff_ (tIndices $ h4Int !! 0) (tSymm $ h4Int !! 0) (tComm $ h4Int !! 0) : theRest)

                    (killed, alives) = List.partition (\x -> x == myH || x == myK) tempTs
                    extractedNames   = List.sort $ filter (\x -> x == h1Name_ || x == h5Eff_) $ fmap tLabel $ concat $ fmap tTensor killed
                --------------
                    
            -- Returns the replaced term by the core energy
            makeEc :: QTerm -> QTerm
            makeEc term = if h1Eff_ `elem` (fmap tLabel $ tTensor term) then replacedTerm else error "makeEc: The given term doesn't have any type1 integral."
              where
                (h1Int, theRest) = List.partition (\x -> tLabel x == h1Eff_) $ tTensor term                
                replacedTerm     = baseTerm ((tNum term)/2.0) (ecEff_ : tCoeff term) theRest 

            -- Returns the replaced term by the core Fock matrix
            makeFc :: QTerm -> QTerm
            makeFc term = if h4Eff_ `elem` (fmap tLabel $ tTensor term) then replacedTerm else error "makeFc: The given term doesn't have any type4 integral."
              where
                (h4Int, theRest) = List.partition (\x -> tLabel x == h4Eff_) $ tTensor term                
                replacedTerm     = baseTerm ((tNum term)/2.0) (tCoeff term) $ baseTensor fcEff_ (tIndices $ h4Int !! 0) (tSymm $ h4Int !! 0) (tComm $ h4Int !! 0): theRest 

            -- Returns True if the term has the h1 integral
            hasH1 :: QTerm -> Bool
            hasH1 term = h1Eff_ `elem` (fmap tLabel $ tTensor term)

            -- Returns True if the term has the h1 integral
            hasH4 :: QTerm -> Bool
            hasH4 term = h4Eff_ `elem` (fmap tLabel $ tTensor term)
            
