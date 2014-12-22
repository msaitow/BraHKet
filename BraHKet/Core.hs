-- ///////////////////////////////////////////////////////////////////////
--  BraHKet :: <H> A Haskell Module For Fermionic Many-Body Operators <H>
--                            Masaaki Saitow
-- ///////////////////////////////////////////////////////////////////////
--
--                                      ,--,         
--                                    ,--.'|         
--                          ,--.   ,--,  | :.--,     
--                         /  /|,---.'|  : '|\  \    
--                        '  / '|   | : _' |` \  `   
--                       /  / / :   : |.'  | \ \  \  
--                      /  / ,  |   ' '  ; :  , \  \ 
--                      \ '\ \  '   |  .'. |  / /` / 
--                       \  \ ' |   | :  | ' ` /  /  
--                        \  . |'   : |  : ;| .  /   
--                         \__\.|   | '  ,/ ./__/    
--                              ;   : ;--'           
--                              |   ,/               
--                              '---'                
--
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
-- ///////////////////////////////////////////////////////////////////////////

----------------------------------
-- Main module -------------------
----------------------------------
module BraHKet.Core
(
  -- Acronyms
  Permut,
  Coeffs,
  QIndices,
  QTensors,
  QTerms,
  
  -- Predefined names for various tensors
  creName_,     -- Creation operators
  desName_,     -- Destruction operators 
  sfGenName_,   -- Spin-free generators
  sfRDMName_,   -- Spin-free density
  sdGenName_,   -- Spin-dependent generator
  kDName_,      -- Kronecker's delta
  h1Name_,      -- One-body int
  h2Name_,      -- Two-body int
  casFockName_, -- CAS Fock Matrix

  -- Functions for QIndex data
  QIndex(..),
  alphaIndex,   -- Alpha-spin orbital index
  betaIndex,    -- Beta-spin orbital index
  sfIndex,      -- Spin-free index
  soIndex,      -- Spin-orbital index   
  
  -- Functions for QTensor data
  QTensor(..),
  baseTensor,  -- Generic tensor
  sfTensor,    -- Spin-free tensor
  baseTensorP, -- Generic tensor
  baseOne,     -- One-body int
  baseKD,      -- Kronecker's delta
  baseERI,     -- Two-body int
  baseSFGen,   -- Spin-free generator
  baseSFRDM,   -- Spin-free density-matrix
  baseCFock,   -- CAS Fock Matrix
  baseSDGen,   -- Spin-dependent generator  
  baseCre,     -- Creation operator
  baseDes,     -- Destruction operator
  tsortIndices,
  tallPermInds,
  tgetConfs,
  
  -- Functions for QTerm data
  QTerm(..),
  baseTerm,            -- Generic term
  getSummedBody,
  isAdditive,
  masquerade,
  allRenames,
  rotateTensors,
  rotateAllIndices,
  generateAllConfs,
  killKDeltas,
  combineTerms,
  canonicalizeTerm,  
  generateInteractions, -- Decompose generic indices into the interactions
  
  ----------------------------------------
  -- Auxiliary data
  ----------------------------------------
  QSpin(..),
  QSum(..),    
  QNature(..),     
  QSpace(..)
  ----------------------------------------
  
) where

import qualified BraHKet.Utils as Utils 
import qualified Data.Map as Map
import qualified Data.Maybe as Maybe
import qualified Data.List as List
--import Control.Applicative

-------------------------------------------------------------------------
-- Acronyms
-------------------------------------------------------------------------
type Permut   = [[Int]]   -- Permutation Symmetry
type Coeffs   = [String]  -- Coefficients
type QIndices = [QIndex]
type QTensors = [QTensor]
type QTerms   = [QTerm]

-------------------------------------------------------------------------
-- @@ Tensor names
-------------------------------------------------------------------------
creName_     = "@Cre" :: String -- Creation operator
desName_     = "@Des" :: String -- Destruction operator 
sfGenName_   = "@E"   :: String -- Spin-free generator
sfRDMName_   = "@D"   :: String -- Spin-free density
sdGenName_   = "@G"   :: String -- Spin-dependent generator
kDName_      = "@kD"  :: String -- Kronecker's delta
h1Name_      = "@h1"  :: String -- One-body int
h2Name_      = "@V2"  :: String -- Two-body int
casFockName_ = "@cF"  :: String -- CAS Fock Matrix

-- ///////////////////////////////////////////////////////////////////////
-- // Fundamental classes to construct the fermionic many-body formulas
-- ///////////////////////////////////////////////////////////////////////

-------------------------------------------------------------------------
-- @@ Data that represent the spin of the orbital index
-------------------------------------------------------------------------
data QSpin = Alpha | Beta | SpinFree | SpinOrbital deriving(Show, Eq, Ord)
  
-------------------------------------------------------------------------
-- @@ Data that represents summation over the index
-------------------------------------------------------------------------
data QSum = Dummy | NonDummy deriving(Show, Eq, Ord)

-------------------------------------------------------------------------
-- @@ Data that represents whether the tensor is an operator 
-------------------------------------------------------------------------
data QNature = Operator | Classical deriving(Show, Eq, Ord)

-------------------------------------------------------------------------
-- @@ Orbital space
-------------------------------------------------------------------------
data QSpace = Core | Active | Virtual | Generic deriving(Show, Eq, Ord)

-------------------------------------------------------------------------             
-- @@ Index class
-------------------------------------------------------------------------
data QIndex = QIndex { iLabel   :: String,
                       iSpace   :: QSpace,
                       iSpin    :: QSpin,
                       isDummy  :: QSum
                     } deriving (Eq, Ord)

instance Show QIndex where
  show ind = "@(" ++ iLabel ind ++ ", " ++ show (iSpace ind) ++ ", " ++ show (iSpin ind) ++ ", " ++ if isDummy ind == Dummy then "Anonymous)" else "Universal)" 

-- Standard constructor for the alpha-spin orbital index
alphaIndex :: String -> QSpace -> QSum -> QIndex
alphaIndex label space dummy = QIndex label space Alpha dummy

-- Standard constructor for the beta-spin orbital index
betaIndex :: String -> QSpace -> QSum -> QIndex
betaIndex label space dummy = QIndex label space Beta dummy

-- Standard constructor for the spin-free orbital index
sfIndex :: String -> QSpace -> QSum -> QIndex
sfIndex label space dummy = QIndex label space SpinFree dummy

-- Standard constructor for the spin-free orbital index
soIndex :: String -> QSpace -> QSum -> QIndex
soIndex label space dummy = QIndex label space SpinOrbital dummy

-------------------------------------------------------------------------
-- @@ Tensor class
-------------------------------------------------------------------------  
data QTensor = QTensor { tLabel   :: String,   -- Name
                         tIndices :: QIndices, -- Indices
                         tSymm    :: Permut,   -- Permutation
                         tComm    :: QNature   -- Whether commutable or not
                       } deriving (Eq, Ord)

showTensor :: QTensor -> String
showTensor ten =
  let
    spinLabels = fmap repSpin $ tIndices ten
    indLabels1 = map (\x -> x ++ ",") $ zipWith (++) (fmap iLabel $ tIndices ten) spinLabels
    indLabels2 = foldl (++) [] indLabels1
    indLabels3 = take ((length indLabels2)-1) $ indLabels2
  in (tLabel ten ++ "(" ++ indLabels3 ++ ")")
     where
       repSpin :: QIndex -> String
       repSpin myInd
         | iSpin myInd == Alpha       = "(a)"
         | iSpin myInd == Beta        = "(b)"
         | iSpin myInd == SpinOrbital = "_"                                     
         | iSpin myInd == SpinFree    = ""
         | otherwise                  = error $ "showTensor: Can't handle this spin-state. >> " ++ (show $ iSpin myInd) ++ " << "

instance Show QTensor where
  show = showTensor

-- Standard constructor
baseTensor :: String -> QIndices -> Permut -> QNature -> QTensor
baseTensor label indices symm iscomm
  | not lenSymm = error $ "QTensor: Lengths of given indices and symm mismatched. >> Label: " ++ label ++ " Indices: " ++ (show $ indices) ++ " << "
  | otherwise = QTensor label indices mySymm iscomm
  where 
    mySymm  = if length symm == 0 then [[0..(length indices)-1]] else symm  
    lenSymm =
          let lenInds = take (length symm) (repeat $ length indices)
          in lenInds == map length symm

-- Standard constructor for the spin-free tensor
sfTensor :: String -> QIndices -> Permut -> QNature -> QTensor
sfTensor label indices symm iscomm = if allSF then baseTensor label indices symm iscomm else error "sfTensor: The orbital indices other than spin-free type detected."
  where
    allSF = foldl (&&) True $ fmap (\x -> iSpin x == SpinFree) indices

-- Standard constructor that calls makePerms internally
baseTensorP :: String -> QIndices -> Permut -> QNature -> QTensor
baseTensorP label indices symm iscomm
  | not lenSymm = error "QTensor: Lengths of given indices and symm mismatched."
  | otherwise = QTensor label indices symm iscomm
  where 
    mySymm  = if length symm == 0 then [[0..(length indices)-1]] else symm  
    lenSymm =   
      let
        allSymm = Utils.makePerms symm
        lenInds = take (length allSymm) (repeat $ length indices)
        lenSymm = map length allSymm
      in lenInds == lenSymm

-- Constructor for the one-body integral tensor
baseOne :: QIndices -> QTensor
baseOne inds
  | length inds /= 2 = error "QTensor: Length of indices for kinetic operator should be 2."
  | otherwise = QTensor h1Name_ inds eriPerms Classical
  where eriPerms =
          let hSource = [[0,1],[1,0]]
          in Utils.makePerms hSource

-- Constructor for the Kronecker's delta
baseKD :: QIndices -> QTensor
baseKD inds
  | length inds /= 2 = error "QTensor: Length of indices for Kronecker's delta shoud be 2."
  | otherwise = QTensor kDName_ inds kDPerm Classical
  where kDPerm = [[0,1],[1,0]]

-- Constructor for the two-body integral tensor
baseERI :: QIndices -> QTensor
baseERI inds
  | length inds /= 4 = error "QTensor: Length of indices for ERI should be 4."
  | otherwise = QTensor h2Name_ inds eriPerms Classical
  where eriPerms =
          let eriSource = [[0,1,2,3],[2,1,0,3],[0,3,2,1],[1,0,3,2]]
          in Utils.makePerms eriSource

-- Constructor for the spin-free unitary group generator
baseSFGen :: QIndices -> QTensor
baseSFGen inds
  | allSF /= True     = error "QTensor: Spin-dependent indices detected in the baseSFGen."
  | odd $ length inds = error "QTensor: Number of indices given to the spin-free generator should be even number."
  | otherwise         = QTensor name inds genPerm Operator
  where
    allSF = foldl (&&) True $ fmap (\x -> iSpin x == SpinFree) inds    
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sfGenName_ ++ (show order)
    genPerm = Utils.makePermSFGen order

-- Constructor for the spin-dependent unitary group generator
baseSDGen :: QIndices -> QTensor
baseSDGen inds
  | hasSF             = error "QTensor: Spin-free indices detected in the baseSDGen."
  | odd $ length inds = error "QTensor: Number of indices given to the spin-dependent generator should be even number."
  | otherwise = QTensor name inds genPerm Operator
  where
    hasSF = foldl (||) False $ fmap (\x -> iSpin x == SpinFree) inds        
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sdGenName_ ++ (show order)
    genPerm = Utils.makePermSFGen order

-- Constructor for the spin-free reduced-density matrix
baseSFRDM :: QIndices -> QTensor
baseSFRDM inds
  | allSF /= True     = error "QTensor: Spin-dependent indices detected in the baseSFGen."
  | odd $ length inds = error "QTensor: Number of indices given to the spin-free density matrix should be even number."
  | otherwise = QTensor name inds genPerm Classical
  where
    allSF = foldl (&&) True $ fmap (\x -> iSpin x == SpinFree) inds        
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sfRDMName_ ++ (show order)
    genPerm = Utils.makePermSFRDM order

-- Constructor for the CAS Fock matrix
baseCFock :: (QIndex, QIndex) -> QTensor
baseCFock (ind1, ind2) = QTensor casFockName_ [ind1, ind2] [[0, 1], [1, 0]] Classical 

-- Constructor for the creation operator
baseCre :: QIndex -> QTensor
baseCre ind = QTensor creName_ [ind] [[0]] Operator

-- Constructor for the destruction operator
baseDes :: QIndex -> QTensor
baseDes ind = QTensor desName_ [ind] [[0]] Operator

-- Returns normal-ordered QTensor
tsortIndices :: QTensor -> QTensor
tsortIndices t = QTensor (tLabel t) sortedIndices (tSymm t) (tComm t)
  where 
    sortedIndices = maximum $ fmap (uniPermRR $ tIndices t) $ tSymm t
    uniPermRR x y = Utils.uniPerm y x

-- Returns all possible indices pattern
tallPermInds :: QTensor -> [QIndices]
tallPermInds t = fmap (Utils.uniPermR $ tIndices t) $ tSymm t

-- Return all possible configuration for this QTensor object
tgetConfs :: QTensor -> QTensors
tgetConfs t = fmap (makeTensor) (tallPermInds t)
  where
    makeTensor x = baseTensor (tLabel t) x (tSymm t) (tComm t)

-------------------------------------------------------------------------
-- @@ Term class
-------------------------------------------------------------------------
data QTerm = QTerm { tNum    :: Double,
                     tCoeff  :: Coeffs,
                     tTensor :: QTensors
                   } deriving (Eq, Ord)

instance Show QTerm where
  show term =
    let tensors = foldl (++) [] $ fmap (\x -> x ++ " ") $ fmap show $ tTensor term
        coeffs  = if length (tCoeff term) == 0 then "" else (foldl (++) [] (fmap (++" ") $ tCoeff term))
        mySign  = if (tNum term) >= 0.0 then "+" else ""
    in mySign ++ (show $ tNum term) ++ " " ++ coeffs ++ take ((length tensors)-1) tensors

-- Constructor for the QTerm class
baseTerm :: Double -> Coeffs -> QTensors -> QTerm
baseTerm num coeff tensors = QTerm num coeff (commutatives ++ incommutatives)
  where
    commutatives   = filter (\x -> (tComm x) == Classical) tensors
    incommutatives = filter (\x -> (tComm x) == Operator ) tensors

-- Returns set of all the indices of the term
getSummedBody :: QTerm -> QIndices
getSummedBody term = List.nub . concat $ map (tIndices) $ (tTensor term)

-- Returns two terms are additive or nor
isAdditive :: QTerm -> QTerm -> Bool
isAdditive kore are = ((tCoeff kore) == (tCoeff are)) && ((tTensor kore) == (tTensor are))

-- Rename all the dummy indices
masquerade :: QTerm -> QTerm
masquerade t = QTerm (tNum t) (tCoeff t) ten
  where
    getSummedInds term =
      let
        thisInds   = getSummedBody term
        mapCore    = mapBase Core    thisInds 0
        mapActive  = mapBase Active  thisInds 0
        mapVirtual = mapBase Virtual thisInds 0
        mapGeneric = mapBase Generic thisInds 0        
      in Map.fromList $ mapCore ++ mapActive ++ mapVirtual ++ mapGeneric

    changedList = getSummedInds t
    changeInds tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
      where
        inds = fmap (replaceInds) $ tIndices tenten
        replaceInds i = if Map.lookup i changedList == Nothing then i else QIndex (Maybe.fromJust $ Map.lookup i changedList) (iSpace i) (iSpin i) (isDummy i)
        
    ten = fmap (changeInds) (tTensor t)

-- Returns all the possible configurations for a given term
allRenames :: QTerm -> QTerms
allRenames term = zipWith (replaceTerms) allPattern dummyTerms
  where
    patterns = [x ++ y ++ z ++ v| x <- cAll, y <- aAll, z <- vAll, v <- gAll]
    allPattern = fmap (Map.fromList) patterns
    cAll = allMapBase Core    $ getSummedBody term
    aAll = allMapBase Active  $ getSummedBody term
    vAll = allMapBase Virtual $ getSummedBody term
    gAll = allMapBase Generic $ getSummedBody term    
    dummyTerms = take (length allPattern) $ repeat term

    replaceTerms pattern t = QTerm (tNum t) (tCoeff t) ten
      where
        changeInds tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
          where
            inds = fmap (replaceInds) $ tIndices tenten
            replaceInds i = if Map.lookup i pattern == Nothing then i else QIndex (Maybe.fromJust $ Map.lookup i pattern) (iSpace i) (iSpin i) (isDummy i)
  
        ten = fmap (changeInds) (tTensor t)

-- Returns all the possibly ordered tensors (Make sure this algorithm swaps the operators as well the c-tensors)    
rotateTensors :: QTerm -> QTerms
rotateTensors term = fmap makeTerms tensors
  where
    classicalTensors = filter (\x -> tComm x == Classical) $ tTensor term
    operatorTensors  = filter (\x -> tComm x == Operator)  $ tTensor term    
    tensors = fmap (\x -> x ++ operatorTensors) $ List.permutations classicalTensors
    makeTerms t = QTerm (tNum term) (tCoeff term) t

-- Rotate all the indices of each tensor in a given term
rotateAllIndices :: QTerm -> QTerms
rotateAllIndices term = fmap makeTerms $ Utils.makeCombi tensors
  where
    tensors = fmap tgetConfs $ tTensor term
    makeTerms t = QTerm (tNum term) (tCoeff term) t

-- Generate all the possible configurations of a given term
generateAllConfs :: QTerm -> QTerms
generateAllConfs term = concat . fmap (rotateAllIndices) $ concat . fmap (rotateTensors) $ allRenames term

-- Kill Kronecker's delta (Still buggy??)
-- In the current implementation, the SpinOrbital indices cannot be handled.
killKDeltas :: QTerm -> Maybe QTerm
killKDeltas term
  | hasSO              = error "killKDeltas: SpinOrbital index detected. Please decompose such indices into the other types of the spin indices first."
  | length kdList == 0 = Just term
  | otherwise          = if foldr (||) False $ fmap isZero kdList then Nothing else Just (QTerm (tNum term) (tCoeff term) exceptkDeltas)
  where
    hasSO = foldl (||) False $ fmap (\x -> iSpin x == SpinOrbital) $ concat . fmap tIndices $ tTensor term
      
    -- First, construct a Map [(killed Index, killer Index)]
    tensors = tTensor term
    kdList  = [x | x <- tensors, tLabel x == kDName_]

    isZero :: QTensor -> Bool
    isZero kd =
      let
        ind1   = (tIndices kd) !! 0
        ind2   = (tIndices kd) !! 1
        space1 = iSpace ind1
        space2 = iSpace ind2
        spin1  = iSpin ind1
        spin2  = iSpin ind2        
      in space1 /= space2 && space1 /= Generic && space2 /= Generic || spin1 /= spin2

    mkMap   = Map.fromList $ fmap killSchedule kdList
    killSchedule kDelta
      | (isDummy ind1 ==    Dummy) && (isDummy ind2 ==    Dummy) && (iSpace ind1 == Generic) = (ind1, ind2)
      | (isDummy ind1 ==    Dummy) && (isDummy ind2 ==    Dummy) && (iSpace ind2 == Generic) = (ind2, ind1)
      | (isDummy ind1 ==    Dummy) && (isDummy ind2 == NonDummy) && (iSpace ind1 == Generic) = (ind1, ind2)
      | (isDummy ind1 == NonDummy) && (isDummy ind2 ==    Dummy) && (iSpace ind2 == Generic) = (ind2, ind1)
      | (isDummy ind1 ==    Dummy) && (isDummy ind2 == NonDummy)                             = (ind1, ind2)
      | (isDummy ind1 == NonDummy) && (isDummy ind2 ==    Dummy)                             = (ind2, ind1)
      | (isDummy ind1 ==    Dummy) && (isDummy ind2 ==    Dummy)                             = (ind1, ind2)
      | (isDummy ind1 == NonDummy) && (isDummy ind2 == NonDummy) && (iSpace ind1 == Generic) = (ind1, ind2)
      | (isDummy ind1 == NonDummy) && (isDummy ind2 == NonDummy) && (iSpace ind2 == Generic) = (ind2, ind2)
      | (isDummy ind1 == NonDummy) && (isDummy ind2 == NonDummy)                             = (ind2, ind2)
      | otherwise = error "killKDeltas: Algorithmic error"
      where
        ind1 = (tIndices kDelta) !! 0
        ind2 = (tIndices kDelta) !! 1
    -- Then, construct new tensors 
    newTensors = fmap (Utils.repeatN replaceTensor (length kdList)) tensors
    replaceTensor tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
      where
        inds = fmap (replaceInds) $ tIndices tenten
        replaceInds i = if Map.lookup i mkMap == Nothing then i  else Maybe.fromJust $ Map.lookup i mkMap
    exceptkDeltas = [x | x <- newTensors, not (tLabel x == kDName_ && (tIndices x) !! 0 == (tIndices x) !! 1)]

-- Function to return the combined terms. This can be implemented by using State monad more elegantly?
combineTerms :: QTerms -> QTerms
combineTerms terms = filter (\x -> (tNum x) /= 0.0) $ fmap (changeFactors) origTerms
  where
    maxTerms = fmap canonicalizeTerm terms
    origTerms = List.nub $ fmap (makeUniTerm) maxTerms

    makeUniTerm :: QTerm -> QTerm 
    makeUniTerm kore = QTerm 1.0 (tCoeff kore) (tTensor kore)

    changeFactors :: QTerm -> QTerm
    changeFactors are = baseTerm myFactor (tCoeff are) (tTensor are)
      where myFactor = collectFactors are maxTerms

    collectFactors :: QTerm -> QTerms -> Double
    collectFactors kore korera = 
      let
        arera = [x | x <- korera, isAdditive kore x]
      in sum $ fmap (tNum) arera

-- Canonicalize the representation of the QTerm object
canonicalizeTerm :: QTerm -> QTerm
canonicalizeTerm konoTerm = if length possibility /= 0 then baseTerm newFactor (tCoeff konoTerm) (tTensor myMax) else konoTerm
  where
    newFactor   = (tNum konoTerm) * (tNum myMax)
    myMax       = maximum possibility
    possibility = normalizeTerms konoTerm

    -- Body of the function
    normalizeTerms :: QTerm -> QTerms
    normalizeTerms thisTerm
      | length koreraTensor == (length $ List.nub koreraTensor) = concat . fmap (rotateAllIndices) $ allRenames canonicalConf
      | otherwise                                               = generateAllConfs canonicalConf
      where
        koreraTensor    = fmap tLabel $ tTensor thisTerm
        canonicalTensor = (List.sortBy (\x y -> (tLabel x) `compare` (tLabel y)) $ filter (\x -> tComm x == Classical) $ tTensor thisTerm) ++ (filter (\x -> tComm x == Operator) $ tTensor thisTerm)
        canonicalConf   = baseTerm (1.0) [] canonicalTensor 
        
---------------------------------------------------------------------------------------------
-- Small utilities
---------------------------------------------------------------------------------------------
    
-- The generalized function for above all the functions
allMapBase :: QSpace -> QIndices -> [[(QIndex, String)]]
allMapBase _ [] = []
allMapBase myLabel inds = zipWith (zipWith (\x y -> (x,y))) allFst allRenames
  where 
    zipFst     = fmap fst $ mapBase myLabel inds 0
    allRenames = List.permutations $ fmap snd $ mapBase myLabel inds 0
    allFst     = take (length allRenames) $ repeat zipFst
    
-- The generalized function for above all the functions
mapBase :: QSpace -> QIndices -> Int -> [(QIndex, String)]
mapBase _ [] _ = []
mapBase mySpace (i:is) num
  | (iSpace i) == mySpace && (isDummy i) == Dummy = [(i, myLabel ++ show num)] ++ (mapBase mySpace is (num+1))
  | otherwise                                     = mapBase mySpace is num
  where
    myLabel = returnLabel mySpace
    
    -- Replacement label for the dummy indices defined here
    returnLabel :: QSpace -> String
    returnLabel n
      | n == Core    = "c"
      | n == Active  = "a"
      | n == Virtual = "v"
      | n == Generic = "g"
      | otherwise    = error $ "mapBase: Can't handle this space >> " ++ (show n) ++ " << "
                       
    
---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Decompose the generic indices into various orbitals indices
generateInteractions :: [QSpace] -> QTerm -> QTerms
generateInteractions targetSpaces thisTerm
  | length genInds == 0 = [thisTerm]
  | otherwise           = fmap replaceSpace allSpaces 
    where
      tensors = tTensor thisTerm
      indices = getSummedBody thisTerm
      genInds = filter (\x -> iSpace x == Generic) indices
      numGens = length genInds

      -- All possible combinations of MO spaces, [[C,C,C], [C,C,A], [C,A,C], ..]
      allSpaces = List.nub $ concat $ fmap List.permutations $ List.nub . filter (\x -> length x == numGens) $ List.subsequences $ concat $ fmap (\x -> take numGens $ repeat x) targetSpaces

      -- Replace the MO space of the term according to the given combinations of the QSpace
      replaceSpace :: [QSpace] -> QTerm
      replaceSpace mySpaces = QTerm (tNum thisTerm) (tCoeff thisTerm) replacedTensors
        where
          killerIndices  = zipWith uniReplace mySpaces genInds
          killSchedule   = Map.fromList $ zip genInds killerIndices -- [(killed generic index, killer index)]
          replacedTensors = fmap replaceTensor tensors

          -- Unit replacement function
          uniReplace :: QSpace -> QIndex -> QIndex
          uniReplace space index = QIndex (iLabel index) space (iSpin index) (isDummy index)

          -- Replace given tensor's indices according to the kill schedule
          replaceTensor :: QTensor -> QTensor
          replaceTensor myTensor = QTensor (tLabel myTensor) replacedInds (tSymm myTensor) (tComm myTensor)
            where
              replacedInds = fmap (replaceIndex) (tIndices myTensor)

              -- Replace the given index if the index should be replaced
              replaceIndex :: QIndex -> QIndex
              replaceIndex givenInd
                | newInd == Nothing = givenInd
                | otherwise         = Maybe.fromJust newInd
                where newInd = Map.lookup givenInd killSchedule


