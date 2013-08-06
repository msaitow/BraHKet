-- ///////////////////////////////////////////////////////////////////////
--  BraHKet :: <H> A Haskell Module For Fermionic Many-Body Operators <H>
--                            Masaaki Saitow
-- ///////////////////////////////////////////////////////////////////////
--------------------------------------------------------------------------
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
--------------------------------------------------------------------------

module BraHKet.Core
(
  -- Predefined names for various tensors
  creName_,   -- Creation operators
  desName_,   -- Destruction operators 
  sfGenName_, -- Spin-free generators
  sfRDMName_, -- Spin-free density
  sdGenName_, -- Spin-dependent generator
  kDName_,    -- Kronecker's delta
  h1Name_,    -- One-body int
  h2Name_,    -- Two-body int

  -- Functions for QSpace data
  QSpace(..),

  -- Functions for QIndex data
  QIndex(..),

  -- Functions for QTensor data
  QTensor(..),
  baseTensor,  -- Generic tensor
  baseOne,     -- One-body int
  baseKD,      -- Kronecker's delta
  baseERI,     -- Two-body int
  baseSFGen,   -- Spin-free generator
  baseSFRDM,   -- Spin-free density-matrix
  baseSDGen,   -- Spin-dependent generator  
  baseCre,     -- Creation operator
  baseDes,     -- Destruction operator
  tsortIndices,
  tallPermInds,
  tgetConfs,

  -- Functions for QTerm data
  QTerm(..),
  baseTerm,      -- Generic term
  getSummedBody,
  isAdditive,
  masquerade,
  allRenames,
  rotateTensors,
  rotateAllIndices,
  generateAllConfs,
  killKDeltas,
  combineTerms,
  normalOrderA, -- NO for Cre/Des operators
  normalOrderG, -- NO for spin-dependent generator
  normalOrderE, -- NO for spin-free generator
  takeVEV
) where

import qualified BraHKet.Utils as Utils 
import qualified Data.Map as Map
import qualified Data.Maybe as Maybe
import qualified Data.List as List
import Control.Applicative

-------------------------------------------------------------------------
-- Acronyms
-------------------------------------------------------------------------
type Permut   = [[Int]]   -- Permutation Symmetry
type Coeffs   = [String]  -- Coefficients
type QIndices = [QIndex]
type QTensors = [QTensor]
type QTerms   = [QTerm]

-- ///////////////////////////////////////////////////////////////////////
-- // Fundamental classes to construct the fermionic many-body formulas
-- ///////////////////////////////////////////////////////////////////////

-------------------------------------------------------------------------
-- @@ Tensor names
-------------------------------------------------------------------------
creName_   = "@Cre" :: String -- Creation operator
desName_   = "@Des" :: String -- Destruction operator 
sfGenName_ = "@E"   :: String -- Spin-free generator
sfRDMName_ = "@D"   :: String -- Spin-free density
sdGenName_ = "@G"   :: String -- Spin-dependent generator
kDName_    = "@kD"  :: String -- Kronecker's delta
h1Name_    = "@h1"  :: String -- One-body int
h2Name_    = "@V2"  :: String -- Two-body int

-------------------------------------------------------------------------
-- @@ Orbital space
-------------------------------------------------------------------------
data QSpace = Core | Active | Virtual | Generic deriving(Show, Eq, Ord)

-------------------------------------------------------------------------                                                        
-- @@ Index class
-------------------------------------------------------------------------
data QIndex = QIndex { iLabel   :: String,
                       iSpace   :: QSpace,
                       isDummy  :: Bool
                     } deriving (Eq, Ord)

instance Show QIndex where
  show ind = "@(" ++ iLabel ind ++ ", " ++ show (iSpace ind) ++ ", " ++ if isDummy ind == True then "Anonymous)" else "Universal)" 

-------------------------------------------------------------------------
-- @@ Tensor class
-------------------------------------------------------------------------  
data QTensor = QTensor { tLabel   :: String,   -- Name
                         tIndices :: QIndices, -- Indices
                         tSymm    :: Permut,   -- Permutation
                         tComm    :: Bool      -- Whether commutable or not
                       } deriving (Eq, Ord)
                        
showTensor :: QTensor -> String
showTensor ten =
  let all1 = map (\x -> x ++ ",") $ foldr (:) [] $ fmap iLabel $ tIndices ten
      all2 = foldl (++) [] all1
      all3 = take ((length all2)-1) $ all2
  in (tLabel ten ++ "(" ++ all3 ++ ")")

instance Show QTensor where
  show = showTensor

-- Standard constructor
baseTensor :: String -> QIndices -> Permut -> Bool -> QTensor
baseTensor label indices symm iscomm
  | not lenSymm = error "QTensor: Lengths of given indices and symm mismatched."
  | otherwise = QTensor label indices symm iscomm
  where lenSymm =
          let
            allSymm = Utils.makePerms symm
            lenInds = take (length allSymm) (repeat $ length indices)
            lenSymm = map length allSymm
          in lenInds == lenSymm

-- Constructor for the one-body integral tensor
baseOne :: QIndices -> QTensor
baseOne inds
  | length inds /= 2 = error "QTensor: Length of indices for kinetic operator should be 2."
  | otherwise = QTensor h1Name_ inds eriPerms True
  where eriPerms =
          let hSource = [[0,1],[1,0]]
          in Utils.makePerms hSource

-- Constructor for the Kronecker's delta
baseKD :: QIndices -> QTensor
baseKD inds
  | length inds /= 2 = error "QTensor: Length of indices for Kronecker's delta shoud be 2."
  | otherwise = QTensor kDName_ inds kDPerm True
  where kDPerm = [[0,1],[1,0]]

-- Constructor for the two-body integral tensor
baseERI :: QIndices -> QTensor
baseERI inds
  | length inds /= 4 = error "QTensor: Length of indices for ERI should be 4."
  | otherwise = QTensor h2Name_ inds eriPerms True
  where eriPerms =
          let eriSource = [[0,1,2,3],[2,1,0,3],[0,3,2,1],[1,0,3,2]]
          in Utils.makePerms eriSource

-- Constructor for the spin-free unitary group generator
baseSFGen :: QIndices -> QTensor
baseSFGen inds
  | odd $ length inds = error "QTensor: Number of indices given to the spin-free generator should be even number."
  | otherwise = QTensor name inds genPerm False
  where
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sfGenName_ ++ (show order)
    genPerm = Utils.makePermSFGen order

-- Constructor for the spin-dependent unitary group generator
baseSDGen :: QIndices -> QTensor
baseSDGen inds
  | odd $ length inds = error "QTensor: Number of indices given to the spin-dependent generator should be even number."
  | otherwise = QTensor name inds genPerm False
  where
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sdGenName_ ++ (show order)
    genPerm = Utils.makePermSFGen order

-- Constructor for the spin-free reduced-density matrix
baseSFRDM :: QIndices -> QTensor
baseSFRDM inds
  | odd $ length inds = error "QTensor: Number of indices given to the spin-free density matrix should be even number."
  | otherwise = QTensor name inds genPerm True
  where
    order = floor $ (fromIntegral $ length inds) / (fromIntegral 2)
    name  = sfRDMName_ ++ (show order)
    genPerm = Utils.makePermSFRDM order

-- Constructor for the creation operator
baseCre :: QIndex -> QTensor
baseCre ind = QTensor creName_ [ind] [[0]] False

-- Constructor for the destruction operator
baseDes :: QIndex -> QTensor
baseDes ind = QTensor desName_ [ind] [[0]] False

-- -- Returns normal-ordered QTensor
-- tsortIndices :: QTensor -> QTensor
-- tsortIndices t = QTensor (tLabel t) sortedIndices (tSymm t) (tComm t)
--   where 
--     sortedIndices = maximum $ fmap (Utils.uniPermR $ tIndices t) $ tSymm t

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

-- Messy -- instance Show QTerm where
-- Messy --   show term =
-- Messy --     let tensors = foldl (++) [] $ fmap (\x -> x ++ " ") $ fmap show $ tTensor term
-- Messy --         coeffs  = if length (tCoeff term) == 0 then "" else (foldl (++) [] (fmap show $ tCoeff term)) ++ " "
-- Messy --     in (show $ tNum term) ++ " " ++ coeffs ++ take ((length tensors)-1) tensors

instance Show QTerm where
  show term =
    let tensors = foldl (++) [] $ fmap (\x -> x ++ " ") $ fmap show $ tTensor term
        coeffs  = if length (tCoeff term) == 0 then "" else (foldl (++) [] (fmap (++" ") $ tCoeff term))
    in (show $ tNum term) ++ " " ++ coeffs ++ take ((length tensors)-1) tensors

-- Constructor for the QTerm class
baseTerm :: Double -> Coeffs -> QTensors -> QTerm
baseTerm num coeff tensors = QTerm num coeff (commutatives ++ incommutatives)
  where
    commutatives   = filter (tComm) tensors
    incommutatives = filter (\x -> not $ tComm x) tensors

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
        mapCore    = cMap thisInds 0
        mapActive  = aMap thisInds 0
        mapVirtual = vMap thisInds 0
        mapGeneric = gMap thisInds 0
      in Map.fromList $ mapCore ++ mapActive ++ mapVirtual ++ mapGeneric

    changedList = getSummedInds t
    changeInds tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
      where
        inds = fmap (replaceInds) $ tIndices tenten
        replaceInds i = if Map.lookup i changedList == Nothing then i else QIndex (Maybe.fromJust $ Map.lookup i changedList) (iSpace i) (isDummy i)
        
    ten = fmap (changeInds) (tTensor t)

-- bug? -- -- Rename all the dummy indices
-- bug? -- masquerade :: QTerm -> QTerm
-- bug? -- masquerade t = QTerm (tNum t) (tCoeff t) ten
-- bug? --   where
-- bug? --     getSummedInds tenten =
-- bug? --       let
-- bug? --         thisInds   = tIndices tenten
-- bug? --         mapCore    = cMap (tIndices tenten) 0
-- bug? --         mapActive  = aMap (tIndices tenten) 0
-- bug? --         mapVirtual = vMap (tIndices tenten) 0
-- bug? --         mapGeneric = gMap (tIndices tenten) 0
-- bug? --       in Map.fromList $ mapCore ++ mapActive ++ mapVirtual ++ mapGeneric
-- bug? --          
-- bug? --     changeInds tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
-- bug? --       where
-- bug? --         allIndices = getSummedInds tenten
-- bug? --         inds = fmap (replaceInds) $ tIndices tenten
-- bug? --         replaceInds i = if Map.lookup i allIndices == Nothing then i else QIndex (Maybe.fromJust $ Map.lookup i allIndices) (iSpace i) (isDummy i)
-- bug? --         
-- bug? --     ten = fmap (changeInds) (tTensor t)

-- Returns all the possible configurations for a given term
allRenames :: QTerm -> QTerms
allRenames term = zipWith (replaceTerms) allPattern dummyTerms
  where
    patterns = [x ++ y ++ z ++ v| x <- cAll, y <- aAll, z <- vAll, v <- gAll]
    allPattern = fmap (Map.fromList) patterns
    cAll = cAllMap $ getSummedBody term
    aAll = aAllMap $ getSummedBody term
    vAll = vAllMap $ getSummedBody term
    gAll = gAllMap $ getSummedBody term
    dummyTerms = take (length allPattern) $ repeat term

    replaceTerms pattern t = QTerm (tNum t) (tCoeff t) ten
      where
        changeInds tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
          where
            inds = fmap (replaceInds) $ tIndices tenten
            replaceInds i = if Map.lookup i pattern == Nothing then i else QIndex (Maybe.fromJust $ Map.lookup i pattern) (iSpace i) (isDummy i)
  
        ten = fmap (changeInds) (tTensor t)

-- Returns all the possibly ordered tensors
rotateTensors :: QTerm -> QTerms
rotateTensors term = fmap makeTerms tensors
  where
    tensors = List.permutations $ tTensor term
    makeTerms t = QTerm (tNum term) (tCoeff term) t

-- Messy! -- -- Rotate all the indices of each tensor in a given term
-- Messy! -- rotateAllIndices :: QTerm -> QTerms
-- Messy! -- rotateAllIndices term = fmap makeTerms $ concat . head . foldl allPossible [] . tensors
-- Messy! --   where
-- Messy! --     tensors = fmap tgetConfs $ tTensor term
-- Messy! --     allPossible (xs:ys:xss) = (++) <$> xs <*> ys : xss
-- Messy! --     makeTerms t = QTerm (tNum term) (tCoeff term) t

-- Rotate all the indices of each tensor in a given term
rotateAllIndices :: QTerm -> QTerms
rotateAllIndices term = fmap makeTerms $ Utils.makeCombi tensors
  where
    tensors = fmap tgetConfs $ tTensor term
    makeTerms t = QTerm (tNum term) (tCoeff term) t

-- Generate all the possible configurations of a given term
generateAllConfs :: QTerm -> QTerms
generateAllConfs term = concat . fmap (rotateAllIndices) $ concat . fmap (rotateTensors) $ allRenames term

-- -- Kill Kronecker's delta
-- killKDeltas :: QTerm -> Maybe QTerm
-- killKDeltas term = Just (QTerm (tNum term) (tCoeff term) exceptkDeltas) 
-- --killKDeltas term = if oldSpace == newSpace then Just (QTerm (tNum term) (tCoeff term) exceptkDeltas) else Nothing
--   where
--     -- First, construct a Map [(killed Index, killer Index)]
--     tensors = tTensor term
--     kdList  = [x | x <- tensors, tLabel x == "kD"]
--     mkMap   = Map.fromList $ fmap killSchedule kdList
--     killSchedule kDelta
--       |     (isDummy ind1) &&     (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
--       |     (isDummy ind1) &&     (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind1)
--       |     (isDummy ind1) && not (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
--       | not (isDummy ind1) &&     (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind1)                                                                               
--       |     (isDummy ind1) && not (isDummy ind2)                             = (ind1, ind2)
--       | not (isDummy ind1) &&     (isDummy ind2)                             = (ind2, ind1)
--       |     (isDummy ind1) &&     (isDummy ind2)                             = (ind1, ind2)
--       | not (isDummy ind1) && not (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
--       | not (isDummy ind1) && not (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind2)
--       | not (isDummy ind1) && not (isDummy ind2)                             = (ind2, ind2)
--       | otherwise = error "Algorithmic error"
--       where
--         ind1 = (tIndices kDelta) !! 0
--         ind2 = (tIndices kDelta) !! 1
--     -- Then, construct new tensors 
--     newTensors = fmap (Utils.repeatN replaceTensor (length kdList)) tensors
--     replaceTensor tenten = QTensor (tLabel tenten) inds (tSymm tenten) (tComm tenten)
--       where
--         inds = fmap (replaceInds) $ tIndices tenten
--         replaceInds i = if Map.lookup i mkMap == Nothing then i  else Maybe.fromJust $ Map.lookup i mkMap
--     -- If index spaces are same return newTensor else all are Nothing
--     oldSpace = fmap (iSpace) $ concat $ map (tIndices) tensors
--     newSpace = fmap (iSpace) $ concat $ map (tIndices) newTensors
--     --exceptkDeltas = [x | x <- newTensors, not (tLabel x == "kD" && (tIndices x) !! 0 == (tIndices x) !! 1)]
--     exceptkDeltas = [x | x <- newTensors]    

-- Kill Kronecker's delta (Still buggy??)
killKDeltas :: QTerm -> Maybe QTerm
--killKDeltas term = Just (QTerm (tNum term) (tCoeff term) exceptkDeltas) 
--killKDeltas term = if oldSpace == newSpace then Just (QTerm (tNum term) (tCoeff term) exceptkDeltas) else Nothing
killKDeltas term = if foldr (||) False $ fmap isZero kdList then Nothing else Just (QTerm (tNum term) (tCoeff term) exceptkDeltas)
  where
    -- First, construct a Map [(killed Index, killer Index)]
    tensors = tTensor term
    kdList  = [x | x <- tensors, tLabel x == kDName_]

    isZero :: QTensor -> Bool
    isZero kd =
      let
        ind1 = (tIndices kd) !! 0
        ind2 = (tIndices kd) !! 1
        space1 = iSpace ind1
        space2 = iSpace ind2
      in space1 /= space2 && space1 /= Generic && space2 /= Generic

    mkMap   = Map.fromList $ fmap killSchedule kdList
    killSchedule kDelta
      |     (isDummy ind1) &&     (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
      |     (isDummy ind1) &&     (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind1)
      |     (isDummy ind1) && not (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
      | not (isDummy ind1) &&     (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind1)                                                                               
      |     (isDummy ind1) && not (isDummy ind2)                             = (ind1, ind2)
      | not (isDummy ind1) &&     (isDummy ind2)                             = (ind2, ind1)
      |     (isDummy ind1) &&     (isDummy ind2)                             = (ind1, ind2)
      | not (isDummy ind1) && not (isDummy ind2) && (iSpace ind1 == Generic) = (ind1, ind2)
      | not (isDummy ind1) && not (isDummy ind2) && (iSpace ind2 == Generic) = (ind2, ind2)
      | not (isDummy ind1) && not (isDummy ind2)                             = (ind2, ind2)
      | otherwise = error "Algorithmic error"
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
    --exceptkDeltas = [x | x <- newTensors]    

-- Function to return the combined terms. This can be implemented by using State monad more elegantly?
combineTerms :: QTerms -> QTerms
combineTerms terms = fmap (changeFactors) origTerms
  where
    maxTerms = fmap (maximum . generateAllConfs) terms
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

---------------------------------------------------------------------------------------------
-- Small utilities
---------------------------------------------------------------------------------------------
    
-- fucntion returns all the renaming combinations for core indices
cAllMap :: QIndices -> [[(QIndex, String)]]
cAllMap [] = []
cAllMap inds = zipWith (zipWith (\x y -> (x,y))) allFst allRenames
  where 
    zipFst     = fmap fst $ cMap inds 0
    allRenames = List.permutations $ fmap snd $ cMap inds 0
    allFst     = take (length allRenames) $ repeat zipFst 

-- Fucntion returns all the renaming combinations for active indices
aAllMap :: QIndices -> [[(QIndex, String)]]
aAllMap [] = []
aAllMap inds = zipWith (zipWith (\x y -> (x,y))) allFst allRenames
  where 
    zipFst     = fmap fst $ aMap inds 0
    allRenames = List.permutations $ fmap snd $ aMap inds 0
    allFst     = take (length allRenames) $ repeat zipFst 

-- Fucntion returns all the renaming combinations for virtual indices
vAllMap :: QIndices -> [[(QIndex, String)]]
vAllMap [] = []
vAllMap inds = zipWith (zipWith (\x y -> (x,y))) allFst allRenames
  where 
    zipFst     = fmap fst $ vMap inds 0
    allRenames = List.permutations $ fmap snd $ vMap inds 0
    allFst     = take (length allRenames) $ repeat zipFst 

-- Fucntion returns all the renaming combinations for generic indices
gAllMap :: QIndices -> [[(QIndex, String)]]
gAllMap [] = []
gAllMap inds = zipWith (zipWith (\x y -> (x,y))) allFst allRenames
  where 
    zipFst     = fmap fst $ gMap inds 0
    allRenames = List.permutations $ fmap snd $ gMap inds 0
    allFst     = take (length allRenames) $ repeat zipFst

---------------------------------------------------------------------------------------------    
---------------------------------------------------------------------------------------------
-- Function to extract dummy core indices
cMap :: QIndices -> Int -> [(QIndex, String)]
cMap [] _ = []
cMap (i:is) num
  | (iSpace i) == Core && (isDummy i) = [(i, "c" ++ show num)] ++ (cMap is (num+1))
  | otherwise                         = cMap is num

-- Function to extract dummy active indices
aMap :: QIndices -> Int -> [(QIndex, String)]
aMap [] _ = []
aMap (i:is) num
  | (iSpace i) == Active && (isDummy i) = [(i, "a" ++ show num)] ++ (aMap is (num+1))
  | otherwise                           = aMap is num

-- Function to extract dummy virtual indices
vMap :: QIndices -> Int -> [(QIndex, String)]
vMap [] _ = []
vMap (i:is) num
  | (iSpace i) == Virtual && (isDummy i) = [(i, "v" ++ show num)] ++ (vMap is (num+1))
  | otherwise                            = vMap is num

-- Functon to extract dummy generic indices
gMap :: QIndices -> Int -> [(QIndex, String)]
gMap [] _ = []
gMap (i:is) num
  | (iSpace i) == Generic && (isDummy i) = [(i, "g" ++ show num)] ++ (gMap is (num+1))
  | otherwise                            = gMap is num

---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
-- Normal ordering function for creation and annihilation operators
normalOrderA :: QTerm -> QTerms
normalOrderA term =
  if length cres == length des then zipWith (makeTerm) signs kDeltas
  else error "normalOrderOp: Numbers of creation and annihilation operators should be equal for the current implementation"
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
      | iSpace ((tIndices op1) !! 0) == Active  || iSpace ((tIndices op2) !! 0) == Active  = error "normalOrderOp: Can't handle index for the active MOs"
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
normalOrderE term = iterateNormalOrder restEops [initTerm]
  where
    tensors   = tTensor term
    operators = [x | x <- tensors, take (length sfGenName_) (tLabel x) == sfGenName_]
    others    = [x | x <- tensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
    otherOps  = [x | x <- others, (take (length sdGenName_) (tLabel x) == sdGenName_) || tLabel x == creName_ || tLabel x == desName_]
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

    -- Contract given two generators -- It's working!!
    contractGen :: QTerm -> QTerms
    contractGen thisTerm
      | (length incommutables) >  2 = error "contractGen: Length of operators given should be lower than 2."
      | (length incommutables) == 2 = concat . fmap (fmap makeTerm) $ fmap (makeContractions) [0..(min order1 order2)]
      | otherwise                   = [thisTerm]
      where
        thisTensors = tTensor thisTerm
        commutables   = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) /= sfGenName_]
        incommutables = [x | x <- thisTensors, take (length sfGenName_) (tLabel x) == sfGenName_]
        e1 = incommutables !! (length incommutables-2)
        e2 = incommutables !! (length incommutables-1)
        order1     = floor $ (fromIntegral $ length $ tIndices e1) / (fromIntegral 2)
        order2     = floor $ (fromIntegral $ length $ tIndices e2) / (fromIntegral 2)
        upperInds1 = take (order1) (tIndices e1)
        upperInds2 = take (order2) (tIndices e2)
        lowerInds1 = reverse $ take (order1) $ reverse (tIndices e1)
        lowerInds2 = reverse $ take (order2) $ reverse (tIndices e2) 

        -- Lower indices of e1 and upper indices of e2 are contracted
        contPairs = [(x, y) | x <- lowerInds1, y <- upperInds2]
        e1PairInds = Map.fromList $ zip lowerInds1 upperInds1  -- -> Lower indices are keys -> Upper indices
        e2PairInds = Map.fromList $ zip upperInds2 lowerInds2  -- -> Upper indices are keys -> Lower indices

        e1Pairs = zip upperInds1 lowerInds1  
        e2Pairs = zip upperInds2 lowerInds2  

        e1Maps = Map.fromList $ zip upperInds1 e1Pairs  -- -> Upper indices are keys -> (Upper, Lower)
        e2Maps = Map.fromList $ zip lowerInds2 e2Pairs  -- -> Lower indices are keys -> (Upper, Lower)

        -- Returns terms with correct prefactors and tensors
        makeTerm :: QTensors -> QTerm
        makeTerm tenten = baseTerm (tNum term) (tCoeff term) (commutables ++ tenten)

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
