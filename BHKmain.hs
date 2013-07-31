
import Data.List
import BraHKet.Core

main = do
  let i = QIndex "I" Active  False
      k = QIndex "K" Active  False
      a = QIndex "A" Virtual False
      c = QIndex "C" Virtual False
      j = QIndex "J" Active True
      l = QIndex "L" Active True      
      b = QIndex "B" Virtual True
      d = QIndex "D" Virtual True           

      v  = baseERI         [c,b,a,d]
      r  = baseSFRDM       [i,j,k,l]      
      t  = baseTensor "T2" [j,l,b,d] [[0,1,2,3],[1,0,3,2]] True
      kd = baseKD [b,d]

      g = baseTerm 1.0 [] [v, t, r, kd]
  print $ "Indices : " ++ (show [i,k,a,c])
  print $ "Tensor  : " ++ (show v)
  print $ "Sorted Tensor : " ++ (show $ tsortIndices v)
  print $ "Perms of t  : " ++  (show $ tSymm t)
  print $ "Perms of v  : " ++  (show $ tSymm v)
  print $ "Perms of r  : " ++  (show $ tSymm r)
  print $ "Perms of kd : " ++  (show $ tSymm kd)  
  --print $ "All possible indices for v : " ++ (show $ tallPermInds v)
  print $ "All possible configurations for v : " ++ (show $ tgetConfs v)  
  print $ "Term1    : " ++ (show g)
  print $ "Term2    : " ++ (show $ masquerade g)
  print $ "Term3    : " ++ (show $ getSummedBody g)
--  print $ "Term4    : " ++ (show $ vAllMap $ getSummedBody g)
--   print $ "Term4    : " ++ (show $ allRenames g)
--   print $ "Term5    : " ++ (show $ rotateTensors g)
--   print $ "Term6    : " ++ (show $ rotateAllIndices g)
--   print $ "Term6'   : " ++ (show $ fmap tgetConfs $ tTensor g)
--   print $ "Term7    : " ++ (show $ generateAllConfs g)
--   print $ "Term7'   : " ++ (show $ nub $ generateAllConfs g)
--   print $ "Term8    : " ++ (show $ maximum $ generateAllConfs g)
--   print $ "Term8    : " ++ (show $ length $ generateAllConfs g)
--   print $ "Term8'   : " ++ (show $ length $ nub $ generateAllConfs g)
  print $ "Term9    : " ++ (show $ killKDeltas g)  
--  print $ "Term7    : " ++ (show $ fmap tgetConfs $ tTensor g)
