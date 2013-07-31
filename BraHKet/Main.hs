
import BraHKet.Core

main = do
  let i = QIndex "I" Active  False
      k = QIndex "K" Active  False
      a = QIndex "A" Virtual True
      c = QIndex "C" Virtual True
      b = QIndex "B" Virtual False
      d = QIndex "D" Virtual False           
      v = QTensor "V" [i,k,a,c] [[1,2,3,4]] True
      t = QTensor "V" [a,c,b,d] [[1,2,3,4]] True
      g = QTerm 1.0 ["Ecas"] [v, t]
  print $ "Indices : " ++ (show [i,k,a,c])
  print $ "Tensor  : " ++ (show t)
  print $ "Term    : " ++ (show g)
  
