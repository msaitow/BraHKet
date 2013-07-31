
import BrahKet.Core

testTensor =
  let i = QIndex "I" Active  False
      k = QIndex "K" Active  False
      a = QIndex "A" Virtual True
      c = QIndex "C" Virtual True
      b = QIndex "B" Virtual False
      d = QIndex "D" Virtual False
      v = baseERI         [i,k,a,c]
      t = baseTensor "T2" [a,c,b,d] [[0,1,2,3],[1,0,3,2]] True
      g = QTerm 1.0 ["Ecas"] [v, t]
  in [v, t, g]
     
