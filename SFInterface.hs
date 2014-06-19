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

-----------------------------------------------------
-- Interface that provides spin-free operators
-----------------------------------------------------
module SFInterface
(
  module BraHKet,
  module Interface.SFOperator,
  module Interface.SFIndex         
) where

import BraHKet
import Interface.SFOperator
import Interface.SFIndex
-----------------------------------------------------
