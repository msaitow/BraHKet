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

---------------------------------
-- Handy interface --------------
---------------------------------
module Interface.SFIndex (

-- Core MOs
i ,
j ,
k ,    
l , 
c0,
c1,
c2,
c3,

-- Active MOs
r ,
s ,
p ,
q ,
t , 
u , 
a0,
a1,
a2,
a3,

-- Virtual MOs
a ,
b ,
c ,
d ,
v0,
v1,
v2,
v3,

-- Generic MO
w, 
x,
y,
z,
g0,
g1,
g2,
g3

) where

---import Fafnir.Core
import BraHKet.Core

-- Declaration for molecular orbital indices
i = sfIndex "i" Core    NonDummy      
j = sfIndex "j" Core    NonDummy           
k = sfIndex "k" Core    NonDummy      
l = sfIndex "l" Core    NonDummy

p = sfIndex "p" Active  NonDummy
q = sfIndex "q" Active  NonDummy
r = sfIndex "r" Active  NonDummy      
s = sfIndex "s" Active  NonDummy      
t = sfIndex "t" Active  NonDummy      
u = sfIndex "u" Active  NonDummy

a = sfIndex "a" Virtual NonDummy
b = sfIndex "b" Virtual NonDummy
c = sfIndex "c" Virtual NonDummy
d = sfIndex "d" Virtual NonDummy

w = sfIndex "w" Generic NonDummy
x = sfIndex "x" Generic NonDummy
y = sfIndex "y" Generic NonDummy
z = sfIndex "z" Generic NonDummy

---------------------------------------
---------------------------------------
c0 = sfIndex "c0" Core    Dummy  
c1 = sfIndex "c1" Core    Dummy 
c2 = sfIndex "c2" Core    Dummy 
c3 = sfIndex "c3" Core    Dummy 

a0 = sfIndex "a0" Active  Dummy  
a1 = sfIndex "a1" Active  Dummy 
a2 = sfIndex "a2" Active  Dummy 
a3 = sfIndex "a3" Active  Dummy 

v0 = sfIndex "v0" Virtual Dummy  
v1 = sfIndex "v1" Virtual Dummy 
v2 = sfIndex "v2" Virtual Dummy 
v3 = sfIndex "v3" Virtual Dummy 

g0 = sfIndex "g0" Generic Dummy  
g1 = sfIndex "g1" Generic Dummy 
g2 = sfIndex "g2" Generic Dummy 
g3 = sfIndex "g3" Generic Dummy 
---------------------------------------
---------------------------------------
