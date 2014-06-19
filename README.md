-------------------------------------------------------------------------------------

                                      ,--,
                                    ,--.'|
                          ,--.   ,--,  | :.--,
                         /  /|,---.'|  : '|\  \
                        '  / '|   | : _' |` \  `
                       /  / / :   : |.'  | \ \  \
                      /  / ,  |   ' '  ; :  , \  \                                                                                                                                                                                                                                    
                      \ '\ \  '   |  .'. |  / /` /
                       \  \ ' |   | :  | ' ` /  /
                        \  . |'   : |  : ;| .  /                                                                                                                                                                                                                                      
                         \__\.|   | '  ,/ ./__/
                              ;   : ;--'                                                                                                                                                                                                                                              
                              |   ,/
                              '---'

    BraHKet :: A Symbolic Manipulation Module for Many-Fermionic Operators
                              Masaaki Saitow

-------------------------------------------------------------------------------------

# Overview
  
**BraHKet is an heir to the symbolic manipulation engine in femto library, which is written in the object-oriented C++.** BraHKet heavily uses the abstract purely-functional programming technique to expand, combine and factorize the many-fermionic formulas. The user specified many-body ansatz (e.g., *Configuration Interaction*, *Coupled-Cluster*, or *Canonical Transformation* theory) is expanded according to the Wick's theorem to give the spin-free tensor contractions. 

## Structure of the code:

- ***BraHKet/Core.hs*** - Main source file in which the fundamental data (e.g., index for molecular orbital, tensor object and the tensor contraction term object) are defined. The second quantization operators such as the creation/annihilation operators, spin-free/dependent excitation operators are defined as a special *uncommutable* type of the tensor data. Note that the facilities for handling the spin-dependent excitations are still under development. The main objects defined in this file are as follows:

	*  **BraHKet.Core.QIndex** - QIndex is a data type that represents the molecular orbital index and is composed of four member variables; the index label (**String** data), orbital type (**QSpace** data), spin-variable (**QSpin** data) and a logical (**QSum**) data. **QSpace** is an enumeration data, which can take a value either of *Core*, *Active*, *Virtual* or *Generic*. **QSum** is also an enumeration that takes value either of *Dummy* or *NonDummy*, former of which is responsible for a summation over the orbital index. **QSpin** specifies the information on the spin-state of the associated orbital index and can take a value either of *Alpha*, *Beta* or *SpinFree*.
	
	* **BraHKet.Core.QTensor** - QTensor represents a tensor data such as not only one- and two-electron integrals but also incommutable second quantized operators. QTensor possesses four member variables; tensor name (**String** data), associate orbital indices (**QIndices** data), sequences of integers (**Permut** data), which represents the permutational symmetry of the indices, and a binary (**QNature**) data that represents whether this tensor is a commutative or not. QNature takes a value either of *Classical* or *Operator*. When *Classical* is set, the tensor is specified as a commutative while *Operator* represents an incommutativity of the tensor.
	
	* **BraHKet.Core.QTerm**   - QTerm is an elementary constitutional unit of the tensor contractions and is composed of three member variables; a numerical factor (**Double** data), a literal coefficient (**Coeffs** data) and a tensor product (**QTensors** data) that constitutes a tensor contraction term.

- ***BraHKet/Wick.hs*** - A module in which several functions are implemented to perform the normal ordering of the various second quantized operators ranging from the usual creation/annihilation operators to the spin-free/dependent unitary group generators. Since the spin-dependent generator takes form as a kind of the anti-symmetric tensor, which however cannot be handled by means of the **QTensor** object so far, the spin-dependent features has yet to be fully functional. A direct-algorithm to carry out the normal ordering of the multiple-commutators such that appears in the *Canonical Transformation* or in the *Full-Internally Contracted Multireference Interaction* theory is available. 

- ***BraHKet/Fock.hs*** - The module provides a series of functions that perform the Fock matrix construction by combining the terms with a specific type of the one- and two-body integrals. Invoking these functions, all the terms appear to be *linked* because the *un-linked* are *closed-loop* contribution are contracted to the active and core Fock matrices.
		
- ***BraHKet/Utils.hs*** - Small utility functions necessary for invoking various algebraic operations such as the normal ordering of the operator and like terms combining are provided. Most of the functions in this source is devoted to the combinatoric algorithms that are not provided in the standard modules of GHC.

As a handy interface to BraHKet, a auxiliary module named ***Interface*** is also available, intended to reduce amount of typing to evaluate the second quantization expressions. Some fundamental orbital indices are defined for core, active and virtual molecular orbitals in ***Interface/SFIndex.hs***. On the basis of this, in ***Interface/SFOperator*.hs**, several auxiliary functions are implemented, for example, for constructing the Hamiltonian or for performing the algebraic manipulations in the correct order.  

## How to compile the examples:

Several sample codes are put in ***Codes/***. To compile these codes, the user needs append ***BraHKet/*** or ***Interface/*** to the standard search path of GHC [(See the manual for GHC)](http://www.haskell.org/ghc/docs/latest/html/users_guide/separate-compilation.html#search-path). For example, compilation of the ***Codes/icci_external_external.hs*** from ***Codes/*** directory requires the following options: 

``bash> ghc -O2 -i:"../" ./icci_external_external.hs``

where ``-O2`` option is added for the rapid execution of the symbolic manipulation. In general, if the compilation of a source code requires the linkage to ***BraHKet***, GHC requires ``-i:"PATH_TO_BraHKet"`` option.

## On extending the capability:

The symbolic manipulation and the optimization of the tensor equations often accompany a substantial amount of the combinatorial operations. To achieve the better performance of the code and the usability to avoid the **bug**, this package is written on the basis of the following strategies:

* The external constructor is provided to define the specific type of the tensor object such as **baseERI** for the symmetric electron-repulsion integrals and **baseSFGen** for the spin-free unitary group generator. Since these tensors are distinguished by the tensor names, which can be extracted by using tLabel record, use of the **QTensor** statement to define such a predefined tensor is what should be avoided. When user defines a tensor whose constructor is not provided, **baseTensor** works as the generic constructor. Also for the other tensor objects, use of the external constructor, whose name usually starts with **base**, is highly recommended. The naive constructors such as **QIndex**, **QTensor** and **QTerm** are public only for the reason that call of them from the other modules achieves much better performance.

* On defining a brand new tensor object, call of the combinatorial functions, such as those defined in the **BraHKet.Core.Utils**, from the constructor should be avoided. If these functions exist in the constructor, a drastic amount of operations arise and the performance of the execution is critically affected. As a design strategy, the fundamental functions like the constructors are preferred to be as simple and light as possible. To analyze performance of the code, GHC, for instance, requires the following options:
``bash> ghc -O2 -i:"../" -prof -auto-all ./icci_external_external.hs``
and also on execution of the compiled binary, ``./icci_external_external +RTS -p -RTS``. The option ``-p`` is responsible for the output of the profiled data into a text file. This helps user to specify the bottleneck of the execution.


# Authorship

Any modifications or use of this software should cite the following papers on which the algorithm implemented in this package is based and femto itself:

 - *M. Saitow, Y. Kurashige and T. Yanai, J. Chem. Phys. **139**, 044118 (2013)*
 - *M. Saitow, Femto :: An Integrated Toolset for the Automated Tensor Generation*

Copyright (C) 2013-2014 by Masaaki Saitow (msaitow514@gmail.com)

This program is free software; you can redistribute it and/or modify                                                                                                                                                                                                                  
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or                                                                                                                                                                                                                     
(at your option) any later version

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                                                                                                                                                        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software                                                                                                                                                                                                                           
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

# Contact

For questions concerning the software or the copyright write to the author:

                            Masaaki Saitow
                        (msaitow514@gmail.com)

