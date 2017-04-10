chumpy
======

Chumpy is a Python-based framework designed to handle the auto-differentiation problem, 
which is to evalute an expression and its derivatives with respect to its inputs, with 
the use of the chain rule. Chumpy is intended to make construction and local
minimization of objectives easier. Specifically, it provides:

- Easy problem construction by using Numpy’s application interface
- Easy access to derivatives via auto differentiation
- Easy local optimization methods (12 of them: most of which use the derivatives) 

Chumpy comes with its own demos, which can be seen by typing the following::

>> import chumpy
>> chumpy.demo() # prints out a list of possible demos

Licensing is specified in the attached LICENSE.txt.

