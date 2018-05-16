Polynomials Library
===================

*This is the OpenMP implementation branch*

This project is a supporting library for [TreElM](https://bitbucket.org/apesteam/treelm).
It does not work on its own, but rather needs to be included in
other projects, which also include TreElM.
The library provides some functionality to deal with Legendre polynomials,
especially transforming between modal and nodal representations.

We are now including FXTPACK, and are working on integrating it to allow a
fast polynomial transformation directly from Legendre modes to Legendre nodes.

License
=======

This library is licensed under the terms of the 2-clause BSD license reproduced below.
This means that polynomials is free software and can be used, reproduced, modified,
distributed and redistributed also for commercial purposes under the conditions
of the BSD license.
The only requirement is that some credit to the authors is given by putting this
copyright notice somewhere in your project.


---
Copyright (C) 2015 University of Siegen.
All rights reserved.
 
Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
 
1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY UNIVERSITY OF SIEGN “AS IS” AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of University of Siegen.
