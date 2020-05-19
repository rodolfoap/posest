    **************************************************************
                   Matlab MEX interface to posest
                    Manolis Lourakis, 2015
    **************************************************************

==================== GENERAL ====================
This directory contains a MEX-file interface for posest, allowing
it to be called directly from matlab.

Tested under matlab version 6.5 (R13) under linux and
version 7.4 (R2007) under Windows


==================== COMPILING ====================
 - On a Linux/Unix system, typing "make" will build the MEX object.

 - Under Windows, use the provided Makefile.w32 as a basis for creating
   your own makefile.


==================== TESTING ====================
At the command line, type
matlab < posest_demo.m
