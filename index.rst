.. rna2D documentation master file, created by
   sphinx-quickstart on Wed Feb 22 16:12:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to rna2D's documentation!
=================================
This module analyze RNA secondary structure in dot-bracket notation, associates groups of nucleotydes with certein structures, builds graph of connections
and writes results to a txt file wich can be easily used as a database. Program works with text files formated as below:

>name 

RNA_SEQUENCE

dot_bracket_notation

RNA sequnce can consist of letters U, A, G, C, u, a, g, c.
Program accepts dot-bracket sequnces with [, ], {, }, >, < as pseudo-knots.

Program will not work if either of those lines is missing.
Output, is stored in "db.txt" file.

Program accepts path to the file as parameter or will ask for a path if parameter is not provided.

It requires Python 2.7 with future module installed or Python 3.5.
 

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   code
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
