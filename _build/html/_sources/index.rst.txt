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

Program will not work if either of those lines is missing.

It requires Python 2.7 or higher.
 

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   code
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
