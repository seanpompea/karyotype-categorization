Cytogenetic karyotype nomenclature parsing and clinical risk categorization
===============

Certain cytogenetic abnormalities are known to play a role in myeloproliferative neoplasm diseases. Automated parsing and categorization of cytogenetic nomenclature datasets can aid in research, phenotype development, and subsequent clinical decision-making.

This pipeline parses a data set of cytogenetic karyotype strings (in ISCN notation) to identify abnormalities and subsequently apply clinical risk categories using the following categorizaton schemes:

* CALGB (Cancer and Leukemia Group B) Risk Classification
* AML MRC 1998 Risk Classification
* AML MRC 2010 Risk Classification
* AML SWOG (Southwest Oncology Group) Risk Classification   
* MDS 5-Group Risk Classification 
* MDS IPSS Risk Classification 
* MF DIPSS+ Risk Classification

## About ISCN notation 

A collection of notes on ISCN karyotype string notation can be found
in the file [knotes.txt](knotes.txt).

## Setup and running using a database

Uses Python 2.7.

Install freetds (or equivalent), which pymssql (which gets installed by the
sqlsrvwrapper lib) uses.
  
    brew install freetds

Create a virtual env and then activate it:

    mkdir venv
    virtualenv -p python2.7 venv
    source venv/bin/activate

Then install dependencies:

    pip install git+https://github.com/wcmc-research-informatics/kickshaws
    pip install cython
    pip install https://github.com/wcmc-research-informatics/sqlsrvwrapper

Create the database config file at this location:

    enclave/mydb-info.json

That file should have the structure/content as documented here:

    https://github.com/wcmc-research-informatics/sqlsrvwrapper/blob/master/README.md

Create the destination database table; see `create-table.sql` for the DDL.

Start a repl:

    python

Then, to run the pipeline, run these commands:

    >>> import bulk as b
    >>> b.trunctable('mydb', 'dm_mpn..karyotype_classifications')
    >>> b.db2db('mydb', 'mydb', 'dm_mpn..karyotype_classifications')

## Setup and running using CSV files

There are additional functions in `bulk.py` for working with CSVs; take a look at the
code for more details.
