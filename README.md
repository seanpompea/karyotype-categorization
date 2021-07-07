karyoclassifier
===============

This pipeline gets data from `dm_mpn..karyotypes`. For each karyotype
string, it processes it for seven risk schemas:

* CALGB (Cancer and Leukemia Group B) Risk Classification
* AML MRC 1998 Risk Classification
* AML MRC 2010 Risk Classification
* AML SWOG (Southwest Oncology Group) Risk Classification   
* MDS 5-Group Risk Classification 
* MDS IPSS Risk Classification 
* MF DIPSS+ Risk Classification 

It loads the resulting dataset into `dm_mpn..karyotype_classifications`.

## About ISCN notation

A collection of notes on ISCN karyotype string notation can be found
in the include file [knotes.txt](knotes.txt).

## Setup and running

Uses Python 2.7.

:alert: The current approach does _not_ use an encrypted database connection. To
ensure an encrypted connection, you'd need to change the database driver/library
to something else; some more details about this topic are here in this
Nexus article:

* https://nexus.weill.cornell.edu/display/ARCH/Best+Practices+for+Establishing+Secure+Database+Connections+Using+Python

Note that the above article focuses on Python 3 (this app is written in
Python 2). 

Ok, moving on to the actual steps.

Install freetds (or equivalent), which pymssql (which gets installed by the
sqlsrvwrapper) uses. There might ben an alternate approach in terms of 
what to install here but this worked for me:
  
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

    enclave/p06-info.json

That file should have the structure/content as documented here:

    https://github.com/wcmc-research-informatics/sqlsrvwrapper/blob/master/README.md

Create the destination database table; see `create-table.sql` for the DDL.

Start a repl:

    python

Then, to run the pipeline, run these commands:

    >>> import bulk as b
    >>> b.trunctable('p06', 'dm_mpn..karyotype_classifications')
    >>> b.db2db('p06', 'p06', 'dm_mpn..karyotype_classifications')

