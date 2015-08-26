============
Installation
============
At the command line::

    $ easy_install riboplot

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv riboplot
    $ pip install riboplot

.. note:: 

    RNA coverage plot requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed. 
    
    This release of riboplot has been tested with bedtools version ``2.17.0``.

    On Ubuntu and derivatives, bedtools can be installed from the repositories using::

        sudo apt-get install bedtools
