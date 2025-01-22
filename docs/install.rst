Installation
============

Requirements
------------

EISPAC depends on a number of Python packages that are commonly used in
scientific and solar research. Normally, the installation process will
automatically check and install missing dependencies, assuming your
environment is configured appropriately. If it does not, you may wish to
try installing the required packages individually first.

-  python >= 3.9

-  numpy >= 1.18

-  scipy >= 1.4

-  matplotlib >= 3.1

-  h5py >= 2.9

-  astropy >= 4.2.1

-  sunpy >= 4.0

-  ndcube >= 2.0

-  pyqt >= 5.9

-  parfive >= 1.5

-  python-dateutil >= 2.8

-  tomli >= 1.1.0 (python < v3.11 only)

.. _sec-install:

Standard Installation
---------------------

EISPAC is available on PyPI. To install, use the following command,

::

   >>> python -m pip install eispac

To upgrade the package, please use:

::

   >>> python -m pip install --upgrade eispac

pip should automatically install all package dependencies. If it does not, please
see the list of required packages above. Note: if you are using conda to manage your
Python packages, you may wish to install or update the dependencies manually first,
before installing eispac using pip (this is by no means required, but it can help
simplify updating packages).

Manual Installation
-------------------

1.  Download or clone "eispac" to a convenient location on your computer (it does not matter where).

::

   >>> git clone https://github.com/USNavalResearchLaboratory/eispac.git

2.  Open a terminal and navigate to the directory
3.  To install:

::

   >>> python -m pip install .

4.  To upgrade:

::

   >>> python -m pip install --upgrade .


The package should then be installed to the correct location for your current Python
environment. You can now import the package using `import eispac`.

Now, you should be all set to do some science!
