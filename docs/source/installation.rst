.. _installation:

##################
Ha3Py installation
##################

`Ha3Py` requires at least Python 3.9. All the required dependencies
will be downloaded and installed during the setup process.

Installing the latest release
#############################

The latest release of Ha3Py is available on the `Python Package
Index <https://pypi.org/project/ha3py/>`__.

You can install it easily through ``pip``:

::

   pip install ha3py

To upgrade from a previously installed version:

::

   pip install --upgrade ha3py

Installing a developer package
##############################

If you want to modify the source code, you should clone the project using ``git``:
::

    git clone https://github.com/JanWiszniowski/ha3py.git

Next, go into the ``ha3py`` main directory and install the code in
"editable mode" by running:

::

    pip install -e .

You can keep your local Ha3Py repository updated by running ``git pull`` commands at regular intervals.
Thanks to ``pip`` editable mode, you don't need to reinstall Ha3Py after each update.
