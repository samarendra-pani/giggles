.. _developing:

Developing
==========

Documentation is inspired from `WhatsHap <https://whatshap.readthedocs.io/en/latest/>`_

Developer's Installation
------------------------

For installation giggles for the purpose of development,
please use a conda environment.

Conda environment can be used using these commands::

    git clone git@github.com:samarendra-pani/giggles.git
    cd giggles
    conda create -n giggles-dev python=3.10
    conda activate giggles-dev
    pip install -e .[dev]

This will install all dependencies related to testing, documentation building, and formatting.

We support testing using `pytest <https://docs.pytest.org/en/stable/>`_ and `tox <https://tox.readthedocs.io/>`_, and formatting using `ruff <https://docs.astral.sh/ruff/>`_.

Adding a new subcommand
-----------------------

For creating a new subcommand under giggles, add a new script under ``giggles/cli/``.
Make sure to follow the same format. There are three functions that need to be added: ``add_arguments()``, ``validate()``, and ``main()``.
Refer to ``giggles/cli/cpptests.py`` which has an extreme barebones structure and refer to ``giggles/cli/genotype.py`` for a descriptive structure.

All new subcommands and associated functions should have corresponding test cases written under ``tests/`` for testing with ``pytest``.

The core C++ classes have a different way of testing which will be detailed under ....

Executing Test Cases
--------------------

To execute all the pytst tests (collected in ``tests``) run::

    pytest

To execute all the core C++ tests run::

    giggles cpptests

Whenever you change any Cython code (``.pyx`` files) or C++ code, you need to re-run the
``pip install -e .`` step in order to compile it.

To generate an html report on the test coverage::

    coverage run -m pytest
    coverage report -m
    coverage html
    firefox  htmlcov/index.html

Code Formatting
---------------

To format using ``ruff``, run::

    ruff format giggles/ test/ setup.py

To check for syntax and style errors using ``ruff``, run::

    ruff check giggles/ test/ setup.py

To check for proper documentation build, proper syntax and styling, and pytest on multiple python versions (installed locally), run::

    tox

Before creating a pull request

#. write test cases for the new code
#. run ``tox`` after fixing formatting and syntax with ``ruff``
#. create the pull request, if the tox check passes with your local python version


Using Pre-Commit
----------------

The tox pipeline has been integrated into GitHub as a CI/CD pipeline.
Instead of doing all the checks locally, the checks can be done online.

There is also `pre-commit <https://pre-commit.com/>`__ support which also allows developers to skip the
formatting and syntax checking using ``ruff``.

To install ``pre-commit``, run::

    pip install pre-commit
    cd <directory of git repository>
    pre-commit install

Note: ``pre-commit`` is not automatically installed in the dev environment. So the above installation has to be performed.

Now you ``pre-commit`` will run everytime you commit to the repository. But the pre-commit
run is restricted to the staged files. As an optional step, you can run ``pre-commit`` on
all the file using::

    pre-commit run --all-files

When the staged files fail the conditions defined in the pre-commit conditions, ``ruff`` formatting
will attempt to fix it. After adding the files updated by ``ruff``, attempt to commit again. If it fails,
then the errors have to be manually fixed.

Writing Documentation
---------------------

The documentation for giggles is written in
`reStructuredText format <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
and is translated by `Sphinx <http://www.sphinx-doc.org/>`_ into HTML format.
The documentation is found under ``docs/``

The documentation is hosted on `Read the Docs <https://readthedocs.org/>`_.

For testing the documentation, run developers installation and under ``docs/``, run ``make html``. This creates the htmls for the
files under ``docs/_build/``.


Wrapping C++ classes
--------------------

Giggles' core genotyping algorithm along with `WhatsHap's <https://github.com/whatshap/whatshap/>`_ 
modified phasing algorithm is written in C++, as are many of the core
data structures such as the “Read” class under ``src/``. Along with that, some external
tools written in various languages (but with C++ bindings) are available under ``external\`` 
(Right now only `WFA2-lib <https://github.com/smarco/WFA2-lib/>`_ is required).

Updating a Giggles C++ Code
+++++++++++++++++++++++++++

Let us look at the “Read” class. The following places in the code may need to
be changed if the Read class is changed or extended:

* ``src/read.cpp``: Implementation of the class (C++).

* ``src/read.h``: Header with the class declaration (also normal C++).

* ``giggles/cpp.pxd``: Cython declarations of the class. This repeats – using
  the Cython syntax this time – a subset of the information from the
  ``src/read.h`` file. These are the C++ functions that need to be accessed 
  inside Python. This duplication is required because Cython
  cannot read ``.h`` files (it would need a full C++ parser for that).

* ``giggles/core.pxd``: This contains declarations of all *Cython* classes
  wrapping C++ classes. Note that the class ``Read`` in this file has the
  same name as the C++ class, but that it is not the same as the C++ one!
  The distinction is made by prefixing the C++ class with ``cpp.``, which is
  the name of the module in which it is declared in (that is, the C++ class
  ``Read`` is declared in ``cpp.pxd``). The wrapping (Cython) class ``Read``
  stores the C++ class in an attribute named ``thisptr``. If you add a new
  class, it needs to be added to this file. If you only modify an existing one,
  you probably do not need to change this file.

* ``giggles/core.pyx``: The Cython implementation of the wrapper classes.
  Again, the name ``Read`` by itself is the Python wrapper class and
  ``cpp.Read`` is the name for the C++ class.


Adding a new Giggles C++ Code
+++++++++++++++++++++++++++++

Let us add a new class to the Giggle core C++ code. Let us add the “GFA” class as a
C++ code which needs to be accessed in Python. The following files need to be created
or updated:

* ``src/gfa.cpp``: Implementation of the class (C++).

* ``src/gfa.h``: Header with the class declaration (also normal C++).

* ``giggles/cpp.pxd``: Cython declarations of the class. This repeats – using
  the Cython syntax this time – a subset of the information from the
  ``src/gfa.h`` file. These are the C++ functions that need to be accessed 
  inside Python. This duplication is required because Cython
  cannot read ``.h`` files (it would need a full C++ parser for that).

* ``giggles/core.pxd``: This contains declarations of all *Cython* classes
  wrapping C++ classes. Here we will need to define the Pythonic side of
  the GFA class using ``cdef class GFA:`` and under it define a pointer to
  the C++ object as ``cdef cpp.GFA *thisptr``.

* ``giggles/core.pyi``: The header for the Cython implementation of 
  the wrapper classes. Here we define the header for Python class ``GFA``
  and all its member functions. These are the functions that Python can
  access.

* ``giggles/core.pyx``: The Cython implementation of the ``GFA`` class will
  be written here. All the functions defined in ``giggles/core.pyi`` for ``GFA``
  will be expanded. This file defines the bridge between C++ and Python. It can
  access the C++ side through ``cpp.GFA *thisptr``.

* ``setup.py``: Now that ``src/gfa.h`` and ``src/gfa.cpp`` have been added, they
  need to be considered for compilation. Hence ``src/gfa.cpp`` will need to be added
  to setup.py under the list ``core_cpp_sources``.


Adding test cases for core C++
------------------------------

If you are creating new C++ classes/functions (including the one's not used in Python), then consider 
writing test cases for them as detailed below (if the function has a Python side, you can integrate them
into ``pytest``). Let us again consider the ``GFA`` class in ``src/gfa.cpp`` for testing.

* ``src/tests/test_gfa.cpp``: Implementation of the class (C++).

* ``src/tests/test_gfa.h``: Header with the class declaration (also normal C++).

* ``src/tests.cpp`` and ``src/tests.h``: Import ``src/tests/test_gfa.h`` and call its main function.

* ``setup.py``: Now that ``src/tests/test_gfa.h`` and ``src/tests/test_gfa.cpp`` have been added, they
  need to be considered for compilation. Hence ``src/tests/test_gfa.cpp`` will need to be added
  to setup.py under the list ``test_sources``.

Note: ``src/tests_data.cpp`` contains the function ``assert_msg(bool, string, string)`` which is used as
follows:

* bool: the condition you are checking.
* string: the first string refers to the class you are evaluating (in this example ``GFA``).
* string: the message describing the condition being checked.

The output is ``[GFA] PASSED: <Condition Message>`` if it passes and ``[GFA] FAILED: <Condition Message>``
if it fails. Failure at any point terminates the testing.

Making a Release
----------------

#. Update ``CHANGES.rst``: Set the correct version number and ensure that
   all nontrivial, user-visible changes are listed.

#. Ensure you have no uncommitted changes in the working copy.

#. Run ``tox``.

#. Tag the current commit with the version number (there must be a ``v`` prefix)::

       git tag -a -m "Version 0.1" v0.1

#. Push the tag::

       git push --tags

#. Wait for the GitHub Action to finish. It will deploy the sdist and wheels to
   PyPI if everything worked correctly.

If something went wrong, fix the problem and follow the above instructions again,
but with an incremented revision in the version number. That is, go from version
x.y to x.y.1. PyPI will not allow you to change a version that has already been
uploaded.