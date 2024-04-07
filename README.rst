=============
tonellishanks
=============

Pure-Python implementation of the `Tonelli-Shanks algorithm <https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm>`__ for calculating a square root modulo a prime.

|pypi|

.. |pypi| image:: https://badge.fury.io/py/tonellishanks.svg#
   :target: https://badge.fury.io/py/tonellishanks
   :alt: PyPI version and link.

Purpose
-------
This library provides a pure-Python implementation of the `Tonelli-Shanks algorithm <https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm>`__. Based on a `simple implementation <https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python>`__, this library prioritizes portability over performance.

Installation and Usage
----------------------
This library is available as a `package on PyPI <https://pypi.org/project/tonellishanks>`__:

.. code-block:: bash

    python -m pip install tonellishanks

The library can be imported in the usual way:

.. code-block:: python

    from tonellishanks import tonellishanks

Development
-----------
All installation and development dependencies are fully specified in ``pyproject.toml``. The ``project.optional-dependencies`` object is used to `specify optional requirements <https://peps.python.org/pep-0621>`__ for various development tasks. This makes it possible to specify additional options (such as ``lint``) when performing installation using `pip <https://pypi.org/project/pip>`__:

.. code-block:: bash

    python -m pip install .[lint]

Testing and Conventions
^^^^^^^^^^^^^^^^^^^^^^^
All unit tests are executed and their coverage is measured when using `pytest <https://docs.pytest.org>`__ (see the ``pyproject.toml`` file for configuration details):

.. code-block:: bash

    python -m pip install .[test]
    python -m pytest

Alternatively, all unit tests are included in the module itself and can be executed using `doctest <https://docs.python.org/3/library/doctest.html>`__:

.. code-block:: bash

    python src/tonellishanks/tonellishanks.py -v

Style conventions are enforced using `Pylint <https://pylint.readthedocs.io>`__:

.. code-block:: bash

    python -m pip install .[lint]
    python -m pylint src/tonellishanks

Contributions
^^^^^^^^^^^^^
In order to contribute to the source code, open an issue or submit a pull request on the `GitHub page <https://github.com/lapets/tonellishanks>`__ for this library.

Versioning
^^^^^^^^^^
The version number format for this library and the changes to the library associated with version number increments conform with `Semantic Versioning 2.0.0 <https://semver.org/#semantic-versioning-200>`__.

Publishing
^^^^^^^^^^
This library can be published as a `package on PyPI <https://pypi.org/project/tonellishanks>`__ by a package maintainer. First, install the dependencies required for packaging and publishing:

.. code-block:: bash

    python -m pip install .[publish]

Ensure that the correct version number appears in ``pyproject.toml``. Create and push a tag for this version (replacing ``?.?.?`` with the version number):

.. code-block:: bash

    git tag ?.?.?
    git push origin ?.?.?

Remove any old build/distribution files. Then, package the source into a distribution archive:

.. code-block:: bash

    rm -rf build dist src/*.egg-info
    python -m build --sdist --wheel .

Finally, upload the package distribution archive to `PyPI <https://pypi.org>`__:

.. code-block:: bash

    python -m twine upload dist/*
