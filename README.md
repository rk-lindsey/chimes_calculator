The Chebyshev Interaction Model for Efficient Simulation (ChIMES) is a machine-learned interatomic potential targeting chemistry in condensed phase systems. ChIMES models are able to approach quantum-accuracy through a systematically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. Though originally developed to enable description of organic molecular materials, ChIMES has successfully been applied to systems spanning ambient water to molten carbon, and leveraged as correction for density functional based tight binding simulations.

The ChIMES calculator comprises a flexible tool set for evaluating ChIMES interactions (e.g. in simulations, single point calculations, etc). Users have the option of directly embedding the ChIMES calculator within their codes (e.g. see ‘’The ChIMES Calculator,’’ in the documentation for advanced users), or evaluating interactions through the beginner-friendly serial interface, each of which have Python, C++, C, and FORTRAN API’s.

Documentation
----------------

[**Full documentation**](https://chimes-calculator.readthedocs.io/en/latest/) is available.

Community
------------------------

Questions, discussion, and contributions (e.g. bug fixes, documentation, and extensions) are welcome. 

Additional Resources: [ChIMES Google group](https://groups.google.com/g/chimes_software).

Contributing
------------------------

Contributions to the ChIMES calculator should be made through a pull request, with ``develop`` as the destination branch. A test suite log file should be attached to the PR. For additional contributing guidelines, see the [documentation](https://chimes-calculator.readthedocs.io/en/latest/contributing.html).

The ChIMES calculator `develop` branch has the latest contributions. Pull requests should target `develop`, and users who want the latest package versions,
features, etc. can use `develop`.

Releases
--------

For most users, we recommend using the ChIMES calculator [stable releases](https://github.com/rk-lindsey/chimes_calculator/releases).

Each ChIMES calculator release series also has a corresponding branch, e.g. `releases/v0.14` has `0.14.x` versions, and `releases/v0.13` has `0.13.x` versions. We back-port important bug fixes to these branches but we do not advance the package versions or make changes that would otherwise change the way ChIMES calculator is used or behaves. So, you can base your ChIMES deployment on a release branch and `git pull` to get fixes, without the continuous changes that comes with `develop`.  The latest release is always available with the `releases/latest` tag.

Authors
----------------

The ChIMES calculator was developed by Rebecca K. Lindsey, Nir Goldman, and Laurence E Fried.

Contributors can be found [here](https://github.com/rk-lindsey/chimes_calculator/graphs/contributors).


Citing
----------------

See [the documentation](https://chimes-calculator.readthedocs.io/en/latest/citing.html) for guidance on referencing ChIMES and the ChIMES calculator in a publication.


License
----------------

The ChIMES calculator is distributed under terms of [GPL-2.0 License](https://github.com/rk-lindsey/chimes_calculator/blob/main/LICENSE).

LLNL-CODE- 817533