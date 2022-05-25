# WEST test-suite

## Executing Tests

In order to manually run the test, simply run `make all` in the test-suite.

## Adding New Tests
Additional tests can be added with the following procedure:

Create a new test directory with:
  - *prepare_inputs.sh*, to generate all inputs (including links to download pseudopotential files)
  - a *ref* subdirectory, that contains reference data
  - python scripts to perform tests

We recommend that you copy and adjust an existing test script.

Add the new test directory to the list in the test-suite/Makefile.

