# WEST test-suite

## Executing Tests

In order to manually run the test, simply run `make all` in the test-suite.

## Adding New Tests
Additional tests can be added with the following procedure:

1) Create a new test directory with:
  - *prepare_inputs.sh*, to generate all inputs (including links to download pseudopotential files)
  - *ref* subdirectory, that contains reference data
  - *Makefile*, that contains instructions about the execution of the codes 

We recommend that you copy and adjust an existing test script.

2) Add the new subdirectory name to the list in the *test-suite/Makefile*.

3) Customize the file *main_test.py*.
