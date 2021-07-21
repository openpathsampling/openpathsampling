The `regressions` directory includes tests that are more complex and either:

1. involve too many parts to be clear exactly which unit test module they belong to;
2. would add too much complication to the unit tests in a module, and therefore are better kept separate.


In general, these will be tests that were discovered by users exposing bugs. Test modules in here may only contain a single (or very few) tests, and will often be in response to a single issue.
