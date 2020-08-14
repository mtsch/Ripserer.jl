# Contributor Guide

All contribution are welcome! Please read these guidelines before starting to work on this
project. Following these guidelines will reduce friction and improve the speed at which your
code gets merged.

## Bug Reports

If you notice code that crashes, is incorrect, or is too slow, please file a bug report. The
report should be raised as a GitHub issue with a minimal working example that reproduces the
condition. The example should include any data needed. If the problem is incorrectness, then
please post the correct result along with an incorrect result.

Please include version numbers of all relevant libraries and Julia itself.

## Documentation Fixes and Improvements

If you find typos, weird wording or missing information in the docs, don't hesistate to open
a pull request or issue. When editing docs, I recommend clicking "Edit on GitHub" and using
GitHub's online editor.

Ideas for examples or tutorials are more than welcome as well!

## Code Style

This project attempts to follow the [BlueStyle](https://github.com/invenia/BlueStyle). In
addition to that, I would like to highlight the following:

* Use descriptive variable names. Exceptions can be made for variables with a short scope.
* Avoid long lines. Line length should be at most 92 charactes.
* Avoid introducing new dependencies if you can.
* Don't over-optimize non-critical code, but avoid doing unnecessary work if you can.
* Open an issue before opening complex PR so we can discuss the changes.
