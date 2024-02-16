We welcome contributions to the **surtvep** package. To submit a contribution:

1. [Fork](https://github.com/UM-KevinHe/surtvep/fork) the repository and make your changes.

2. Submit a [pull request](https://help.github.com/articles/using-pull-requests).

3. Ensure that you have signed the contributor license agreement. After you sign, you can click the "Recheck" link in that comment and the check
   will flip to reflect that you've signed.

## How to make changes

Before you submit a pull request, please do the following:

* Add an entry to NEWS.md concisely describing what you changed.

* If appropriate, add unit tests in the tests/ directory.

* Run Build->Check Package in the RStudio IDE, or `devtools::check()`, to make sure your change did not add any messages, warnings, or errors.

Doing these things will make it easier for the **surtvep** development team to evaluate your pull request. Even so, we may still decide to modify your code or even not merge it at all. Factors that may prevent us from merging the pull request include:

* breaking backward compatibility
* is hard to understand
* is hard to maintain in the future
* is computationally expensive
* is not intuitive for people to use

We will try to be responsive and provide feedback in case we decide not to merge your pull request.


## Filing issues

If you find a bug in **surtvep**, you can also [file an issue](https://github.com/UM-KevinHe/surtvep/issues/new). Please provide as much relevant information as you can, and include a minimal reproducible example if possible.
