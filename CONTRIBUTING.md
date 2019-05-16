Contributing to QIIME 2
-----------------------

QIIME 2 is an open-source project, and we are very interested in contributions from the community. Please get in touch through the [QIIME 2 Forum](https://forum.qiime2.org) if you'd like to get involved.

## How to Contribute in Issues

For any issue, there are fundamentally three ways an individual can
contribute:

1. By opening the issue for discussion: if you believe that you
   have uncovered a bug in qiime2, you can create a new issue in the `qiime2/qiime2`
   issue tracker to report it.
2. By helping to triage the issue: This can be done either by providing
   supporting details (a test case that demonstrates a bug), or providing
   suggestions on how to address the issue.
3. By helping to resolve the issue: Typically this is done either in the form
   of demonstrating that the issue reported is not a problem after all, or more
   often, by opening a Pull Request that changes some bit of something in
   `qiime2/qiime2` in a concrete and reviewable manner.

## Triaging a Bug
Contributors are encouraged to help one another make forward progress as much
as possible, empowering one another to solve issues collaboratively. If you
choose to comment on an issue that you feel either is not a problem that needs
to be fixed, or if you encounter information in an issue that you feel is
incorrect, explain *why* you feel that way with additional supporting context,
and be willing to be convinced that you may be wrong. By doing so, we can often
reach the correct outcome much faster.

## Resolving a Bug 

In the vast majority of cases, issues are resolved by opening a Pull Request.
The process for opening and reviewing a Pull Request is similar to that of
opening and triaging issues, but carries with it a necessary review and approval
workflow that ensures that the proposed changes meet the minimal quality and
functional guidelines of the qiime2 project.

here’s how to submit a pull request in github:
Fork the repository and clone it locally. Connect your local to the original “upstream” repository by adding it as a remote. Pull in changes from “upstream” often so that you stay up to date so that when you submit your pull request, merge conflicts will be less likely. (See more detailed instructions here.)
Create a branch for your edits.
Reference any relevant issues or supporting documentation in your PR (for example, “Closes #37.”)
Include screenshots of the before and after if your changes include differences in HTML/CSS. Drag and drop the images into the body of your pull request.
Test your changes! Run your changes against any existing tests if they exist and create new ones when needed. Whether tests exist or not, make sure your changes don’t break the existing project.
Contribute in the style of the project to the best of your abilities. This may mean using indents, semi-colons or comments differently than you would in your own repository, but makes it easier for the maintainer to merge, others to understand and maintain in the future.


### QIIME 2 Users

Check out the [User Docs](https://docs.qiime2.org) - there are many tutorials,
walkthroughs, and guides available. If you still need help, please visit us at
the [QIIME 2 Forum](https://forum.qiime2.org/c/user-support).
* [Tutorials](https://docs.qiime2.org/2019.4/tutorials/)
 * [Concepts](https://docs.qiime2.org/2019.4/concepts/)
 * [Technical Support](https://forum.qiime2.org/c/technical-support)
 Just have a question? Please ask it in our [forum](https://forum.qiime2.org/c/user-support).
  * [Basic Command Line Interface](https://forum.qiime2.org/t/tutorial-basic-command-line-interface-for-beginners-to-qiime2/6804)

### QIIME 2 Developers

Check out the [Developer Docs](https://dev.qiime2.org) - there are many introductions for
architecture, testing, plugin, api reference, internal details and guides available. If you still need help, please
visit us at the [QIIME 2 Forum](https://forum.qiime2.org/c/dev-discussion).
Want to discuss with other developer? please join us [forum](https://forum.qiime2.org/c/dev-discussion)

This document is based heavily on the following:
https://github.com/atom/atom/blob/master/CONTRIBUTING.md
https://github.com/nodejs/node/edit/master/doc/guides/contributing/issues.md
https://opensource.guide/how-to-contribute/
