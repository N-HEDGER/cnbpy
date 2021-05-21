# CNBPY

Package for interacting with ABIDE data.

## Pre-requisites:

1. Install git-annex
2. pip install datalad (I work with 0.14.4)
3. Datalad has an awkward interaction with jupyter notebooks, so to resolve this you will also need to install nest_asyncio (1.5.1).

## For development:

1. *Fork* to make your own version.
2. *Git clone* into a local location on your own system.
3. to install to your environment, use *python setup.py develop*
4. Create a branch relevant to what you are working on e.g. *git branch fmriprep*
5. *git checkout fmriprep*
6. Make changes
7. Test things in notebooks stored in /cnbpy/test/develop, with this at the top:

```python
%load_ext autoreload
%autoreload 2
```

This will allow you to test code you have just modified without restarting the kernel. 

8. Commit changes when happy
9. *git push origin fmriprep:fmriprep*
10. Make a pull request.


## How to use:

Functionality so far is shown in /cnbpy/test/develop

The structure is modular:

1. *datalad.py* - Class for interacting with datalad data.
2. *bids.py* - Class for interacting with BIDS data.
3. *scratchpad.py* - defunct stuff that I don't yet want to get rid of.

Everything is well-annotated so far.

So there will also be:

1. *fmriprep.py* - tools for fmriprepping data.
2. *analyze.py* - tools for analysing the data.
3. *visualise.py* - tools for visualising the data.

## Principles

1. Try to get this to work for in principle any BIDS dataset. Write classes suitable for any BIDS data - but then superclasses to sit on top that are specific for ABIDE.
2. Keep things modular and object based.
3. Annotate things in the same style that I have.
4. Try to depend on as few external packages as possible. Try not to deviate from standard scientific python modules.



