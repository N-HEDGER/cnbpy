# CNBPY

Package for interacting with ABIDE data.

## Pre-requisites:

1. Install git-annex
2. pip install datalad (I work with 0.14.4)
3. Datalad has an awkward interaction with jupyter notebooks, so you will also need to install nest_asyncio (1.5.1).

## For development:

1. *Fork* to make your own version.
2. *Git clone* into a local location on your own system.
3. to install to your environment, use *python setup.py develop*
4. Create a branch relevant to what you are working on e.g. *git branch fmriprep*
5. *git checkout fmriprep*
6. Make and commit new changes
7. Test things in notebooks stored in cnbpy/test/develop, with this at the top:

%load_ext autoreload
%autoreload 2

8. *git push origin fmriprep:fmriprep*
9. Make a pull request.


## How to use

Functionality so far is shown in cnbpy/test/develop

The basic structure so far is modular:

*datalad.py* - Class for interacting with datalad data.
*bids.py* - Class for interacting with BIDS data.
*scratchpad.py* - defunct stuff that I don't yet want to get rid of.

So there will also be:

*fmriprep.py*
*analyze.py*
*visualise.py*



