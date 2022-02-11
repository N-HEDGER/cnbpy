#!/bin/bash

export FS_LICENSE=---fsli---
export SUBJECTS_DIR=---sdir---

mri_surf2surf --srcsubject ---sub--- --srcsurfval ---src--- --trgsubject ---trg--- --trgsurfval ---trgsurf--- --hemi ---hemi---

