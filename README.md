## MINUANO | Notts JWST Pipeline Config

NIRCAM BASIC REDUCTION

# 0. Configurations
    
## 0.1 CRDS - Reference Files used by JWST pipeline

    export CRDS_PATH=<your-path-to-crds-cache>/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

    export CRDS_CONTEXT='jwst_0942.pmap'

    double check if CRDS_CONTEXT is the latest here: https://jwst-crds.stsci.edu/ 

## 0.2 Conda Environment

        conda create -n <env_name> python
        conda activate <env_name>
        pip install jwst==1.5.2

        Running 1.5.2 (instead of 1.5.3) due to a bug in SkyMatchStep: https://github.com/spacetelescope/jwst/issues/6920

## 0.3 Folder Structure
        /ASNs -> json with the groupings necessary for mosaic production (ceers5_f115w.json, for example)
        /uncal -> raw files
        /calibrated -> stage 2 cal files
        /skymatch -> stage 2 with custom sky subtraction (ASN for single sky subtraction)
        /science -> stage 3 mosaics

# 1. NIRCAM BASIC
 
## 1.1 STAGE1: Detector reduction *uncal -> *rate, *rateint. 

    find uncals/ -name '*uncal.fits' | xargs -n1 -P24 -I {} strun STAGE1_params.asdf {}

    PERFORMANCE for CEERS Pointing 5 (90 uncals):

        real    19m50.942s
        user    297m27.575s
        sys     67m7.145s


## 1.2 STAGE2: WISP Removal

    Run wisp_removal.py in wisp_templates/ but visually checking the SW bands.

## 1.3 STAGE2: IMAGE Reduction 2

    time find calibrated/ -name '*rate.fits' | xargs -n1 -P24 -I {} strun STAGE2_params.asdf {}

    real    6m45.502s
    user    95m26.210s
    sys     41m47.545s


## 1.3 STAGE2: 1/f Willott Removal

    find calibrated/ -name '*bg.fits' |  xargs -n1 -P24 -I {} python scripts/batch_1overf.py {}


## 1.4 STAGE2: Flat Median Sky Subtraction

    find calibrated/ -name '*cal.fits' | xargs -n1 -P24 -I {} python scripts/batch_bgsub_flat.py {}

### 1.6 STAGE2: 2D BG Removal (agressive)

    find calibrated/ -name '*1overf.fits' | xargs -n1 -P24 -I {} python scripts/remove_2Dbg.py {}

## 1.7 STAGE2: Astrometric Alignment TWEAKWCS TO GAIA.


## 1.4 STAGE3: IMAGE Reduction 3 - Mosaics

    strun STAGE3_params30mas.asdf config_asn.json

    Can also run in parallel for each filter, but beware of resource consumption:
        For ex: cat runSTAGE3 | xargs -n1 -P4 -I {} sh -c {}

    where runSTAGE3 is a list of commands

    if 1.3 is used, .asdf file needs to be updated to skip SkyMatchStep if there is no overlap
    between the cal files. In case of overlap, don't need to skip it but change skymethod to match+global
