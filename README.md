Code and data repository for geometric morphometrics project on *Neonympha mitchellii* wings.


Major commits:

1. 2017-04-27: Initial commit of `RMarkdown` document rendered with `knitr`.
1. 2018-04-24: Initial commit of `jupyter` notebook and `keras` image classification.
1. 2019-09-19: Using `pytorch` and `fastai` now. Have data sets arranged and categorized. Got the data imported and a `resnet34` model working.
1. 2019-09-20: Fixed the driver issue, purrs like a kitten now. The 4% error model mistakes Nar for Nhe and vice versa. That is awesome.
1. 2019-09-23: Ran the `resnet50` model and got error to 2.5%. Both `34` and `50` only confuse * Neonympha areolata* and *N. helicta*. Started branch `google_images`, added images from net and personal collection
1. 2019-09-24: Further refining the models, tweaking image size and validation set size. Added a bash script to launch the `docker` image.
1. 2019-09-25: Cropped some images to remove borders and ruler edges (which were present at higher frequencies in *N. helicta* images).
1. 2019-10-02: Beginning progressive augmentation ala Jeremy Howard.
1. 2019-10-03: With progression resizing data augmentation get the error down to 9% with the error coming between *N. areolata* and *N. helicta*.
1. 2019-11-03: Branching off `google_images` to make wire frames for the PCA plots. Woops, alrady did it. But I found that the procrustes anova code needs to be updated.  Working on `resnet50` with the resampling paradigm.
1. 2019-11-04: Cleaned up the `R` code to work with the latest `geomorph` and saved the packages and state in `misc/Neonympha_geommorph_packages.txt`.
1. 2019-11-06: Removing the google images doesn't seem to effect the accuracy of the models.
1. 2020-01-07: Split the file structure to allow one directory with the museum only images and another with the supplemented images. Added `museum_images.ipynb` to hold the code from the museum image only analysis.
1. 2020-01-22: Used `Augmentor` package to augment the museum images. 
1. 2020-01-22: Used `Augmentor` package to augment the museum images.
1. 2020-06-09: Updated paths for museum images.
1. 2020-06-16: Having issues with the `ClassConfusion` widget.
1. 2023-06-21: Back at it. Will update the `fastai` code to `v2`.
1. 2023-06-24: Added a `tf` model.