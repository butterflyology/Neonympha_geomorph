## Notes for project:

Nmi = *Neonympha mitchellii mitchellii*
Nfr = *Neonympha mitchelli francisi*
Nhe = *Neonympha helicta*
Nar - *Neonympha areolata*

Often there will be a `_NJ` signifying state by two-letter code.

For the `docker` image:

    docker pull nvcr.io/nvidia/pytorch:19.09-py3
    nvidia-docker run -it --rm --ipc=host nvcr.io/nvidia/pytorch:19.09-py3
    docker ps
    docker exec -it <container-id> bash
    docker commit -m "fastai_jupyter" <image_id> butterflyology/fastai_jupyter
    docker push butterflyology/fastai_jupyter

This will run the `docker` image and point the volume to the repo:

    nvidia-docker run --rm --ipc=host -p 8888:8888 -v ~/Projects/Neonympha_classification:/workspace butterflyology/fastai_jupyter jupyter notebook

- 2019-09-19: Got the `docker` file working and the data imported and the silly thing crashes with the `resnet34` model because of some driver problem. Adding a flag (--ipc=host) fixed that. Should still look into the drive issue.
- 2019-09-23: Inspired by `pyimagesearch` I wanted to increase the sample of photos in my data set. For the pull from `images.google.com` images my searches were:
  - "Neonympha areolata" -fracisci -mitchellii -helicta -
  Saved the output as `Nar_photos.txt`
  - "Neonympha mitchellii francisci"  -areolata -helicta -
  Saved output as `Nfr_photos.txt`
  - "Neonympha mitchellii mitchellii"  -areolata -helicta -francisci - Saved output as `Nmi_photos.txt`
  - "Neonympha helicta"  -areolata -mitchellii -francisci - Saved output as `Nhe_photos.txt`

- 2019-09-25: Should consider reaching out to people for more images to test. Or, maybe this is better, make an app people can check photos on.
- 2019-09-25: Progressive augmentation beginning with 64 pixel images is done, next 128 and the 256.
- 2019-11-03: Making the `resnet50` model. Makes the same mistakes as `resnet34`.
