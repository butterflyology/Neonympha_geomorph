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


    nvidia-docker run --rm -p localhost:8888 -v  butterflyology/fastai_jupyter jupyter notebook
