# V2DL3 Docker file for Eventdisplay

## Using V2DL3 with the docker image

Docker requires explicit binding of paths into the image. Replace `data_path` in the following command accordingly.

```
docker run --rm -it -v "<data_path>:/data vts-v2dl3 python /V2DL3/pyV2DL3/script/v2dl3_for_Eventdisplay.py --help
```

Run the image and provide a bash environment:
```
docker run --rm -it -v "<data_path>:/data vts-v2dl3 bash
```

## Building the Image

Note that this Dockerfile is prepared for a Github Action workflow. V2DL3 source code is expected to be in the current directory.

```
docker build -t vts-v2dl3 .
```
