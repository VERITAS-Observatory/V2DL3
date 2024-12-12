# V2DL3 Docker file for Eventdisplay

## Using V2DL3 with docker

Docker requires explicit binding of paths into the image. Replace `data_path`
in the following command accordingly.

```bash
docker run --rm -it -v "<data_path>":/data vts-v2dl3 \
           v2dl3-eventdisplay --help
```

Run the image and provide a bash environment:

```bash
docker run --rm -it -v "<data_path>:/data vts-v2dl3 bash
```

## Building the Image

Docker packages are prepared for each release with a Github Action workflow, see the package registry [here](https://github.com/VERITAS-Observatory/V2DL3/pkgs/container/v2dl3).

To build images locally:

```bash
docker build -t vts-v2dl3 .
```
