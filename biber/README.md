# Description
This repository contains a docker image that allows running BEB application.

#  BIBER File Directory structure

```
biber
│   README.md
|   Dockerfile
|   requirements
│   biber.py
|   df_helper.py
|   python_encodings.csv
|       
│
└───data
│   │
│   └───2022
│       │   metadata_beb_02082022.csv
|       |   counts_beb_17032022.csv
│       │   counts_beb_17032022.xlsx
│       │   metadata_beb_17032022.xlsx
│   
└───production/
│   counts.csv
│   metadata.csv
│   
│   
└───images
    │   biber.png
    │   ibio.png
```

# Prerequisites
-----

We assume you have installed Docker and it is running.

See the [Docker website](http://www.docker.io/gettingstarted/#h_installation) for installation instructions.
-----
# Build


Steps to build a Docker image:

# Clone this repo
```
git clone https://github.com/bioquimico/biber.git
```

# Go into the directory
 ``` 
 cd biber
 ```
        
# Build docker image
You can build this docker image from a dockerfile using this command
```
docker build -t biberapp:latest .
```

# Run the docker container
Simply enter the following command to run biber application with example data
```
docker run --name biberapp -p 8501:8501 -ti --rm biberapp:latest
```
or in detached mode or in the background

```
docker run --name biberapp -p 8501:8501 -d --rm biberapp:latest
```


**Local development**

 - copy new counts file into a running container
 ```
 docker cp $(pwd)/data/COUNTS_FILE.csv biberapp:/app/data/production/counts.csv 
```
 - copy new metadata file into a running container
 ```
 docker cp $(pwd)/data/METADATA_FILE.csv biberapp:/app/data/production/metadata.csv 
```

 - Mount your working data folder in the container 
  ```
  docker run -p 8501:8501 -ti --rm -v $(pwd)/data/production:/app/data/production biberapp:latest
  ```
to mount ``` counts.csv ``` and ```metadata.csv``` files from the  ```$(pwd)/data/production``` directory.


# Other useful instructions
## Stop container

    docker stop biberapp
    
## Start container

    docker start biberapp
    
## Inspect container

    docker exec -t -i biberapp /bin/bash

# How to extend this image

Here is an example of a custom Dockerfile which replaces the default metadata.csv by a custom metadata:

    FROM biber:latest
    MAINTAINER adm@canessalab.org
    COPY ./metadata.csv /app/data/production/

