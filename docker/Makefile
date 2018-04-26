.PHONY: devel build run stop clean

DATA_PATH = /Users/lucas/data/

all: run

devel:
	docker build --compress -t fastore:devel -f Dockerfile.devel ..

build: devel
	docker build --compress -t fastore:latest -f Dockerfile.build ..

run: build
	docker build --compress -t fastore:run -f Dockerfile.run ..
	docker run -i --name=fastore-test -v $(DATA_PATH):/data -t fastore:run

stop:
	docker stop fastore-test

clean: stop
	docker rm fastore-test