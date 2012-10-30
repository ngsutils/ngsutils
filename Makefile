all: init
init:
	./init.sh
clean:
	find . -name '*.pyc' --exec rm \{\} \;
	rm -rf ./env

