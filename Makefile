all: init
init:
	./init.sh
test:
	./test.sh
clean:
	find . -name '*.pyc' -exec rm \{\} \;
	rm -rf ./env

