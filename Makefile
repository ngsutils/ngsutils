init:
	./init.sh

update:
	git pull; ./init.sh -q
