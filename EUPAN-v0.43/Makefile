export SRC_DIR	= src
SUBDIRS = $(SRC_DIR)/bam2cov \
		  $(SRC_DIR)/ccov
all: 
	@for subdir in $(SUBDIRS); \
	do \
	(cd $$subdir && make); \
	done
.PHONY: all

clean:
	@for subdir in $(SUBDIRS); \
	do \
	(cd $$subdir && make clean); \
	done