#!/bin/bash

SRC_DIRS="antrag gen falcon ntru"
OUT_DIR=supercop-build

rm -rf supercop-build
mkdir supercop-build
for SRC_DIR in $SRC_DIRS
do
	for SRC in $(ls $SRC_DIR)
	do
		sed -r \
			-e "s/(#include \")([a-z][a-z.]+)/\1${SRC_DIR}_\2/" \
			-e 's/(#include ")..\/([a-z]+)\//\1\2_/' \
			$SRC_DIR/$SRC > $OUT_DIR/${SRC_DIR}_$SRC
	done
done
rm $OUT_DIR/antrag_main.c # not useful; entrypoint is supercop/sign.c
rm $OUT_DIR/antrag_randombytes.c # provided by supercop itself + supercop/sign.c

cp supercop/api.h supercop/sign.c $OUT_DIR

echo "done."
