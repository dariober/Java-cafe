#!/bin/bash

shopt -s expand_aliases

alias markdups='java -jar ~/Tritume/MarkDupsByStartEnd.jar'

echo "TEST: HELP"
markdups -h 

echo "TEST: NOTHING MARKED OR CHANGED"
markdups -i test_data/testfile-1.sam > out.sam
diff test_data/testfile-1.sam out.sam # diff should return nothing

echo "TEST: "
