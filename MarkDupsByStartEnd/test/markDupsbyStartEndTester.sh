#!/bin/bash

shopt -s expand_aliases

alias markdups='java -jar ~/Tritume/MarkDupsByStartEnd.jar'

echo "TEST: CAN SHOW HELP"
markdups -h 

echo "TEST: DOES NOTHING MARKED OR CHANGED"
markdups -i test_data/testfile-1.sam > out.sam
diff test_data/testfile-1.sam out.sam # diff should return nothing

echo "TEST: IGNORE READ GROUP"
markdups -i test_data/testfile-1.sam -rg > out.sam
# ... TODO


echo "TEST: SAM OUTPUT PASSES Picard Validation"
# ...

echo "TEST: TAG FIELDS ARE ALL RETURNED"
# ...

echo "TEST: CAN READ FROM STDIN"
# ...