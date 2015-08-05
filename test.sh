#!/bin/sh -
echo '------------ TEST START ------------' >> binit.log
tail -f binit.log &
trap 'kill $!' EXIT
time python binitGemMap.py --log=DEBUG -m 2 -t 4 -o TEST/output \
	TEST/input/*map.gz
