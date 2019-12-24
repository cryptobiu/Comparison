#! /bin/bash
for i in `seq $1 1 $2`;
do	
	./Comparison  -partyID $i -partiesNumber $3    -fieldType $4 -partiesFile Parties.txt  -internalIterationsNumber 1 numThreads 1 -lambda $5 -delta $6 &
	echo "Running $i..."
done
