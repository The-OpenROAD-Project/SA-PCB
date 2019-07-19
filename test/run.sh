#!/bin/bash

cd ..
make
cd test
cp ../bin/annealer .
./annealer -p ./apte -i 200 -j 10 -t 0.025 -r 0
rm ./annealer
