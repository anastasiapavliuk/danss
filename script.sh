#!/bin/bash

M=0
root6 'Nevents.cc("99.root")' >> out.txt & echo .q &> /dev/null

I=$(cat out.txt)

let M=$M+$(echo $I | rev | cut -d "." -f1 | rev)

echo $M
