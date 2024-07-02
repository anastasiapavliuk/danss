#!/bin/bash

NAME="example"

while [ -n "$1" ]
do
case "$1" in
-n) echo " Name "$2" is used" & NAME=$2 ;;
esac
shift
done

for ((i= 1;i<=20; i++))
do
mkdir dir_$i
echo "This is file inside dir_"$i"...bla-bla-bla" > dir_$i/${NAME}
done

