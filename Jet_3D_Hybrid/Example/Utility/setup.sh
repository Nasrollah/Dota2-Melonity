#!/bin/bash

touch loop.log
cat /dev/null > loop.log
chmod 644 loop.log
echo " " > loop.log
ls -al loop.log
max=64
for i in $(seq 0 1 $max)
do
  touch loop.$i.out loop.$i.err
  echo " " > loop.$i.out 
  echo " " > loop.$i.err
  chmod 644 loop.$i.out loop.$i.err
  ls -al loop.$i.out loop.$i.err
done
exit
