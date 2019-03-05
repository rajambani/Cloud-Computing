#!/bin/bash

name=$(echo $1 | cut -d '-' -f1 -)
node1=$(echo $1 | cut -d '-' -f2 - | tr -d '[')
node2=$(echo $1 | cut -d '-' -f3 - | tr -d ']')


if [ "$(hostname)" == "$name-$node1" ]
then
    iperf3 -s
else
    iperf3 -c  $(hostname) -P 8 -u
fi

