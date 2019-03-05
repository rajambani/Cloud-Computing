#!/bin/bash

name=$(echo $1 | cut -d '-' -f1 -)
node1=$(echo $1 | cut -d '-' -f2 - | tr -d '[')
node2=$(echo $1 | cut -d '-' -f3 - | tr -d ']')


if [ "$(hostname)" == "$name-$node1" ]
then
    echo " $(hostname) " > server.txt
else
    echo "I am the client, and the server is $name-$node1." > $(hostname).txt
fi

