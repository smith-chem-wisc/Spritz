#!/bin/bash

function validate(){
if [[ `wget --spider $1 2>&1 | grep 'exists'` ]]; then return 0; else return 1; fi
}

if validate $1; then echo "true"; else echo "false"; fi
