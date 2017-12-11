#!/bin/bash 
ps aux | grep ${1} | grep -v kill_process | grep -v grep | sh
