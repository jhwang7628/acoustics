#!/bin/bash

# add all cpp, h, hpp extensions and CMakeLists.txt in all directories
git add *.cpp *.h *.hpp *.sh '*CMakeLists.txt'
git reset src/modal_model/ModalAnalysis.*
