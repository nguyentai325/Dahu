#!/bin/zsh
if [[ $# -eq 2 ]] { dot -Tpdf $1 -o $2 } elif [[ $# -eq 1 ]] {  dot -Tpdf $1 -o ${1%dot}pdf } else { echo "usage:" $0 "in.dot [out.pdf]" }
