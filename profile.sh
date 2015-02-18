#!/usr/bin/env sh

python -m cProfile -o profile.out ./show_me_pdbs.py demo_data/0001 -qf              
runsnake profile.out 

