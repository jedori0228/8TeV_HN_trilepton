#!/bin/bash
## (sample, PU, cut, muonidsel) ##

for sample in 40 50 60
  do
  for cut in 22
  do
    # loose
    root -l -q -b "chi2_lambda_stack_mumumu.C+($sample, 1, $cut, 0)"
  done
done
