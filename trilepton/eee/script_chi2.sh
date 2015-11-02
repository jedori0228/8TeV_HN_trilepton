#!/bin/bash
## (sample, PU, cut, muonidsel) ##

for sample in 40 50 60
  do
  for cut in 0 1 2 3 4
  do
    # loose
    root -l -q -b "chi2_lambda_stack_eee.C+($sample, 1, $cut, 0)"
  done
done
