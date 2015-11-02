#!/bin/bash
perl -pi -e 's/Users/home/g' *.C
perl -pi -e 's/new_sample\///g' *.C
perl -pi -e 's/plots_withcut/plots/g' *.C
mkdir plots
for cut in "before_dR" "dR" "dR_chi2" "dR_chi2_lambda" "dR_W" "control"
do
  mkdir plots/$cut
done

for cut in "before_dR" "dR" "dR_chi2" "dR_chi2_lambda" "dR_W" "dR_nbjets"
do
  mkdir plots/control/$cut
done