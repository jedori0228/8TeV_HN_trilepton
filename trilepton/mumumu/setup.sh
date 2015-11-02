#!/bin/bash
perl -pi -e 's/Users/home/g' *.C
perl -pi -e 's/new_sample\///g' *.C
perl -pi -e 's/plots_withcut/plots/g' *.C
