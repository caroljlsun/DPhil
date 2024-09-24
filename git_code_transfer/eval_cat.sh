#!/usr/bin/bash 

# read number of guide regions to design


ROUND_PATH="/fsx/home_dirs/czjs/doctoral_work/2023/oligo/attempt_left_right"


GUIDE_FILE_PATH="$ROUND_PATH""/leftright_coords_trial.txt"

GUIDES_LEFT=$(cat $GUIDE_FILE_PATH)

echo $GUIDES_LEFT