#!/bin/bash
./src/LSSCode/bin/LSS_BSK_calc -input $1 -output $2 -beta $3 -printinfo False -numNNB 300
cp ./xdl_beta_skeleton/*.BSKIndex ./data/
rm -r ./xdl_beta_skeleton