#!/bin/sh
# convert global dam locations to regional map
# IX = IXX - offsetX
# IY = IYY - offsetY

# global map name
GLB="glb_06min"

# regional map name
REG="conus_06min"

# input dam location file
INDAM="../dat/damloc_${GLB}.txt"

# output dam location file
OUTDAM="../dat/damloc_${REG}.txt"

CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

#
python src/regionalize_damloc.py $GLB $REG $INDAM $OUTDAM $CaMa_dir

wait