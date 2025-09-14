
rsync -rav --exclude='.git' --exclude='data/' --exclude='tests/' --exclude='results/' \
      --exclude='rsync.sh' --exclude='output/' --exclude='external/'  --exclude='*.RData' \
      --exclude='.Rproj.user' --exclude='*.rds' --exclude='*.csv' \
       ../Immune-Cell-Datamining/ mi542876@stokes.ist.ucf.edu:nguyen_lab/Immune-Cell-Datamining/
