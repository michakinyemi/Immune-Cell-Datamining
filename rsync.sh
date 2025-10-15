
rsync -rav --exclude='.git' --exclude='data/' --exclude='tests/' --exclude='results/' --exclude="renv/" \
      --exclude='rsync.sh' --exclude='output/' --exclude='external/'  --exclude='*.RData' \
      --exclude='.Rproj.user' --exclude='*.rds' --exclude='*.csv' --exclude='init.R' \
      --exclude='init.R' \
       ../Immune-Cell-Datamining/ mi542876@stokes.ist.ucf.edu:nguyen_lab/Immune-Cell-Datamining/
