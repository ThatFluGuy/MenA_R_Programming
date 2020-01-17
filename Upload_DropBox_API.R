### Try creating program to upload PSA files to dropbox
### via API instead of website

# Preliminary step: log in to dropbox, go to developer mode.
# Created a new app called "MJ_R_app" and put in
# the following URLs:

# From this app, get the key and secret used below to
# create authentication token.

# Added the following packages beyond the usual:
# assertive and all assertive.$$$$$ derivations
# httpuv, promises, later

library(rdrop2)
library(httpuv)

deliv.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Deliverables/Deliverables 2019/"

# First use, need to get token by logging in to website:
db.token <- drop_auth(new_user=TRUE, key="po6we5g0lllrga1", 
                      secret="35am0yvfael0l7q")
saveRDS(db.token, paste(deliv.dir, "db.token.RDS", sep=""))

# Later uses
db.token <- readRDS(paste(deliv.dir, "db.token.RDS", sep=""))

# All files in deliverables directory, need to subset to
# only PSA files
files.v <- list.files(path=deliv.dir, pattern="stochastic_burden_est") 


for (f in 2:length(files.v)){
  drop_upload(file=paste(deliv.dir, files.v[f], sep="/"),
              path="request/KjNv7oMaN9VqdyObg1Vr",
              dtoken=db.token)
}

# Uploaded first file as a test run on 2020.01.16
# Waiting for confirmation from Dinithi