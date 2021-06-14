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


#for (f in 1:length(files.v)){
for (f in 128:130){
  drop_upload(file=paste(deliv.dir, files.v[f], sep="/"),
              #path="request/KjNv7oMaN9VqdyObg1Vr",
              path="VIMC/Upload_2019",
              dtoken=db.token)
}

# Uploaded first file as a test run on 2020.01.16
# Waiting for confirmation from Dinithi

# Uploading directly to VIMC dropbox link doesn't work.
# Upload to personal account, then from outside firewall tranfer from within dropbox
# As of 2020/02/01, uploaded files   1 -   5 in the files.v list
# As of 2020/02/14, uploaded files   6 -  99 in the files.v list
# As of 2020/02/17, uploaded files 100 - 130 in the files.v list