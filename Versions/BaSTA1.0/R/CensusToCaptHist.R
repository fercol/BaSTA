CensusToCaptHist <-
function(ID,d,dformat="dd/mm/yyyy"){

# Check data
if(class(ID) != "character") ID = as.character(ID)


if(class(d) == "Date") {yr = as.numeric(format(d, format = "%Y")) }



if(class(d) != "Date"){
    # Deduce the d format, if it is not already a d object
    st = strsplit(dformat, split = NULL)
    y  = which(st[[1]]=="y")

    # Extract the year
    yr = as.numeric(substr(d,y[1],y[length(y)]))
    }


# Add a dummy individual seen in all years, 
# to cope with years where there are no recaptures
dyr = min(yr):max(yr)
ID  = c(ID,rep("XdummyX",length(dyr)))
yr  = c(yr,dyr)

mat        = as.matrix(table(ID,yr))
mat[mat>0] = 1

#Remove the dummy row
mat        = mat[-which(rownames(mat)=="XdummyX"),]
return(mat)
}

