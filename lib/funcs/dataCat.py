#!/usr/bin/python

##########################################################
#
#    Package to catologue file names, graphs, etc.
#    
#   When a new shot is added, include its details here 
#   where necessary (e.g., vpolID for bad pertrubation data
# 
##########################################################


def CatalogueCall(shotid,timeid,runid):
    fileCat={}
    fileCat["GTEDGEsupp"]="Inputs/GTEDGEsupp_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3Consts"]="Inputs/GT3Consts_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3Profs"]="Inputs/GT3Profs_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3NBI"]="Inputs/GT3NBI_"+str(shotid)+"_"+str(timeid)+".csv"

    fileCat["Erspl"]="p"+str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_ersplrhob_sc.txt"
    fileCat["toplasma"]=str(shotid)+"_"+str(timeid)+"_toplasma"
    fileCat["bpol"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_bpol.txt"
    fileCat["btor"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_btor.txt"
    fileCat["R"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_Z.txt"
    fileCat["BDRY"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_BDRY.txt"
    fileCat["EPOTEN"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_EPOTEN.txt"
    fileCat["lim"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_lim.txt"
    fileCat["psirz"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_psirz.txt"
    fileCat["Z"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_Z.txt"
    
    vTorPoints={(166606,1950,"j8099"):[198],
            (144977,3000,"j3000"):[193],
            (118890,1515,"r90"):[181],
            (118890,1560,"r90"):[191]}
    
    try:
        fileCat["vtorID"]=vTorPoints[(shotid,timeid,runid)]
    except:
        pass
    return fileCat
    
    

        
        