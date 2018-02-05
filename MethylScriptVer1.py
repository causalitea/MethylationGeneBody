# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 00:31:41 2017

Finds gene body methylation on the genome

@author: Daniel
"""
import os

from datetime import datetime

import csv
print("Modules imported")

mapfile='biomartNoRef.csv'
methylfile='Hyper_bed_CSV.csv'

startTime= datetime.now()
line=0
lines=28156
methylcount=0
notfoundcount=0
didintfindgene=0 #keep track of how many genes dont get found

# makes the list of chromosome possibilitities
listofchr=[]
listofchrnotappended=[]
print("lists made")
# get chr names from locus file
with open(mapfile,'r') as locusstart:
    locusstartRead= csv.reader(locusstart)
    for row3 in locusstartRead: 
        if row3[4] not in listofchr:
            listofchr.append(row3[4])   
print("list one edited")
# get chr names from methyl file
#also count which chr are not in both files
with open(methylfile,'r') as hyperlist:
    hyperlistRead= csv.reader(hyperlist)
    for row4 in hyperlistRead:
        if row4[0] not in listofchr:
            listofchr.append(row4[0])
            listofchrnotappended.append(row4[0])

print("list two edited, starting file creation function")
# takes list and makes a csv for each and fills the csv with its data

for string in listofchr:
    name=string+"file"
    with open(string+"CSV"+".csv",'a+') as name:
        nameWriter= csv.writer(name,lineterminator="\n")
        with open(mapfile,'r') as locusstart:
            locusstart2Read= csv.reader(locusstart)
            for x in locusstart2Read:
                chrm=x[4]
                if chrm == string:
                    nameWriter.writerow(x)
    print("Created",string,"CSV file from methylation file and gene map...")
print("lists created in",datetime.now()-startTime)
print("Calling methyliation file")        

#calls the hypermethylation file        
with open(methylfile,'r') as hyper:
    hyperRead= csv.reader(hyper)
    for row in hyperRead:
        line=line+1
        chromosome=row[0]
        methposition=float(row[1])
        methposition2=float(row[2])
        percent=row[3]
        print("Checking line",line, "in Hypermethylation")
        print("Chromosome",chromosome)
        print("At position",methposition,methposition2)
        filename=chromosome+"CSV.csv"
        print("opening:",filename)
        print("in directory",os.getcwd())
        with open(filename,'r') as name:
            print("opened filename")

            
            # now were going to check our Name CSV map file for a match with Hyper                                
            nameREAD= csv.reader(name)
            for row2 in nameREAD:
                chromosome2=row2[4]
                start=float(row2[0])
                end=float(row2[1])
                strand=row2[2]
                genename=row2[3]
                length=end-start
                notfoundcount=notfoundcount+1
#                print("Matching meth",methposition,"chromosome",chromosome,"mapchr",chromosome2,"start",start,"end",end)
                if methposition >=start and methposition <= end:
                    print("line",line,"Found in chromosome",chromosome)
                    notfoundcount=0 #reset not found count
                    methylcount=methylcount+1
                    with open("CONDENSEDOUTPUT.csv","a+") as condensed:
                        condensedWrite=csv.writer(condensed,lineterminator="\n")
                        condensedWrite.writerow([percent,genename,chromosome,methposition,length,percent])
                        print(" ")
                        print(100*line/lines,"%","line",line,"of",lines)
                        print(" ")
                        print("elapsed",datetime.now()-startTime)
                        break
                            
            else:
                print("Did not find position",chromosome,methposition,methposition2)
                didintfindgene=didintfindgene+1
                notfoundcount=0
                with open("CONDENSED_NOT_FOUND.csv","a+") as notcondensed:
                    NotcondensedWrite=csv.writer(notcondensed,lineterminator="\n")
                    NotcondensedWrite.writerow([percent,chromosome,methposition,methposition2])                        
                    

# We are going to find what percent methylation each gene has
#for q in listofchr:
#    with open(q+"_OUTPUT"+methylfile,r):
#        
        
    
    





print("Time for MISSING READS")


Errorfile="Chr_Not_Found.csv"
 # Lets find all the methylation reads that are not in 
with open(methylfile,'r') as methylfileerror, open(Errorfile, 'a+') as errorfile:
    errorcount=0
    methylerrorRead=csv.reader(methylfileerror)
    methylerrorWriter=csv.writer(errorfile,lineterminator="\n")
    for something in methylerrorRead:
        if something[0] in listofchrnotappended: #if methylation chromosome is in the missing chr list then report this line as missing
            errorcount=errorcount+1
            methylerrorWriter.writerow(something)
            print("Missing read found")
            
            
             

# with open('a', 'w') as a, open('b', 'w') as b:
#    do_something()           
            
endTime=datetime.now()                    
print(".............................")                    
print("........Done")
print(" ")
print("Started",startTime)
print("Ended",endTime)
print("Total")
print(endTime-startTime)
print(" ")
print("These genes could not be found in map file and were not accounted for:")
print(listofchrnotappended)
print("Refer to"+ Errorfile,"for missing reads")
print("Number of missing reads =",errorcount)
print("Percent of methylation file reads not counted",100*errorcount/lines)
print("This attempt found",methylcount,"methylation reads out of",lines)
print("This is",100*(lines-methylcount)/lines,"%")
print("The number of genes that could not be found is", didintfindgene)

