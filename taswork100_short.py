#This file is a translation from perl to python of taswork.pl
#Translator: Ngoc Duong
#Date: 01/16/2019

import sys
import os
from subprocess import PIPE, Popen, check_call, STDOUT
import operator
import random #to test, should del
#argv[1]=xfile argv[2]=keyfile argv[3]=workingdirectory argv[4]=similarity minimum

#global variables
BASES = ['A','G','C','T']


#Help function to find Levenshtein distaince
def levenshtein(s1, s2):
	if len(s1) < len(s2):
		return levenshtein(s2, s1)
	if len(s2) == 0:
		return len(s1)
	previous_row = range(len(s2) + 1)
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			insertions = previous_row[j + 1] + 1
			deletions = current_row[j] + 1
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]

#for each tag
def getconsensus(testfile):#remove test file
        file = open(testfile, 'a+')
	
	in5 = open(wdir + "/" + rname + ".d/" + holdtag + ".oneline", 'r')
        auth = in5.readline().rstrip()
        authprts = auth.split('\t')
        authname = authprts[0].split("_")
        consens2 = [authname[0] + "_" + authname[2] + "\n" + authprts[1] + "\n"]
        authseq = authprts[1]
        
	#find positions of am codes in auth seq
	amcode_positions = {}
	for s in range(len(authseq)):
		if(authseq[s] != "-" and not authseq[s].upper() in BASES and authseq[s].upper() != "Z"):
			amcode_positions[s] = [0,0,0,0,0,0]
	cx = in5.readline().rstrip()
        while(cx):
        	subseq = cx.split("\t")
        	inseq = subseq[1]
        	consens2.append(subseq[0] + "\n" + subseq[1] + "\n")

        	for s,list in amcode_positions.items():
			if(inseq[s].upper() == "A"):
				amcode_positions[s][1] += 1
               		if(inseq[s].upper() == "C"):
				amcode_positions[s][2] += 1
               		if(inseq[s].upper() == "G"):
				amcode_positions[s][3] += 1
               		if(inseq[s].upper() == "T"):
				amcode_positions[s][4] += 1
               		if(inseq[s].upper() == "-"):
				amcode_positions[s][5] += 1
		cx = in5.readline().rstrip()
	in5.close()

	consensline = ""
        for pl in range(len(authseq)):
        	if(not pl in amcode_positions.keys()):
			consensline = consensline + "."
		else:
			het_dict = {"A":amcode_positions[pl][1], "C":amcode_positions[pl][2], "G":amcode_positions[pl][3], "T":amcode_positions[pl][4], "d":amcode_positions[pl][5]}
        		het = max(het_dict.iteritems(), key=operator.itemgetter(1))[0]
 
			allelecnt = sum([1 for key,value in het_dict.items() if value > 0])
        			
       		 	if(allelecnt > 2):
				het = "N"
        		elif(allelecnt >1):
        			if(het_dict['A'] and het_dict['C']):
					het = "M"
        			if(het_dict['A'] and het_dict['G']):
					het = "R"
        			if(het_dict['A'] and het_dict['T']):
					het = "W"
        			if(het_dict['d'] and (het_dict['A'] or het_dict['G'] or het_dict['C'] or het_dict['T'])):
					het = "X"
        			if(het_dict['C'] and het_dict['G']):
					het = "S"
   		     		if(het_dict['C'] and het_dict['T']):
					het = "Y"
        			if(het_dict['G'] and het_dict['T']):
					het = "K"
        
        		consensline = consensline + het
        consensline = consensline + "\n"
        consstring = authname[0] + "_" + rname + "\n" + consensline
        consens2 = [consstring] + consens2
        #bubble sort consens2 on similarity to key file and get hit count
	
	totalreadcnt = 0
        len_ = len(consens2)
	align = holdtag + "_" + rname + ".fas"
        if(len_ >= 3):
	        out2 = open(wdir + "/" + rname + ".d/" + align, 'w')
	        if(len_ == 3):
			prt1x = consens2[2].split("_")
			totalreadcnt = prt1x[1]
		else:
        		for c3x in range(2, len_):
				prt1x = consens2[c3x].split("_")
				totalreadcnt += int(prt1x[1])
						
		for l in consens2:
			out2.write(l)
        	if(rname != "" and holdtag != "" and totalreadcnt > 0):
			hits.write(rname + "\t" + holdtag + "\t" + str(totalreadcnt) + "\n")
        	out2.close()
        	
		consensline = consensline.rstrip()
        	contrk = 0
		postrk = -1
        	for conprt in consensline:
        		postrk += 1
        		if(authseq[postrk] != "-"):
				contrk += 1
        		if(conprt != "." and postrk < len(consensline)-1):
        			tagpos = holdtag + "_" + str(contrk)
				rept.write(rname + "\t" + tagpos + "\t" + holdtag + "\t" + str(contrk) + "\t" + conprt + "\n")
	os.system("cp " + wdir + "/" + rname + ".d/" + align + " " + wdir + "/Analysis/Markers/" + holdtag + ".d/" + rname + "_" + holdtag + ".fas")
	file.close()
     
#This function verifies that a given file is a fasta or fastq
#exit the program if the file is not fasta or fastq
def verify_file_type(file_path):
	infile = open(file_path, 'r')
	firstline = infile.readline()
	typchk = firstline[0]
	if (typchk == '>'):
		return 'a'
	elif (typchk == '@'):
		return 'q'
	else:
		print ("\n\n FILE MUST BE FASTA OR FASTQ, STARTING WITH > OR @. FIRST CHARACTER OF " + rname + " is " + str(typchk) + "\n\n")
		sys.exit()



#This function makes a directory if it does not exits yet
def make_dir(path):
	if not os.path.exists(path):
			os.makedirs(path)


#This function prints out a help message then exits the program
def print_help():
	if (len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help"):
		print "How to run: python taswork.py [xfile] [cleankey file] [directory of raw data file] [mininal similarity]"
		print "Example   : python taswork.py   x1         cleankey             PRO366                   1.0" 
		sys.exit()


#This function filters out seqs whose length < threshold bps and makes a new file to store the resulting seqs
def filter_out_short_seqs(infile, outfile, threshold):
	ftyp = verify_file_type(infile)
	inraw = open(infile, 'r')
	outraw = open(outfile, 'w')
	rwcnt = 0
	intag = inraw.readline().rstrip()
	while(intag):
		rwcnt += 1
		inrd = inraw.readline().rstrip()
		if (len(inrd) >=threshold):
			outraw.write(">" + str(rwcnt) + "_" + rname + "\n" + inrd + "\n")
		if (ftyp == 'q'):
			inrd3 = inraw.readline()
			inrd4 = inraw.readline()
		intag = inraw.readline()
	inraw.close()
	outraw.close()


#This function returns a tag5p and a tag3p for a given seq
def get_tag5p_tag3p(seq):
	cnt5p = 0	#number of bases before the ambiguity code
	end5p = 0
	end3p = 0
	for in5 in seq:
		if not in5 in BASES:
			end5p = cnt5p	#end5p = index of the ambiguity code
			break
		else:
			cnt5p += 1
	lastbase = len(seq) - 1
	countdown = lastbase #countdown starts at the last index of the sequence
			
	if (end5p == countdown):
		with open(wdir + "/Analysis/Error.log", 'a') as f:
			f.write(inkeyid + " lacks ambiguity codes\n")
		return None, None
	
	while True:
		if not seq[countdown] in BASES:
			end3p = countdown+1
			countdown = 0
		else:
			countdown -= 1
		if (countdown == 0):
			break
			
	return seq[0:end5p], seq[end3p:]



if __name__ == "__main__":
	#print out help message
	print_help()
		
	#get some command arguments
	minsim = float(sys.argv[4])
	xfile = open(sys.argv[1], 'r')
	wdir = sys.argv[3]
	
	#for each read file
	fid = xfile.readline()
	while(fid):
		fid = fid.rstrip()
		rname = os.path.basename(fid)	
		
		make_dir(wdir + "/" + rname + ".d")
		
		#verify file type
		ftyp = verify_file_type(fid)
		
		#filter out seqs whose length < 37 bps
		filter_out_short_seqs(wdir + "/" + rname, wdir + "/" + rname + ".d/" + rname, 37)
		
		#clusterize seqs
		with open(os.devnull, 'wb') as devnull:
			check_call(['cd-hit-est','-M','4096', '-c','1','-r','0','-d','0','-i',wdir + '/' + rname + '.d/' + rname,'-o', wdir + '/' + rname + '.d/' + rname + '.out'], stdout=devnull, stderr=STDOUT)
		os.system("sed -e 's/" + rname + "/0/g' " + wdir + "/" + rname + ".d/\\" + rname + ".out>" + wdir + "/" + rname + ".d/" + rname + ".fou")
		
		
		#####grep out valid seqs######
		inkeys = open(wdir + "/hpcjcl/cleankey", 'r')
		
		rept = open(wdir + "/" + rname + ".d/" + rname + ".snps.txt", 'w')
		hits = open(wdir + "/" + rname + ".d/" + rname + ".hits.txt", 'w')
		
		inkeyid = inkeys.readline()
		while(inkeyid):
			inkeyid = inkeyid.rstrip()
			inkeyseq = inkeys.readline().rstrip()
			#get tag5p, tag3p
			tag5p, tag3p = get_tag5p_tag3p(inkeyseq)		
			if (tag5p == None):
				inkeyid = inkeys.readline()
				continue
			
			holdtag = inkeyid[1:]
			make_dir(wdir + "/Analysis/Markers/" + holdtag + ".d")
			
			with open(wdir + "/" + rname + ".d/" + holdtag + ".hits.tmp", 'w') as f:
				f.write(inkeyid + "\n" + inkeyseq + "\n")

			
			os.system("grep -B1 " + tag5p + " " + wdir + "/" + rname + ".d/" + rname + ".fou| grep -B1 " + tag3p + " |grep -v ^--$ | sed 's#^.*" + tag5p + "#" + tag5p + "#' |sed 's#" + tag3p + ".*$#" + tag3p + "#'>>" + wdir + "/" + rname + ".d/" + holdtag + ".hits.tmp")
			
			#############################
			#collapse the file to unique reads with new counts	
			holdids = []
			holdseqs = []
			intmp = open(wdir + "/" + rname + ".d/" + holdtag + ".hits.tmp", 'r')
			tagid = intmp.readline().rstrip()
			tagrd = intmp.readline().rstrip()
			#count the number of am. codes in this key seq
			ambcnt = 0	
			for t in tagrd:
				if not t in BASES:
					ambcnt += 1
			#move to next key seq if no am. code
			if (ambcnt == 0):
				inkeyid = inkeys.readline()
				continue

			#compare read seqs against the key seq
			inid = intmp.readline()
			while (inid):
				inid = inid.rstrip()
				inseq = intmp.readline().rstrip()
				comdis = levenshtein(tagrd,inseq)-ambcnt
				inid2 = inid + "_" + str(comdis)
				if (comdis == 0 and not inseq in holdseqs):
					holdids.append(inid2)
					holdseqs.append(inseq)
				inid = intmp.readline()
			
			seqcnt = len(holdseqs)
			totalcount = 0
			
			#get total read count
			trkseqs = 0
			fnlseqscnt = len(holdseqs)
			
			while (trkseqs < fnlseqscnt):
				qseq = holdseqs[trkseqs].rstrip()
				gpcnt = Popen(["grep", "-c", qseq, wdir + "/" + rname + ".d/" + rname], stdout = PIPE)
				gpcnt = int(gpcnt.communicate()[0].split()[0])
				totalcount += gpcnt
				idbits = holdids[trkseqs].split('_')
				holdids[trkseqs] = idbits[0] + "_" + str(gpcnt) + "_" + idbits[2]
				trkseqs += 1
					
			trkids = -1
			for r1 in holdids:
				trkids += 1
				idbits = r1.split("_")
				if (len(idbits) > 1 and float(idbits[1]) <= totalcount*0.02 ):
					holdids[trkids] = ""
					holdseqs[trkids] = ""

			holdids.insert(0, tagid)
			holdseqs.insert(0, tagrd)


			howmanyaligns = len(holdids)
			if(howmanyaligns == 1):
				inkeyid = inkeys.readline()
				continue
			
			outseq = open(wdir + "/" + rname + ".d/" + holdtag + ".hits2", 'w')
			for prntseq in range(seqcnt+2):
				if(prntseq < len(holdseqs) and holdseqs[prntseq] != ""):
					outseq.write(holdids[prntseq] + "\n" + holdseqs[prntseq] + "\n")
			intmp.close()
			outseq.close()
			
			#process one tag
			
			rawalignname = rname + "_" + holdtag + ".RAW" + ".align" + ".fasta"
			with open(os.devnull, 'wb') as devnull:
				check_call(['muscle', '-in', wdir + '/' + rname + '.d/' + holdtag + '.hits2', '-out', wdir + '/' + rname + '.d/' + rawalignname], stdout=devnull, stderr=STDOUT)
					
			in4 = open(wdir + "/" + rname + ".d/" + rawalignname, 'r')
			out1 = open(wdir + "/" + rname + ".d/" + holdtag + ".onelineraw", 'w')
			
			line1 = in4.readline().rstrip()
			if(line1[1:] == holdtag):
				line1 = line1 + "_1000000_key"
				out1.write(line1 + "\t")
			else:
				out1.write(line1 + "\t")
			
			lin4 = in4.readline()
			while(lin4):
				lin4 = lin4.rstrip()
				chkstrt = lin4[0]
				if(chkstrt == ">"):
					if(lin4[1:] == holdtag):
						lin4 = lin4 + "_1000000_key"
						out1.write("\n" + lin4 + "\t")
					else:
						out1.write("\n" + lin4 + "\t")
				else:
					out1.write(lin4)
				lin4 = in4.readline()
			
			out1.close()
			in4.close()
			
			os.system("sort -t_ -k2,2rn " + wdir + "/" + rname + ".d/" + holdtag + ".onelineraw>" + wdir + "/" + rname + ".d/" + holdtag + ".oneline")
			getconsensus(wdir + "/" + rname + ".d/test")
			
			inkeyid = inkeys.readline() #end of while(inkeyid)
		
		hits.close()
		rept.close()
			
		os.system("rm " + wdir + "/" + rname + ".d/*.tmp " + wdir + "/" + rname + ".d/*.onelin* " + wdir + "/" + rname + ".d/*.hits2 " + wdir + "/" + rname + ".d/*.fasta " + wdir + "/" + rname + ".d/" + rname + ".fou " + wdir + "/" + rname + ".d/" + rname + ".out " + wdir + "/" + rname + ".d/" + rname + ".out.clstr ")
		os.system("mv " + wdir + "/" + rname + ".d " + wdir + "/Analysis/Genotypes/")
		
		fid = xfile.readline()
#####getconsensus function/definition is moved to top because of python order of execution
