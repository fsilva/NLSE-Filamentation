#!/usr/bin/env python
#
# 


#TODO TODO TODO add version check

import sys
import os
import commands
import time


def main():
	start = time.clock()
#Get list of directories with results
	if(len(sys.argv) < 2):
		results = GetResultsList("./SimulacaoNLSEv4/")
	else:
		results = GetResultsList(sys.argv[1])

	if(len(results[0]) == 0):
		print results[1][0]

	commands.getoutput('mkdir Results')

	for i in xrange(len(results[0])):
		try:
			f = open("../Results/"+results[0][i]+".png",'r')
			f.close()
			print 'Directory ' + results[0][i] + ' already processed.'
		except:
			import NLSEv4
			NLSEv4.Process(results[0][i],results[1][i])
		

	print str(time.clock()-start)+" Seconds Elapsed"





def GetResultsList(directory):
	os.chdir(directory)

	ls = commands.getoutput('ls')

	results1 = ls.split('\n')
	results2 = []
	dirList = []

	for res in results1:
		if(res.find('Results_') != -1):
			if(res.find('png') == -1):
				dirList.append(res)
				data = res.split('_')
				results2.append(data[1]+"/"+data[2]+"/"+data[3]+" - "+data[6]+":"+data[7]+":"+data[8])		

	if(len(results2) == 0):
		results2.append('No Results found in target directory')
	return [dirList,results2]




main()
