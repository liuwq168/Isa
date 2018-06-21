# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 15:16:15 2016

@author: Lwq
"""

import time
import os,sys
from optparse import OptionParser
import pandas as pd
import random
from scipy import stats
import multiprocessing

"""
# =============================================================================

           S T A R T   D E F I N I T I O N S 

# =============================================================================
"""
# =============================================================================
def prepare_optparser():
	usage = "usage: python %prog -f <genome> -o <out_dir> [options]"
	description = "Use varscan to create a pipeline"
	optparser = OptionParser(version="%prog -- 1.0",description=description,\
	usage=usage,add_help_option=False)
	optparser.add_option("-f","--Genome",dest="Reference",type="string",\
	default=os.path.abspath("~/reference/"),
	help="Input reference genome file",metavar="<file>")
	optparser.add_option("-o","--out_dir",dest="out_dir",type="string",
                         default=os.path.abspath("./CNV"),
                         help="Output directory,default=current_directory/bioSV",\
                         metavar="<directory>")
	optparser.add_option("-t","--tum bam",dest="tum_bam",type="string",\
	help="Input tumor bam file",metavar="<Tbam>")
	optparser.add_option("-n","--nor bam",dest="nor_bam",type="string",\
	help="Input normal bam file",metavar="<Nbam>")
	optparser.add_option("-c","--config",dest="config",type="string",\
	help="Input config file",metavar="<config>")
	optparser.add_option("--MC",dest="MC",default=10,type="int",
                         help="Minimum coverage for duplication \
                         calling [10]",metavar="<float>")
	optparser.add_option("--MVF",dest="MVF",default=0.08,type="float",
                         help="Minimum variant frequency for duplication \
                         calling [0.08]",metavar="<float>")
	optparser.add_option("--SV",dest="SV",default=0.05,type="float",
                         help="Somatic value for duplication \
                         calling [0.05]",metavar="<float>")
	optparser.add_option("--samtools",dest="samtools",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/samtools"),
                        help="samtools path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/samtools")),metavar="<file>")
	optparser.add_option("--varscan",dest="varscan",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/VarScan.v2.3.7.jar"),
                        help="varscan path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/VarScan.v2.3.7.jar")),metavar="<file>")
	optparser.add_option("--lumpyexpress",dest="lumpyexpress",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/lumpyexpress"),
                        help="lumpyexpress path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/lumpyexpress")),metavar="<file>")
	optparser.add_option("--extractSplitReads_BwaMem",dest="extractSplitReads_BwaMem",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/lumpy-sv/scripts/extractSplitReads_BwaMem"),
                        help="extractSplitReads_BwaMem path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/lumpy-sv/scripts/extractSplitReads_BwaMem")),metavar="<file>")
	optparser.add_option("--pairend_distro",dest="pairend_distro",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/lumpy-sv/scripts/pairend_distro.py"),
                        help="pairend_distro.py path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/lumpy-sv/scripts/pairend_distro.py")),metavar="<file>")
	optparser.add_option("--lumpy",dest="lumpy",type="string",
                         default=os.path.abspath(sys.path[0]+"/tools/lumpy-sv/bin/lumpy"),
                        help="lumpy path [%s]"%(os.path.abspath(sys.path[0]\
                        +"/tools/lumpy-sv/bin/lumpy")),metavar="<file>")
	optparser.add_option("-h","--help",action="help",\
    	help="Show this help message and exit")
    	return optparser
# =============================================================================
def pyshell(cmd):
        e=os.system(cmd)
        if e!=0:
                raise RunError("Wrong![CMD]\n%s\nRun failed"%cmd)
# =============================================================================
def BICseq2(Tbam,Nbam,Outfinal):
	os.mkdir(Outfinal+'/Tseq')
	os.mkdir(Outfinal+'/Nseq')
	#Generate_seg(Tbam,Nbam):
	Name=['chr'+str(i) for i in (list(range(1,23))+['X','Y'])]
	#BICseq2(Tbam,Nbam):
	pyshell('./samtools view -U BWA,%s/Tseq/,N,N %s'%(Outfinal,Tbam))
	pyshell('./samtools view -U BWA,%s/Nseq/,N,N %s'%(Outfinal,Nbam))
	config1=Outfinal+'/config.tu.norm'
	config2=Outfinal+'/config.no.norm'
	config3=Outfinal+'/config.seg'
	SeqOut1=Outfinal+'/Tseq'
	SeqOut2=Outfinal+'/Nseq'
	typ1='tumo'
	typ2='norm'
	typ3='tumo'
	typ4='norm'
	#Generate_config_norm(config,SeqOut,typ):
	with open(config1,'w') as f:
		f.write('chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n')
		for i in range(len(Name)):
			f.write(Name[i]+'\tRef/'+Name[i]+'.fa\thg19.CRG.50bp/hg19.50mer.CRC.'+Name[i]+'.txt\t'+SeqOut1+'/'+Name[i]+'.seq\t'+Outfinal+'/'+Name[i]+'.'+typ1+'.bin\n')
	with open(config2,'w') as f:
		f.write('chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n')
		for i in range(len(Name)):
			f.write(Name[i]+'\tRef/'+Name[i]+'.fa\thg19.CRG.50bp/hg19.50mer.CRC.'+Name[i]+'.txt\t'+SeqOut2+'/'+Name[i]+'.seq\t'+Outfinal+'/'+Name[i]+'.'+typ2+'.bin\n')
	#Generate_config_seg(config,typ1,typ2):
	pyshell('NBICseq-norm_v0.2.4/NBICseq-norm.pl %s/config.tu.norm %s/OutputTu'%(Outfinal,Outfinal))
	pyshell('NBICseq-norm_v0.2.4/NBICseq-norm.pl %s/config.no.norm %s/OutputNo'%(Outfinal,Outfinal))
	with open(config3,'w') as f:
		f.write('chromName\tbinFileNorm.Case\tbinFileNorm.Control\n')
		for i in range(len(Name)):
			f.write(Name[i]+'\t'+Outfinal+'/'+Name[i]+'.'+typ3+'.bin\t'+Outfinal+'/'+Name[i]+'.'+typ4+'.bin\n')
	pyshell('NBICseq-seg_v0.7.2/NBICseq-seg.pl %s/config.seg %s/OutputSeg'%(Outfinal,Outfinal))
# =============================================================================
def ControlFreec(Tbam,Nbam,ConfigFILE,Outfinal):
	#Generate_config(Tbam,Nbam,ConfigFILE):
	os.mkdir(Outfinal+'/ControlFreec')
	f=open(ConfigFILE,'w')
	f.write('[general]\n')
	f.write('breakPointThreshold = 0.5\n\n')
	f.write('chrFiles = /media/_EXTend2016/LiuWT/CNV/reference/UCSC_Chromosomes\n')
	f.write('chrLenFile = /media/_EXTend2016/LiuWT/CNV/reference/hg19.len\n\n')
	f.write('coefficientOfVariation = 0.04\n\n')
	f.write('degree = 3\n\n')
	f.write('forceGCcontentNormalization = 0\n\n')
	f.write('minCNAlength = 1\n\n')
	f.write('minMappabilityPerWindow = 0.20\n\n')
	f.write('minimalSubclonePresence = 20\n\n')
	f.write('readCountThreshold = 8\n\n')
	f.write('maxThreads = 5\n\n')
	f.write('outputDir = %s/ControlFreec\n\n'%(Outfinal))
	f.write('ploidy = 2\n\n')
	f.write('printNA = FALSE\n\n')
	f.write('readCountThreshold = 8\n\n')
	f.write('step = 100\n\n')
	f.write('window = 100\n\n')
	f.write('[sample]\n\n')
	f.write('mateFile = '+os.path.abspath(Tbam)+'\n')
	f.write('inputFormat = BAM\n')
	f.write('mateOrientation = 0\n\n')
	f.write('[control]\n\n')
	f.write('mateFile = '+os.path.abspath(Nbam)+'\n')
	f.write('inputFormat = BAM\n')
	f.write('mateOrientation = 0\n\n')
	f.close()
	#Calculate_freec(ConfigFILE):
	pyshell("freec -conf %s"%(ConfigFILE))
# =============================================================================
def cnMOPS(Tbam,Nbam):
	pyshell("Rscript cn.mops.modify.1230.R %s %s"%(Tbam,Nbam))
# =============================================================================
def saasCNV(Tbam,Nbam,Outfinal):
	Outmpileup=Outfinal+'/'+Tbam.replace("bam","")+"normal.mpileup"
    	Outvcf=Outmpileup.replace("mpileup","vcf")
    	OutFormatvcf=Outvcf.replace("vcf","format.vcf")
	#Calculate_saasCNV(Tbam,Nbam,Outmpileup,Outvcf,OutFormatvcf):
	pyshell("samtools mpileup -Bf ../reference/ucsc.hg19.fa %s %s > %s"%(Tbam,Nbam,Outmpileup))
    	pyshell("java -jar VarScan.v2.3.7.jar mpileup2snp %s --min-var-freq 0 --p-value 0.99 --min-avg-qual 1 --min-coverage 10 --output-vcf 1 > %s"%(Outmpileup,Outvcf))
    	df=pd.read_csv(Outvcf,skiprows=23,sep='\t')
    	df1=df[df['Sample1'].apply(lambda x:'./.:.:' not in x)]
    	df2=df1[df1['Sample2'].apply(lambda x:'./.:.:' not in x)]
    	df2.rename(columns={'#CHROM':'CHROM'},inplace=True)
    	df2['MQ']=df['QUAL'].apply(lambda x:x)
    	df2['Tumor.GT']=df2['Sample1'].apply(lambda x:x.split(':')[0])
    	df2['Tumor.REF.DP']=df2['Sample1'].apply(lambda x:x.split(':')[4])
    	df2['Tumor.ALT.DP']=df2['Sample1'].apply(lambda x:x.split(':')[5])
    	df2['Normal.GT']=df2['Sample2'].apply(lambda x:x.split(':')[0])
    	df2['Normal.REF.DP']=df2['Sample2'].apply(lambda x:x.split(':')[4])
    	df2['Normal.ALT.DP']=df2['Sample2'].apply(lambda x:x.split(':')[5])
    	df3=df2[list(df2.columns[:6])+list(df2.columns[-7:])]
    	df3.to_csv(OutFormatvcf,header='infer',index=None,sep='\t')
    	pyshell("Rscript saasCNV.new.1213.R %s"%(OutFormatvcf))
# =============================================================================
def Lumpy(Tbam,Nbam,Outfinal):	
	Tdubam=Outfinal+'/'+out_dir+"/"+Tbam.replace("bam","discordants.unsorted.bam")
	Tsubam=Outfinal+'/'+out_dir+"/"+Tbam.replace("bam","splitters.unsorted.bam")
	Ndubam=Outfinal+'/'+out_dir+"/"+Nbam.replace("bam","discordants.unsorted.bam")
	Nsubam=Outfinal+'/'+out_dir+"/"+Nbam.replace("bam","splitters.unsorted.bam")
	Td=Outfinal+'/'+out_dir+"/"+Tbam.replace("bam","discordants")
	Ts=Outfinal+'/'+out_dir+"/"+Tbam.replace("bam","splitters")
	Nd=Outfinal+'/'+out_dir+"/"+Nbam.replace("bam","discordants")
	Ns=Outfinal+'/'+out_dir+"/"+Nbam.replace("bam","splitters")
	Tb=Tbam.split('.')
	Nb=Nbam.split('.')
	Outvcf=Outfinal+'/'+out_dir+'/'+Tb[0]+'_'+Tb[1]+'_'+Nb[1]+'.vcf'
	#runPreLumpy(Tbam,Nbam,Tdubam,Tsubam,Ndubam,Nsubam,Td,Ts,Nd,Ns,Outvcf):
	pyshell("samtools view -b -F 1294 -@ 16 %s > %s"%(Tbam,Tdubam))
	pyshell("samtools view -h -@ 16 %s | "%(Tbam)+extractSplitReads_BwaMem+" -i stdin | samtools view -Sb -> %s"%(Tsubam))
	pyshell("samtools sort -m 5G -T PREFIX -@ 10 %s -o %s"%(Tdubam,Td))
	pyshell("samtools sort -m 5G -T PREFIX -@ 10 %s -o %s"%(Tsubam,Ts))
	pyshell("samtools view -b -F 1294 -@ 16 %s > %s"%(Nbam,Ndubam))
	pyshell("samtools view -h -@ 16 %s | "%(Nbam)+extractSplitReads_BwaMem+" -i stdin | samtools view -Sb -> %s"%(Nsubam))
	pyshell("samtools sort -m 5G -T PREFIX -@ 10 %s -o %s"%(Ndubam,Nd))
	pyshell("samtools sort -m 5G -T PREFIX -@ 10 %s -o %s"%(Nsubam,Ns))
	pyshell("lumpyexpress -B %s,%s -S %s,%s -D %s,%s -o %s"%(Tbam,Nbam,Ts,Ns,Td,Nd,Outvcf))
# =============================================================================
"""
# =============================================================================

           E N D    O F    D E F I N I T I O N S 

# =============================================================================
"""
if __name__ == "__main__":
	sys.stdout.write(time.strftime('%Y-%m-%d-%H:%M',time.localtime())+'\n')
    	optparser=prepare_optparser()
    	(options,args) = optparser.parse_args()	
	if not (options.out_dir):
        	optparser.print_help()
        	sys.exit("-o must be assigned")
    	if not os.path.exists(options.out_dir):
        	os.mkdir(options.out_dir)
	if not (options.tum_bam):
        	optparser.print_help()
        	sys.exit("-t must be assigned")	
	out_dir=options.out_dir
    	#samtools=options.samtools
	extractSplitReads_BwaMem=options.extractSplitReads_BwaMem
	#lumpyexpress=options.lumpyexpress
	Tbam=options.tum_bam
	Nbam=options.nor_bam
	Outfinal=Tbam.split('.')[0].split('T')[0]+'_FinalOut'
	os.mkdir(Outfinal)
	os.mkdir(Outfinal+'/'+out_dir)
	configFILE=options.config
	ConfigFILE=Outfinal+'/'+configFILE
	#=============================================================================
	p1 = multiprocessing.Process(target = BICseq2, args = (Tbam,Nbam,Outfinal,))
	p2 = multiprocessing.Process(target = ControlFreec, args = (Tbam,Nbam,ConfigFILE,Outfinal,))
	p3 = multiprocessing.Process(target = cnMOPS, args = (Tbam,Nbam,))
	p4 = multiprocessing.Process(target = saasCNV, args = (Tbam,Nbam,Outfinal,))
	p5 = multiprocessing.Process(target = Lumpy, args = (Tbam,Nbam,Outfinal,))
	p1.start()
	p2.start()
	p3.start()
	p4.start()
	p5.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()
	p5.join()
	#=============================================================================
	sys.stdout.write(time.strftime('%Y-%m-%d-%H:%M',time.localtime())+'\n')
    	sys.stdout.write("Run sucessfully\n")

#python2.7 Integrate.CNV.1213.py -o lumpy_out -t 101T.add.sort.bam -n 101N.add.sort.bam -c config.txt
