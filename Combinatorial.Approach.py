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
def pyshell(cmd):
        e=os.system(cmd)
        if e!=0:
                raise RunError("Wrong![CMD]\n%s\nRun failed"%cmd)

def Rename():
	pyshell('cp OutputSeg bicseq2.vcf')
	pyshell('cp CNMOPS.cnvs.vcf cnmops.vcf')
	pyshell('cp tumor_5X.sort.bam_CNVs freec.vcf')
	pyshell('cp tumor_5X_sort_sort.vcf lumpy.vcf')
	pyshell('cp SAASCNV.new.seq.cnv.txt saascnv.vcf')

def Ratio2CN(Ratio):
	return int(round(2**(1+Ratio)))
	
def trans(x):
	if x=='chrX':
		return 23
	elif x=='chrY':
		return 24
	elif x=='chrM':
		return 25
	else:
		return int(x.split('chr')[-1])	
		
def Score3(S1,E1,S2,E2):
    Overlap=abs(max(S1,S2)-min(E1,E2))
    Det=E2-S2
    if Det==0:
         return 100.
    else:
         return 100.*Overlap/Det


def Judge(st):
	if st<2:
		return 'loss'
	else:
		return 'gain'		 
		 
def Deal_info(info):
    s1=info.split(';;')
    ans=[]
    for i in s1:
            if 'CNV' in i:
                    ans.append(i.split('|')[-3].split('=')[-1])
            else:
                    ans.append(str(i.split('|')[-2].split('=')[-1]))
    return ';'.join(ans)

def Set_type(x):
    if ';' not in x:
        return True
    else:
        t=re.findall(r"\d+\.?\d*",x)
        if len(set(t))==1:
            x1=int(t[0])
            if x1>2 and 'gain' in x:
                return True
            elif x1<2 and 'loss' in x:
                return True
            else:
                return False
        else:
            x1=[int(i) for i in t]
            if max(x1)<2 and 'loss' in x:
                return True
            elif min(x1)>2 and 'gain' in x:
                return True
            else:
                return False
				
def CNV(x):
    if ';' not in x:
        if 'gain' in x:
            return 'gain'
        elif 'loss' in x:
            return 'loss'
        else:
            x1=int(x)
            if x1>2:
                return 'gain'
            else:
                return 'loss'
    else:
        if 'gain' in x:
            return 'gain'
        elif 'loss' in x:
            return 'loss'
        else:
            x1=[int(i) for i in x.split(';')]
            if max(x1)<2:
                return 'loss'
            else:
                return 'gain'

def Tools(st):
    st1=st.split(',')
    if (st1.count('freec')>=1 and st1.count('cnmops')>=1) or (st1.count('freec')>=1 and st1.count('saascnv')>=1) or (st1.count('cnmops')>=1 and st1.count('saascnv')>=1) or (st1.count('freec')>=1 and st1.count('cnmops')>=1 and st1.count('saascnv')>=1):
        return True
    else:
        return False

def Score31(S1,E1,S2,E2):
    Overlap=abs(max(S1,S2)-min(E1,E2))
    Det=E2-S2
    if Det==0:
         return 100.
    else:
         return 100.*Overlap/Det		
		 
def bicseq2(bs_vcf,bs_bed):
	df=pd.read_csv(bs_vcf,sep='\t')
	df.head()
	df['CN']=df['log2.copyRatio'].apply(lambda x:Ratio2CN(x))
	df1=df[df['log2.copyRatio'].apply(lambda x:abs(x)>=0.2)]	
	df2=df1[[0,1,2,-1]]
	df3=df2[df2['CN']!=2]
	df3.to_csv(bs_bed,index=None,header=None,sep='\t')

def cnmops(cp_vcf,cp_bed):
	df=pd.read_csv(cp_vcf,sep='\t')
	df.head()
	df['CN']=df['CN'].apply(lambda x:int(x.split('CN')[-1]))
	df1=df[df['CN']!=2]
	df2=df1[[0,1,2,-1]]
	df2.to_csv(cp_bed,index=None,header=None,sep='\t')

def freec(fc_vcf,fc_bed):
	df=pd.read_csv(fc_vcf,header=None,sep='\t')
	df.columns=['chr','start','end','CN','cnv']
	df.head()
	df['chr']=df['chr'].apply(lambda x:"chr"+str(x))
	df1=df[df['CN']!=2]
	df2=df1[[0,1,2,3]]
	df2.to_csv(fc_bed,index=None,header=None,sep='\t')

def lumpy(ly_vcf,ly_bed):
	df=pd.read_csv(ly_vcf,skiprows=32,sep='\t')
	df.head()
	df1=df[df['tumor_10X'].apply(lambda x:'./.:0:0:0' not in x)]
	df2=df1[df1['normal_10X'].apply(lambda x:'./.:0:0:0' in x)]
	df3=df2[df2['#CHROM'].apply(lambda x:len(x)<=5)]
	df4=df3[df3['ALT'].apply(lambda x:len(x)<=5)]
	df4.head()
	df4['END']=df4['INFO'].apply(lambda x:int(x.split(';')[3].split('=')[-1]))
	df4['SU.PE.SR']=df4['tumor_10X'].apply(lambda x:';'.join(x.split(':')[-3:]))
	df4.head()
	df4['ALT']=df4['ALT'].apply(lambda x:x.replace('<INS>','<DUP>'))
	filter=['<DUP>','<DEL>']
	df5=df4[df4['ALT'].apply(lambda x:x in filter)]
	df5['CNV']=df5['ALT'].apply(lambda x:'gain' if x=='<DUP>' else 'loss')
	df5.head()
	df6=df5[[0,1,-3,-1,-2]]
	df6['key']=df6['#CHROM'].apply(lambda x:trans(x))
	df7=df6.sort_values(['key','POS'])
	df8=df7[[0,1,2,3,4]]
	df8.head()
	df8.to_csv(ly_bed,index=None,header=None,sep='\t')

def saasCNV(sc_vcf,sc_bed):
	df=pd.read_csv(sc_vcf,sep='\t')
	df.head()
	filter=['LOH','gain','loss']
	df1=df[df['CNV'].apply(lambda x:x in filter)]
	df1['CNV']=df1['CNV'].apply(lambda x:x.replace('LOH','loss'))
	df1['ratio.mbaf.p.value']=list(map(lambda x:';'.join([str(x[0]),str(x[1]),str(x[2])]),df1[['log2ratio.p.value','log2mBAF.p.value','p.value']].values))
	df1.head()
	df2=df1[[0,1,2,-2,-1]]
	df3=df2[df2['ratio.mbaf.p.value'].apply(lambda x:'0.0;0.0;0.0' in x)]
	df3.to_csv(sc_bed,index=None,header=None,sep='\t')		

def bicseq2_com(bc_vcf,bc_bed):
	df=pd.read_csv(bc_vcf,sep='\t')
	df.head()
	df['CN']=df['log2.copyRatio'].apply(lambda x:Ratio2CN(x))
	df1=df[df['CN']!=2]
	df2=df1[[0,1,2,-1]]
	df2.to_csv(bc_bed,index=None,header=None,sep='\t')

def cnmops_com(cc_vcf,cc_bed):
	df=pd.read_csv(cc_vcf,sep='\t')
	df.head()
	df['CN']=df['CN'].apply(lambda x:int(x.split('CN')[-1]))
	df1=df[df['CN']!=2]
	df2=df1[[0,1,2,-1]]
	df2.to_csv(cc_bed,index=None,header=None,sep='\t')

def freec_com(fc_vcf,fc_bed):
	df=pd.read_csv(fc_vcf,header=None,sep='\t')
	df.columns=['chr','start','end','CN','cnv']
	df.head()
	df['chr']=df['chr'].apply(lambda x:"chr"+str(x))
	df1=df[df['CN']!=2]
	df2=df1[[0,1,2,3]]
	df2.to_csv(fc_bed,index=None,header=None,sep='\t')

def lumpy_com(lc_vcf,lc_bed):
	df=pd.read_csv(lc_vcf,skiprows=32,sep='\t')
	df.head()
	df1=df[df['tumor_10X'].apply(lambda x:'./.:0:0:0' not in x)]
	df2=df1[df1['normal_10X'].apply(lambda x:'./.:0:0:0' in x)]
	df3=df2[df2['#CHROM'].apply(lambda x:len(x)<=5)]
	df4=df3[df3['ALT'].apply(lambda x:len(x)<=5)]
	df4.head()
	df4['END']=df4['INFO'].apply(lambda x:int(x.split(';')[3].split('=')[-1]))
	df4['SU.PE.SR']=df4['tumor_10X'].apply(lambda x:';'.join(x.split(':')[-3:]))
	df4.head()
	df4['ALT']=df4['ALT'].apply(lambda x:x.replace('<INS>','<DUP>'))
	filter=['<DUP>','<DEL>']
	df5=df4[df4['ALT'].apply(lambda x:x in filter)]
	df5['CNV']=df5['ALT'].apply(lambda x:'gain' if x=='<DUP>' else 'loss')
	df5.head()
	df6=df5[[0,1,-3,-1,-2]]
	df6['key']=df6['#CHROM'].apply(lambda x:trans(x))
	df7=df6.sort_values(['key','POS'])
	df8=df7[[0,1,2,3,4]]
	df8.head()
	df8.to_csv(lc_bed,index=None,header=None,sep='\t')

def saascnv_com(sc_vcf,sc_bed):
	df=pd.read_csv(sc_vcf,sep='\t')
	df.head()
	filter=['LOH','gain','loss']
	df1=df[df['CNV'].apply(lambda x:x in filter)]
	df1['CNV']=df1['CNV'].apply(lambda x:x.replace('LOH','loss'))
	df2=df1[list(map(lambda x:x[0]<=0.05 and x[1]<=0.01 and x[2]<=0.001,df1[['log2ratio.p.value','log2mBAF.p.value','p.value']].values))]
	df2.head()
	df2['ratio.mbaf.p.value']=list(map(lambda x:';'.join([str(x[0]),str(x[1]),str(x[2])]),df2[['log2ratio.p.value','log2mBAF.p.value','p.value']].values))
	df2.head()
	df3=df2[[0,1,2,-2,-1]]
	df3.to_csv(sc_bed,index=None,header=None,sep='\t')

def golden_com(gc_vcf,gc_bed):
	df=pd.read_csv(gc_vcf,sep='\t')
	df.head()
	df['chr']=df['chr'].apply(lambda x:"chr"+str(x))
	df['CN']=df['Copy_Number'].apply(lambda x:x+2)
	df1=df[[0,1,2,-1]]
	df2=df1[df1['CN']!=2]
	df2.head()
	df2.to_csv(gc_bed,index=None,header=None,sep='\t')	
	
def Query():
	pyshell("bedtools intersect -wa -wb -a bicseq2.bed -b gold.10X.bed > bicseq2.query.gold.10X.bed")
	pyshell("bedtools intersect -wa -wb -a cnmops.bed -b gold.10X.bed > cnmops.query.gold.10X.bed")
	pyshell("bedtools intersect -wa -wb -a freec.bed -b gold.10X.bed > freec.query.gold.10X.bed")
	pyshell("bedtools intersect -wa -wb -a lumpy.bed -b gold.10X.bed > lumpy.query.gold.10X.bed")
	pyshell("bedtools intersect -wa -wb -a saascnv.bed -b gold.10X.bed > saascnv.query.gold.10X.bed")

def bicseq2_deal(bd_vcf,be_bed):
	df=pd.read_csv(bd_vcf,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:Judge(x[0])==Judge(x[1]),df[[3,7]].values))]
	df1.head()
	df1[8]=list(map(lambda x:Score3(x[0],x[1],x[2],x[3]),df1[[1,2,5,6]].values))
	df1.head()
	df2=df1[df1[8]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3])
	df4=df2.drop_duplicates([4,5,6,7])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3]]
	df5.to_csv('bicseq2.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[4,5,6,7]]
	df6[8]=df6.loc[:,[4,5,6,7]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[8]]
	df7.to_csv(be_bed,index=None,header=None,sep='\t')

def cnmops_deal(cd_vcf,cd_bed):
	df=pd.read_csv(cd_vcf,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:Judge(x[0])==Judge(x[1]),df[[3,7]].values))]
	df1.head()
	df1[8]=list(map(lambda x:Score3(x[0],x[1],x[2],x[3]),df1[[1,2,5,6]].values))
	df1.head()
	df2=df1[df1[8]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3])
	df4=df2.drop_duplicates([4,5,6,7])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3]]
	df5.to_csv('cnmops.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[4,5,6,7]]
	df6[8]=df6.loc[:,[4,5,6,7]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[8]]
	df7.to_csv(cd_bed,index=None,header=None,sep='\t')

def freec_deal(fd_vcf,fd_bed):
	df=pd.read_csv(fd_vcf,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:Judge(x[0])==Judge(x[1]),df[[3,7]].values))]
	df1.head()
	df1[8]=list(map(lambda x:Score3(x[0],x[1],x[2],x[3]),df1[[1,2,5,6]].values))
	df1.head()
	df2=df1[df1[8]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3])
	df4=df2.drop_duplicates([4,5,6,7])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3]]
	df5.to_csv('freec.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[4,5,6,7]]
	df6[8]=df6.loc[:,[4,5,6,7]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[8]]
	df7.to_csv(fd_bed,index=None,header=None,sep='\t')

def lumpy_deal(ld_vcf,ld_bed):
	df=pd.read_csv(ld_vcf,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:x[0]==Judge(x[1]),df[[3,8]].values))]
	df1.head()
	df1[9]=list(map(lambda x:Score3(x[0],x[1],x[2],x[3]),df1[[1,2,6,7]].values))
	df1.head()
	df2=df1[df1[9]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3,4])
	df4=df2.drop_duplicates([5,6,7,8])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3,4]]
	df5.to_csv('lumpy.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[5,6,7,8]]
	df6[9]=df6.loc[:,[5,6,7,8]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[9]]
	df7.to_csv(ld_bed,index=None,header=None,sep='\t')

def saascnv_deal(sd_vcf,sd_bed):
	df=pd.read_csv(sd_vcf,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:x[0]==Judge(x[1]),df[[3,8]].values))]
	df1.head()
	df1[9]=list(map(lambda x:Score3(x[0],x[1],x[2],x[3]),df1[[1,2,6,7]].values))
	df1.head()
	df2=df1[df1[9]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3,4])
	df4=df2.drop_duplicates([5,6,7,8])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3,4]]
	df5.to_csv('saascnv.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[5,6,7,8]]
	df6[9]=df6.loc[:,[5,6,7,8]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[9]]
	df7.to_csv(sd_bed,index=None,header=None,sep='\t')

def FiveToolsOverlap():
	pyshell("multiIntersectBed -header -i bicseq2.bed cnmops.bed freec.bed lumpy.bed saascnv.bed -names bicseq2 cnmops freec lumpy saascnv > Five.overlap.bed")

def bicseq2_info(bi_vcf,bi_bed):
	df=pd.read_csv(bi_vcf,header=None,sep='\t')
	df.head()
	df[4]='bicseq2'
	df[5]=list(map(lambda x:"chr="+str(x[0])+"|"+"start="+str(x[1])+"|"+"end="+str(x[2])+"|"+"CN="+str(x[3])+"|"+"Name="+str(x[4]),df[[0,1,2,3,4]].values))
	df.head()
	df1=df[[0,1,2,5]]
	df1.to_csv(bi_bed,index=None,header=None,sep='\t')

def cnmops_info(ci_vcf,ci_bed):
	df=pd.read_csv(,header=None,sep='\t')
	df.head()
	df[4]='cnmops'
	df[5]=list(map(lambda x:"chr="+str(x[0])+"|"+"start="+str(x[1])+"|"+"end="+str(x[2])+"|"+"CN="+str(x[3])+"|"+"Name="+str(x[4]),df[[0,1,2,3,4]].values))
	df.head()
	df1=df[[0,1,2,5]]
	df1.to_csv(,index=None,header=None,sep='\t')

def freec_info(fi_vcf,fi_bed):
	df=pd.read_csv(,header=None,sep='\t')
	df.head()
	df[4]='freec'
	df[5]=list(map(lambda x:"chr="+str(x[0])+"|"+"start="+str(x[1])+"|"+"end="+str(x[2])+"|"+"CN="+str(x[3])+"|"+"Name="+str(x[4]),df[[0,1,2,3,4]].values))
	df.head()
	df1=df[[0,1,2,5]]
	df1.to_csv(,index=None,header=None,sep='\t')

def lumpy_info(li_vcf,li_bed):
	df=pd.read_csv(,header=None,sep='\t')
	df.head()
	df[5]='lumpy'
	df[6]=list(map(lambda x:"chr="+str(x[0])+"|"+"start="+str(x[1])+"|"+"end="+str(x[2])+"|"+"CNV="+str(x[3])+"|"+"SU.PE.SR="+str(x[4])+"|"+"Name="+str(x[5]),df[[0,1,2,3,4,5]].values))
	df.head()
	df1=df[[0,1,2,6]]
	df1.to_csv(,index=None,header=None,sep='\t')

def saascnv_info(si_vcf,si_bed):
	df=pd.read_csv(,header=None,sep='\t')
	df.head()
	df[5]='saascnv'
	df[6]=list(map(lambda x:"chr="+str(x[0])+"|"+"start="+str(x[1])+"|"+"end="+str(x[2])+"|"+"CNV="+str(x[3])+"|"+"ratio.mbaf.p.value="+str(x[4])+"|"+"Name="+str(x[5]),df[[0,1,2,3,4,5]].values))
	df.head()
	df1=df[[0,1,2,6]]
	df1.to_csv(,index=None,header=None,sep='\t')
	
def MergeInfoBed():
	pyshell("cat bicseq2.info.bed cnmops.info.bed freec.info.bed lumpy.info.bed saascnv.info.bed > Info.merge.bed")

def Finalout(fl_over_bed,fl_info_bed,final_bed,filter_bed):
	df1=pd.read_csv(fl_over_bed,sep='\t')
	df2=pd.read_csv(fl_info_bed,header=None,sep='\t')
	df2.columns=['chrom','start','end','info']
	df1.head()
	df2.head()
	df1['info']=list(map(lambda x:';;'.join(df2[df2['chrom']==x[0]][list(map(lambda y:((x[1]>=y[0] and x[1]<=y[1]) or (x[2]>=y[0] and x[2]<=y[1])) and y[0]!=x[2] and y[1]!=x[1],df2[df2['chrom']==x[0]][['start','end']].values))]['info']),df1[['chrom','start','end']].values))
	df1.to_csv(final_bed,index=None,sep='\t')
	df1['cnv_type']=df1['info'].apply(lambda x:Deal_info(x))	
	df1.head()
	set(df1['cnv_type'])
	df1['set_type']=df1['cnv_type'].apply(lambda x:';'.join(list(set(x.split(';')))))
	df1.head()
	set(df1['set_type'])
	df1['set_type'].value_counts()
	df2=df1[df1['set_type'].apply(lambda x:'gain' not in x or 'loss' not in x)]
	df2.head()
	df2['set_type'].value_counts()
	df3=df2[df2['set_type'].apply(lambda x:Set_type(x))]			
	df3.head()
	df3['set_type'].value_counts()
	len(df3.index)
	df3['cnv']=df3['set_type'].apply(lambda x:CNV(x))			
	df3.head()
	df3['cnv'].value_counts()
	df3.to_csv(filter_bed,index=None,header=None,sep='\t')
	
def Sort_filter(sf_out,sf_sout):
	df=pd.read_csv(sf_out,header=None,sep='\t')
	df.head()
	df[14]=df[0].apply(lambda x:trans(x))
	df1=df.sort_values([14,1])
	df1.head()
	df2=df1[list(df1.columns)[:-1]]
	df2.to_csv(sf_sout,index=None,header=None,sep='\t')
	pyshell("bedtools merge -i Our.method.origin.filter.sort.out -c 4,5,11,13,14 -o max,distinct,distinct,distinct,distinct > Our.method.origin.filter.merge.out")

def merge_filter(mer_out,final_out):
	df=pd.read_csv(mer_out,header=None,sep='\t')
	df.head()
	df[4]=df[4].apply(lambda x:','.join(list(set(x.split(',')))))
	df.head()
	df[5]=df[5].apply(lambda x:x.replace(',chr',';;chr'))
	df.head()
	df[5]=df[5].apply(lambda x:';;'.join(list(set(x.split(';;')))))
	df.head()
	df[6]=df[6].apply(lambda x:x.replace(',',';'))
	df[6]=df[6].apply(lambda x:';'.join(list(set(x.split(';')))))
	df.head()	
	df1=df[df[4].apply(lambda x:Tools(x))]		
	df1.head()
	df1[4].value_counts()
	df1.to_csv(final_out,index=None,header=None,sep='\t')
	pyshell("bedtools intersect -wa -wb -a Our.method.origin.filter.final.out -b gold.10X.bed > Our.method.origin.filter.final.query.gold.bed")

def OurMethod(query_bed,over_bed):	
	df=pd.read_csv(query_bed,header=None,sep='\t')
	df.head()
	df1=df[list(map(lambda x:x[0]==Judge(x[1]),df[[7,11]].values))]
	df1.head()
	df1[12]=list(map(lambda x:Score31(x[0],x[1],x[2],x[3]),df1[[1,2,9,10]].values))
	df1.head()
	df2=df1[df1[12]>=50.]
	df2.head()
	len(df2.index)
	df3=df2.drop_duplicates([0,1,2,3,4,5,6,7])
	df4=df2.drop_duplicates([8,9,10,11])
	len(df3.index)
	len(df4.index)
	df5=df3[[0,1,2,3,4,5,6,7]]
	df5.to_csv('our.method.found.overlap.with.gold.bed',index=None,header=None,sep='\t')
	df6=df4[[8,9,10,11]]
	df6[12]=df6.loc[:,[8,9,10,11]].apply(lambda x:'|'.join([str(i) for i in x]),axis=1)
	df7=df6[[12]]
	df7.to_csv(over_bed,index=None,header=None,sep='\t')	
	
	
def main():
	Rename()
	bicseq2('bicseq2.vcf','bicseq2.bed')
	cnmops('cnmops.vcf','cnmops.bed')
	freec('freec.vcf','freec.bed')
	lumpy('lumpy.vcf','lumpy.bed')
	saascnv('saascnv.vcf','saascnv.bed')
	bicseq2_com('bicseq2.vcf','bicseq2.origin.bed')
	cnmops_com('cnmops.vcf','cnmops.origin.bed')
	freec_com('freec.vcf','freec.origin.bed')
	lumpy_com('lumpy.vcf','lumpy.origin.bed')
	saascnv_com('saascnv.vcf','saascnv.origin.bed')
	golden_com('gold.10X.vcf','gold.10X.bed')
	Query()
	bicseq2_deal('bicseq2.query.gold.10X.bed','bicseq2.overlap.with.gold.bed')
	cnmops_deal('cnmops.query.gold.10X.bed','cnmops.overlap.with.gold.bed')
	freec_deal('freec.query.gold.10X.bed','freec.overlap.with.gold.bed')
	lumpy_deal('lumpy.query.gold.10X.bed','lumpy.overlap.with.gold.bed')
	saascnv_deal('saascnv.query.gold.10X.bed','saascnv.overlap.with.gold.bed')
	FiveToolsOverlap()
	bicseq2_info('bicseq2.bed','bicseq2.info.bed')
	cnmops_info('cnmops.bed','cnmops.info.bed')
	freec_info('freec.bed','freec.info.bed')
	lumpy_info('lumpy.bed','lumpy.info.bed')
	saascnv_info('saascnv.bed','saascnv.info.bed')
	Finalout('Five.overlap.bed','Info.merge.bed','Five.mult.tools.overlap.final.out.bed','Our.method.origin.filter.out')
	Sort_filter('Our.method.origin.filter.out','Our.method.origin.filter.sort.out')	
	merge_filter('Our.method.origin.filter.merge.out','Our.method.origin.filter.final.out')
	OurMethod('Our.method.origin.filter.final.query.gold.bed','our.method.overlap.with.gold.bed')