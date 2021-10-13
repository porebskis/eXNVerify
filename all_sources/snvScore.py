#!/usr/local/bin/python
import argparse
import os
import sys
import pandas as pd
import numpy as np
import time
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(prog='snvScore')
parser.add_argument('SampleBED',type=str,help='Path to the mosdepth per-base BED output')
parser.add_argument('SNVGermlineTXT',type=str,help='Path to Clivar-generated table with pathogenic germline SNVs')
parser.add_argument('SNVSomaticTXT',type=str,help='Path to Clivar-generated table with pathogenic somatic SNVs')
parser.add_argument('Threshold',type=int,nargs='?',help='SNV coverage quality threshold (optional, positive)',default=0)
args = parser.parse_args()

sample_name = args.SampleBED
while sample_name.find('/')!=-1:
    sample_name = sample_name[sample_name.find('/')+1:]

def snv_coverage(snv,chrom_cover):
    snv = snv.dropna()
    snv['coverage']=0.0
    snv=snv.drop_duplicates()
    snv = snv.reset_index(drop=True)
    cover_reg = chrom_cover[(chrom_cover.end>snv.position.iloc[0]) & (chrom_cover.start<=snv.position.iloc[-1])]
    cover_reg = cover_reg.reset_index(drop=True)
    for ind in snv.index:
        buf = cover_reg[(cover_reg.end>snv.position[ind]) & (cover_reg.start<=snv.position[ind])]
        snv.coverage[ind] = buf.coverage
    return snv

def CatchChromoRegs(BED_fname,chrom_names):
    BED = open(BED_fname, 'rt')
#     chrom_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
#                    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
#                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
#                    'chrX', 'chrY','chrM']
    chrom_start_pos = np.zeros(len(chrom_names)+1,dtype='int32')
    line_num = 0
    for chrom,i in zip(chrom_names,np.arange(len(chrom_names))):
        pos_catched = False
        while not pos_catched:
            line = BED.readline()
            line = line[:line.find('\t')]
            if line == chrom:
                pos_catched = True
                chrom_start_pos[i] = line_num
            line_num+=1
    while line =='chrM':
        line = BED.readline()
        line = line[:line.find('\t')]
        line_num+=1
    chrom_start_pos[-1]=line_num-1
    return chrom_start_pos   

def ExecuteClinicalCoverageDepthCalc(chrom_names,SNVG,SNVS,SampleBED):
    snv_cov = pd.DataFrame(columns=['chr','position','coverage','type'])
    all_cov = np.array([])
    # start = time.time()
    res = CatchChromoRegs(SampleBED,chrom_names)
    rows = ['' for i in range(24)]
    for chrom,chr_num in zip(chrom_names[:-1],np.arange(24)):
    # for chrom,chr_num in zip(chrom_names[:3],np.arange(3)):
        chrom_cover = pd.read_csv(SampleBED,delimiter='\t',header=None,names=['chr','start','end','coverage'],skiprows=res[chr_num],nrows=res[chr_num+1]-res[chr_num])
        all_cov = np.append(all_cov,chrom_cover.coverage.values,axis=0)
        snvg_part = SNVG[SNVG.chr==chrom]
        snvs_part = SNVS[SNVS.chr==chrom]
        if snvg_part.size>0:
            snvg_part = snv_coverage(snvg_part,chrom_cover)
            snvg_part['type'] = 'germline'
            snv_cov=pd.concat([snv_cov,snvg_part])
            germ_row = '%8.0f %10.0f %6.0f %6.0f %6.0f '%(len(snvg_part),
                                                            np.median(snvg_part.coverage),
                                                            np.std(snvg_part.coverage),
                                                            np.min(snvg_part.coverage),
                                                            np.max(snvg_part.coverage))
        else:
            germ_row = '%8.0f %10.0f %6.0f %6.0f %6.0f '%(0,0,0,0,0)
        if snvs_part.size>0:
            snvs_part=snv_coverage(snvs_part,chrom_cover)
            snvs_part['type'] = 'somatic'
            snv_cov=pd.concat([snv_cov,snvs_part])
            soma_row = '%8.0f %10.0f %6.0f %6.0f %6.0f'%(len(snvs_part),
                                                            np.median(snvs_part.coverage),
                                                            np.std(snvs_part.coverage),
                                                            np.min(snvs_part.coverage),
                                                            np.max(snvs_part.coverage))
        else:
            soma_row = '%8.0f %10.0f %6.0f %6.0f %6.0f'%(0,0,0,0,0)
        rows[chr_num] = '%6s'%(chrom)+' '+germ_row+soma_row
    # end=time.time()
    germ_cov = snv_cov[snv_cov.type=='germline']
    soma_cov = snv_cov[snv_cov.type=='somatic']
    above_thres = [np.sum(germ_cov.coverage>=args.Threshold),np.sum(soma_cov.coverage>=args.Threshold)]
    out_table = open('output/'+sample_name+'.snvScore.txt','w')
    out_table.write('SNV coverage report - %s\n\n'%(sample_name))
    if args.Threshold>0:
        out_table.write('%2.0f%% of all pathogenic germline SNVs and %2.0f%% of all pathogenic somatic SNVs are covered above threshold (%2.0f)\n\n'%(above_thres[0]/len(germ_cov)*100,above_thres[1]/len(soma_cov)*100,args.Threshold))
        print('%2.0f%% of all pathogenic germline SNVs and %2.0f%% of all pathogenic somatic SNVs are covered above threshold (%2.0f)\n\n'%(above_thres[0]/len(germ_cov)*100,above_thres[1]/len(soma_cov)*100,args.Threshold))
    out_table.write('Whole genome coverage:\n')
    out_table.write('%6s %4s %4s %12s %12s %4s %7s\n'%('median',
                                                       'mean',
                                                       'std',
                                                       '1st quartile',
                                                       '3rd quartile',
                                                       'min',
                                                       'max'))
    out_table.write('%6.0f %4.0f %4.0f %12.0f %12.0f %4.0f %7.0f\n\n'%(np.median(all_cov),
                                                   np.mean(all_cov),
                                                   np.std(all_cov),
                                                   np.quantile(all_cov,0.25),
                                                   np.quantile(all_cov,0.75),
                                                   np.min(all_cov),
                                                   np.max(all_cov)))
    out_table.write('Pathogenic (G - germline, S - somatic) SNV coverage:\n')
    out_table.write('(\'count\' is the number of variants in a given region)\n\n')
    out_table.write('%6s %8s %10s %6s %6s %6s %8s %10s %6s %6s %6s\n'%('region','count(G)','median(G)','std(G)','min(G)','max(G)','count(S)','median(S)','std(S)','min(S)','max(S)'))
    out_table.write('----------------------------------------------------------------------------------------\n')
    out_table.write('%6s %8.0f %10.0f %6.0f %6.0f %6.0f %8.0f %10.0f %6.0f %6.0f %6.0f\n'%('ALL',len(germ_cov),np.median(germ_cov.coverage),np.std(germ_cov.coverage),np.min(germ_cov.coverage),np.max(germ_cov.coverage),len(soma_cov),np.median(soma_cov.coverage),np.std(soma_cov.coverage),np.min(soma_cov.coverage),np.max(soma_cov.coverage)))
    out_table.write('----------------------------------------------------------------------------------------\n')
    for row in rows:
        out_table.write(row+'\n')
    out_table.write('----------------------------------------------------------------------------------------\n')
    out_table.close()
    # print('Elased time is %3.2f sec.'%(end-start))

chrom_names = ['chr1','chr2','chr3','chr4','chr5',
               'chr6','chr7','chr8','chr9','chr10',
               'chr11','chr12','chr13','chr14','chr15',
               'chr16','chr17','chr18','chr19','chr20',
               'chr21','chr22','chrX','chrY','chrM']

SNVG = pd.read_csv(args.SNVGermlineTXT,sep='\t')[['GRCh38Chromosome','GRCh38Location']]
SNVG.columns = ['chr','position']
SNVS = pd.read_csv(args.SNVSomaticTXT,sep='\t')[['GRCh38Chromosome','GRCh38Location']]
SNVS.columns = ['chr','position']
# SNVG.chr.iloc[SNVG.chr=='MT'] = 'M'
SNVG.chr = 'chr'+SNVG.chr
# SNVS.chr.iloc[SNVS.chr=='MT'] = 'M'
SNVS.chr = 'chr'+SNVS.chr

ExecuteClinicalCoverageDepthCalc(chrom_names,SNVG,SNVS,args.SampleBED)



