#!/usr/local/bin/python
import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument('SampleBED',type=str,help='Path to the mosdepth per-base BED output')
parser.add_argument('RefExomeBED',type=str,help='Path to the all exons BED file')
parser.add_argument('SNVGermlineTXT',type=str,help='Path to Clivar-generated table with pathogenic germline SNVs')
parser.add_argument('SNVSomaticTXT',type=str,help='Path to Clivar-generated table with pathogenic somatic SNVs')
parser.add_argument('CoverThreshold',type=int,help='Coverage quality threshold')
parser.add_argument('GeneName_s',nargs='+',help = 'Gene name(s)')

args = parser.parse_args()

def snv_coverage(snv,chrom_cover):
    snv = snv.reset_index(drop=True)
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

def generate_figure(plot_num,snvg,snvs,plot_coverage,exon_cover,gene,chrom,output_name,threshold):
    if plot_num==1:
        plt.figure(figsize=(30,10))
    cp1 = ['red','blue','#B0B4FF']
    plt.subplot(2,1,plot_num)
    if plot_num==1:
        plt.plot(plot_coverage.position,plot_coverage.coverage,'-',c='grey',lw=0.5)
    else:
        plt.plot([0,exon_cover.end.iloc[-1]],[threshold,threshold],'--r',alpha=0.5)
    for ind in exon_cover.index:
        buf = exon_cover.loc[ind]
        if (buf.end-buf.start)>1:
            plt.plot([buf.start,buf.end-1],[buf.coverage,buf.coverage],'-',c=cp1[2],lw=4)
            plt.plot([buf.start,buf.end-1],[0,0],'-',c=cp1[2],lw=4)
        else:
            plt.plot([buf.start],[buf.coverage],'s',c=cp1[2],ms=2.5)
            plt.plot([buf.start],[0],'s',c=cp1[2],ms=2.5)
        if buf.coverage<threshold and plot_num==2:
            if (buf.end-buf.start)>1:
                plt.plot([buf.start,buf.end-1],[buf.coverage,buf.coverage],'-',c='r',alpha=0.2,lw=7)
            else:
                plt.plot([buf.start],[buf.coverage],'s',c='r',alpha=0.2,ms=4)
    if snvg.size>0:
        plt.plot(snvg.position,np.ones((snvg.position.size)),'.',c=cp1[0],ms=3.5)
        plt.plot(snvg.position,snvg.coverage+0.1,'.',c=cp1[0],ms=3.5)
    if snvs.size>0:
        plt.plot(snvs.position,-np.ones((snvs.position.size)),'.',c=cp1[1],ms=3.5)
        plt.plot(snvs.position,snvs.coverage-0.1,'.',c=cp1[1],ms=3.5)
    plt.grid(axis='y')
    if plot_num==1:
        plt.title(gene+' coverage | '+chrom,fontdict = {'fontsize' : 17})
        plt.xticks([plot_coverage.position.iloc[0],plot_coverage.position.iloc[-1]],
                   [plot_coverage.position.iloc[0],plot_coverage.position.iloc[-1]])
    else:
        plt.title(gene+' coverage (only exons) | '+chrom,fontdict = {'fontsize' : 17})
        plt.xticks([0,exon_cover.end.iloc[-1]],
                  [plot_coverage.position.iloc[0],plot_coverage.position.iloc[-1]])
    if plot_num==1:
        legend_entries = [Line2D([0],[0],color='grey',lw=0.5),
                          Line2D([0],[0],color=cp1[2],lw=4),
                          Line2D([0],[0],marker='.',color=cp1[0],lw=0,markersize=3.5),
                          Line2D([0],[0],marker='.',color=cp1[1],lw=0,markersize=3.5)]
        legend_labels = ['sequence coverage',
                         'exon coverage',
                         'germline SNV',
                         'somatic SNV']
        plt.legend(legend_entries,legend_labels,loc='upper left')
    if plot_num==2:
        legend_entries = [Line2D([0],[0],color='red',ls='--',alpha=0.5),
                          Line2D([0],[0],color=cp1[2],lw=4),
                          Line2D([0],[0],marker='.',color=cp1[0],lw=0,markersize=3.5),
                          Line2D([0],[0],marker='.',color=cp1[1],lw=0,markersize=3.5)]
        legend_labels = ['coverage threshold',
                        'exon coverage',
                        'germline SNV',
                        'somatic SNV',
                        ]
        plt.xlim(-0.05*exon_cover.end.iloc[-1],1.05*exon_cover.end.iloc[-1])
        plt.legend(legend_entries,legend_labels,loc='upper left')
        plt.savefig(output_name,bbox_inches='tight')

        
def ReadBEDPart(BED_fname,chrom):
    BED = open(BED_fname, 'rt')
    line_num = 0
    while True:
        line = BED.readline()
        line_chr = line[:line.find('\t')]
        if line_chr==chrom:
            start = line_num
            break
        else:
            line_num +=1
    while True:
        line = BED.readline()
        line_chr = line[:line.find('\t')]
        if line_chr!=chrom:
            stop = line_num
            break
        else:
            line_num +=1
    BED.close()
    BEDPart = pd.read_csv(BED_fname,delimiter='\t',header=None,names=['chr','start','end','coverage'],skiprows=start,nrows=stop-start+1)
    return BEDPart


def ChromosomeRegionViewGeneration(exons,snvgs,snvss,sample_bed_path,gene,threshold):
    
    sample_name = sample_bed_path
    while sample_name.find('/')!=-1:
        sample_name = sample_name[sample_name.find('/')+1:]
    rows = ['' for i in range(len(exons.chr.unique()))]
    for chrom,chr_num in zip(exons.chr.unique(),range(len(exons.chr.unique()))):
        chrom_cover = ReadBEDPart(sample_bed_path,chrom)
        chrom_exons = exons[exons.chr==chrom].reset_index(drop=True)
        snvg = snvgs[snvgs.chr==chrom[3:]]
        snvs = snvss[snvss.chr==chrom[3:]]
        exon_cover = pd.DataFrame(columns=['chr','start','end','coverage','exnum'])
        plot_coverage = pd.DataFrame(columns=['position','coverage'])

        reg_starts = chrom_cover.index[chrom_cover.start<chrom_exons.start.iloc[0]][-1]
        reg_ends = chrom_cover.index[chrom_cover.end>chrom_exons.end.iloc[-1]][0]
        cover_reg = chrom_cover[(chrom_cover.index>=reg_starts) & (chrom_cover.index<=reg_ends)]
        cover_reg = cover_reg.reset_index(drop=True)
        del reg_starts,reg_ends

        for ind in chrom_exons.index:
            buf = cover_reg[(cover_reg.end>chrom_exons.start[ind])&(cover_reg.start<chrom_exons.end[ind])]
            buf.start.iloc[0] = chrom_exons.start[ind]
            buf.end.iloc[-1] = chrom_exons.end[ind]
            buf['exnum']=ind
            exon_cover = pd.concat([exon_cover,buf])
        for ind in cover_reg.index:
            buf = pd.DataFrame(data={'position':[cover_reg.start[ind],cover_reg.end[ind]-1],'coverage':[cover_reg.coverage[ind],cover_reg.coverage[ind]]})
            plot_coverage = pd.concat([plot_coverage,buf])
        exon_cover = exon_cover.reset_index(drop=True)
        plot_coverage = plot_coverage.reset_index(drop=True)
        if snvg.size>0:
            snvg = snv_coverage(snvg,chrom_cover)
            at_germ = np.sum(snvg.coverage>=threshold)/len(snvg)*100
            print('%3.0f%% of pathogenic germline SNVs for %5s in %5s are covered above threshold (%2.0fx)'%(
            at_germ,gene,chrom,threshold))
            germ_row = '%8.0f %10.0f %6.0f %6.0f %6.0f %4.0f%%'%(len(snvg),
                                                        np.median(snvg.coverage),
                                                        np.std(snvg.coverage),
                                                        np.min(snvg.coverage),
                                                        np.max(snvg.coverage),
                                                        at_germ)
        else:
            germ_row = '%8.0f %10.0f %6.0f %6.0f %6.0f %4.0f%%'%(0,0,0,0,0,0)
        if snvs.size>0:
            snvs = snv_coverage(snvs,chrom_cover)
            at_soma = np.sum(snvs.coverage>=threshold)/len(snvs)*100
            print('%3.0f%% of pathogenic  somatic SNVs for %5s in %5s are covered above threshold (%2.0fx)'%(
            at_soma,gene,chrom,threshold))
            soma_row = '%9.0f %10.0f %6.0f %6.0f %6.0f %4.0f%%'%(len(snvs),
                                                        np.median(snvs.coverage),
                                                        np.std(snvs.coverage),
                                                        np.min(snvs.coverage),
                                                        np.max(snvs.coverage),
                                                        at_soma)
        else:
            soma_row = '%9.0f %10.0f %6.0f %6.0f %6.0f %4.0f%%'%(0,0,0,0,0,0)
        
        rows[chr_num] = '%6s %6s'%(gene,chrom)+' '+germ_row+soma_row+'\n'
        output_name = 'output/'+gene+'.'+chrom+'.'+sample_name+'.pdf'
        generate_figure(1,snvg,snvs,plot_coverage,exon_cover,gene,chrom,output_name,threshold)

        first_pos = exon_cover.start.iloc[0]
        exon_cover.start -= first_pos 
        exon_cover.end -= first_pos 
        snvs.position -= first_pos
        snvg.position -= first_pos
        for num in exon_cover.exnum.unique()[:-1]:
            last = exon_cover[exon_cover.exnum==num].end.iloc[-1]
            mask = exon_cover.exnum>num
            first = exon_cover[mask].start.iloc[0]
            bufor = exon_cover[mask]
            bufor.start = bufor.start - (first-last-70)
            bufor.end = bufor.end - (first-last-70)
            exon_cover[mask] = bufor
            mask = snvs.position>=bufor.start.iloc[0]
            bufor_snv = snvs[mask]
            bufor_snv.position = bufor_snv.position - (first-last-70)
            snvs[mask] = bufor_snv
            mask = snvg.position>=bufor.start.iloc[0]
            bufor_snv = snvg[mask]
            bufor_snv.position = bufor_snv.position - (first-last-70)
            snvg[mask] = bufor_snv
        generate_figure(2,snvg,snvs,plot_coverage,exon_cover,gene,chrom,output_name,threshold)
    return rows
exon0 = pd.read_csv(args.RefExomeBED,delimiter='\t',header=None,names=['chr','start','end','gene'])
snvg0 = pd.read_csv(args.SNVGermlineTXT,sep='\t')
snvs0 = pd.read_csv(args.SNVSomaticTXT,sep='\t')

exon0=exon0.dropna()
exon0=exon0.drop_duplicates()
exon0 = exon0.reset_index(drop=True)
snvg0 = snvg0[['Gene(s)','GRCh38Chromosome','GRCh38Location']]
snvg0.columns = ['gene','chr','position']
snvs0 = snvs0[['Gene(s)','GRCh38Chromosome','GRCh38Location']]
snvs0.columns = ['gene','chr','position']
        
genes = args.GeneName_s
sample_name = args.SampleBED
while sample_name.find('/')!=-1:
    sample_name = sample_name[sample_name.find('/')+1:]
logs = []
for gene in genes:    
    exons = exon0[exon0.gene==gene][['chr','start','end']].reset_index(drop=True)
    snvgs = snvg0[snvg0.gene==gene][['chr','position']]
    snvss = snvs0[snvs0.gene==gene][['chr','position']]
    if np.any(exon0.gene==gene):
        log=ChromosomeRegionViewGeneration(exons,snvgs,snvss,args.SampleBED,gene,args.CoverThreshold)
        logs=np.append(logs,log,axis=0)
    else:
        print('Gene %s was not found in Exome reference list.'%(gene))
out_table = open('output/'+sample_name+'.geneCoverage.txt','wt')
out_table.write('geneCoverage - %s\n\n'%(sample_name))
out_table.write('Input genes: %s\n\n'%(genes))
out_table.write('Pathogenic (G - germline, S - somatic) SNV coverage:\n')
out_table.write('- \'count\' is the number of variants in a given region\n- \'AT\' is the percentage of SNV coverage above threshold (%2.0fx)\n\n'%(args.CoverThreshold))
out_table.write('%6s %6s %8s %10s %6s %6s %6s %5s %8s %10s %6s %6s %6s %5s\n'%('gene','chr','count(G)','median(G)','std(G)','min(G)','max(G)','AT(G)','count(S)','median(S)','std(S)','min(S)','max(S)','AT(S)'))
out_table.write('-----------------------------------------------------------------------------------------------------------\n')
for log in logs:
    out_table.write(log)
out_table.write('-----------------------------------------------------------------------------------------------------------\n')
out_table.close()
        

