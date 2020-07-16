#!/usr/bin/env python3
'''
SNPs and INDELs Calling modules
'''

import time
import os
import sys
import uuid
import json
import logging
import tempfile
import argparse
import subprocess as sp
from cosmos.api import Cosmos

#dbsnp='/ref/hg38/dbSNP/00-common_all.vcf.gz'
dbsnp=None
#'/home/ubuntu/prj/snpexample/ref/hg38/dbSNP/00-common_all.vcf.gz'
#dbsnp='/ref/dbsnp_146.hg38.vcf.gz'
biallelic_vcf="/ref/truthSets/population_biallelic_snvs/COLO829BL_population_biallelic.vcf.gz"

#default_truth_sets={ 'GM12878':'/local/HG002_v3.3.2_highconf_triophased.vcf.gz'}
default_truth_sets={ 'GM12878':'/local/HG002_v3.3.2_highconf_triophased.vcf.gz'}

#default_truth_sets={ 'GM12878':'/ref/truthSets/GM12878/illumina_platinum/pg2017/hg38/small_variants/NA12878/NA12878.vcf.gz',
#            'NA12878':'/ref/truthSets/GM12878/illumina_platinum/pg2017/hg38/small_variants/NA12878/NA12878.vcf.gz',
#            'HCC1187':'/ref/truthSets/HCC1187/somatic/small_sv/HCC1187_somatic_pass.vcf.gz',
#            'COLO829':'/ref/truthSets/COLO829/hg38/small_variants/COLO829.vcf.gz'
#          }

ConfidentRegions="/local/HG002_v3.3.2_highconf_noinconsistent.bed"
#ConfidentRegions="/ref/truthSets/GM12878/illumina_platinum/pg2017/hg38/small_variants/ConfidentRegions.sorted.bed"
#ConfidentRegions="/ref/wgs_calling_regions.hg38.interval_list"

popu_germline_AF='/ref/hg38/gatk_resources/af-only-gnomad.hg38.vcf.gz'
pon='/ref/hg38/gatk_resources/1000g_pon.hg38.vcf.gz'
default_ref = '/ref/hg38/GRCh38.p12.fa' 
#default_ref = '/ref/Homo_sapiens_assembly38.fasta' 

'''
---------------------------------------
        *** HELPER FUNCTIONs ***

    GATK Docker container manipulation
---------------------------------------
'''

hex_uid = lambda: uuid.uuid4().hex[:6].upper()
base_no_ext = lambda f, e: os.path.basename(f).split('.'+e)[0]

def gatkContainer():
    return '_'.join(['gatk', hex_uid()])


def runDocker(container):
    return os.system("docker run -d -t --name %s -v ~:/gatk/$HOME -d broadinstitute/gatk:4.1.2.0" % container)


def stopDocker(container):
    return r"""
    docker stop {container}; \
    docker rm {container}
    """.format(container=container)


'''
    -----------------------------------
    Files and directories Manipulations
    -----------------------------------
'''
def tmpDir(dir=os.getcwd):
    tmp_dir = 'tmp'+ hex_uid()
    while os.path.isdir(os.path.join(dir,tmp_dir)):
        tmpDir = 'tmp'+hex_uid()
    tmp_dir= os.path.join(dir,tmp_dir)
    os.system('mkdir -p %s' % tmp_dir)
    return tmp_dir


def readJSON(filename):
    try:
        with open(filename, "r") as jsonfile:
            return json.load(jsonfile)
    except Exception as e:
        logging.error("Error" % e)
    return None


'''
--------------------------------------------------------
        ***  WORKFLOW UNDERLYING FUNCTIONS ***

    Split reference genome to non-overlapping intervals
--------------------------------------------------------
'''
def get_intervals(in_interval_list, out_interval_list):
    if os.path.isfile(in_interval_list):
        return os.system(f"cp {in_interval_list} {out_interval_list}")

    intervals = in_interval_list.strip().split(',')
    with open(out_interval_list, 'w') as fp:
        for interval in intervals:
            print(interval, file=fp)


def split_intervals_byWindow(ref, out_interval_list, window_size):
    index_file = ref+'.fai'
    if not os.path.exists(index_file):
        logging.error("%s does not exist." % index_file)

    with open (index_file, 'r') as f:
        sizes = [line.strip().split() for line in f]
 
    with open (out_interval_list, 'w') as bed:
        w = window_size
        for _chr, _len, *rest in sizes:
            chr_size = int(_len)
            p1=1
            p2=p1+w-1
            while (p1 < chr_size):
                p2 = chr_size if chr_size-p2 < w/2 else p2
                print("%s:%d-%d" % (_chr, p1, p2), file=bed)
                p1=p2+1
                p2=p1+w-1
    pass

def split_intervals_byNs(ref, out_interval_list):
    default_intervals= "/ref/hg38/ACGT_intervals.txt"
    if ref == default_ref:
        logging.info(f"using {default_intervals}")
        return os.system(f"cp {default_intervals} {out_interval_list}")

    out_dir=os.path.dirname(out_interval_list)
    os.system(r"""
    gatk ScatterIntervalsByNs \
    -O {out_dir}/skip_Ns.list \
    -R {ref} \
    --MAX_TO_MERGE 1000000 \
    --OUTPUT_TYPE ACGT \
    --VERBOSITY ERROR ; \
    gatk IntervalListToBed \
    -I {out_dir}/skip_Ns.list \
    -O {out_dir}/skip_Ns.bed ; \
    awk '{{print $1, $2, $3, $3-$2}}' {out_dir}/skip_Ns.bed > {out_dir}/tmp.bed ; \
    sort -k4nr {out_dir}/tmp.bed  | cut -f 1-3 | awk '{{print $1":"$2+1"-"$3}}' > {out_interval_list}; \
    rm -f {out_dir}/skip_Ns.list {out_dir}/skip_Ns.bed {out_dir}/tmp.bed; \
    """.format(ref=ref, out_dir=out_dir, out_interval_list=out_interval_list))
    

def split_to_intervals(ref, interval_list, intervals=None, window=0):
    if intervals is not None:
        return  get_intervals(intervals, interval_list)
    window = abs(window)
    if window == 0:
        return split_intervals_byNs(ref, interval_list)
    else:
        return split_intervals_byWindow(ref, interval_list, window)

'''
    ----------------------    
    Pre-processing modules
    ---------------------- 
'''
def baseRecalibrator(ref, in_bam, out_prefix, interval, memory):
    #hapmap=     "/ref/truthSets/GM12878/gatk/hapmap_3.3.hg38.vcf.gz"
    #g1000=      "/ref/truthSets/GM12878/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    
    #hapmap = "/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    #g1000="/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

    return r"""
    gatk  --java-options -Xmx{memory} \
    BaseRecalibrator -R {ref} \
    -I {in_bam} \
    -L {interval} \
    -O {out_prefix}.table \
    --known-sites  {dbsnp} \
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, interval=interval, memory=memory, dbsnp=dbsnp)

#    --known-sites {hapmap} \
#    --known-sites {g1000}

def applyBqsrFilter(ref, in_bam, out_prefix, interval, memory):
    return r"""
    gatk --java-options -Xmx{memory} \
    ApplyBQSR -R {ref} \
    -I {in_bam} \
    -L {interval} \
    -O {out_prefix}.bam \
    --bqsr-recal-file {out_prefix}.table
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, interval=interval, memory=memory)


'''
    --------------------------------
    SNPs and INDELs calling modules
    --------------------------------
'''
def haplotypeCaller(ref, in_bam, out_prefix, sample, interval, memory):
    tmp_dir = tmpDir(os.path.dirname(out_prefix))
    tmp_prefix = os.path.join(tmp_dir, os.path.basename(out_prefix))
    user = os.environ["USER"]
    os.system(f"mkdir -p /{tmp_dir}/tmp")    

    return r"""
    gatk --java-options -Xmx{memory} \
    HaplotypeCaller -R {ref} \
    -I {in_bam} \
    -O {tmp_prefix}.g.vcf.gz \
    -bamout {out_prefix}.bamout.bam \
    -isr INTERSECTION \
    -ip 1000 \
    -L {interval} \
    -ERC GVCF \
    --native-pair-hmm-threads 1 \
    --max-alternate-alleles 3 \
    -contamination 0 \
    --tmp-dir /{tmp_dir}/tmp \
    --QUIET ; \
    gatk  SelectVariants -R {ref} \
    -V {tmp_prefix}.g.vcf.gz \
    -L {interval}\
    -O {out_prefix}.g.vcf ;  \
    bgzip -f  {out_prefix}.g.vcf ; \
    tabix -f -p vcf {out_prefix}.g.vcf.gz; \
    rm -rf {tmp_dir}
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, tmp_dir=tmp_dir, tmp_prefix=tmp_prefix, interval=interval, user=user, memory=memory)


def mutect2_tumor(ref, in_bam, out_prefix, sample, interval,memory):
    tmp_dir = tmpDir(os.path.dirname(out_prefix))
    tmp_prefix = os.path.join(tmp_dir, os.path.basename(out_prefix))

    return r"""
    gatk --java-options -Xmx{memory} \
    Mutect2 -R {ref} \
    -I {in_bam} \
    -tumor {sample} \
    -O {tmp_prefix}.vcf.gz \
    -L {interval} \
    -bamout {out_prefix}_bamout.bam \
    -isr INTERSECTION \
    -ip 1000 \
    --annotation GenotypeSummaries \
    --annotation ChromosomeCounts \
    --germline-resource  {popu_germline_AF} \
    --panel-of-normals   {pon} \
    --native-pair-hmm-threads 1 \
    --QUIET ; \
    gatk  SelectVariants -R {ref} \
    -V {tmp_prefix}.vcf.gz \
    -L {interval} \
    -O {out_prefix}.vcf ; \
    bgzip -f  {out_prefix}.vcf; \
    tabix -f -p vcf {out_prefix}.vcf.gz; \
    cp {tmp_prefix}.vcf.gz.stats  {out_prefix}.vcf.gz.stats ; \
    #rm -rf {tmp_dir}
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, sample=sample, interval=interval,
               tmp_dir=tmp_dir, tmp_prefix=tmp_prefix, popu_germline_AF=popu_germline_AF,
                pon=pon, memory=memory) 



def mutect2_tumor_n_normal(in_tbam, in_nbam, out_prefix, ref, t_sample, n_sample, interval, memory):
    tmp_dir = tmpDir(os.path.dirname(out_prefix))
    tmp_prefix = os.path.join(tmp_dir, os.path.basename(out_prefix))

    return r"""
    gatk Mutect2 --java-options -Xmx{memory} \
    -R {ref} \
    -I {in_tbam} \
    -tumor {t_sample} \
    -I {in_nbam} \
    -normal {n_sample} \
    -O {tmp_prefix}.vcf.gz \
    -L {interval} \
    -bamout {out_prefix}_bamout.bam \
    -isr INTERSECTION \
    -ip 1000 \
    --annotation GenotypeSummaries \
    --annotation ChromosomeCounts \
    --germline-resource  {popu_germline_AF} \
    --panel-of-normals   {pon} \
    --native-pair-hmm-threads 1 \
    --QUIET ; \
    gatk  SelectVariants -R {ref} \
    -V {tmp_prefix}.vcf.gz \
    -L {interval} \
    -O {out_prefix}.vcf ; \
    bgzip -f  {out_prefix}.vcf; \
    tabix -f -p vcf {out_prefix}.vcf.gz; \
    cp {tmp_prefix}.vcf.gz.stats  {out_prefix}.vcf.gz.stats ; \
    #rm -rf {tmp_dir}
    """.format(ref=ref, in_tbam=in_tbam, in_nbam=in_nbam, out_prefix=out_prefix,
                t_sample=t_sample, n_sample=n_sample, interval=interval, 
                tmp_dir=tmp_dir, tmp_prefix=tmp_prefix, pon=pon,
                popu_germline_AF=popu_germline_AF, memory=memory)



'''
    -------------------------------------------------------------------------
     Post-precessing modules, includs genotyping and calls filtering modules 
    -------------------------------------------------------------------------
'''
def genotypeGvcf(in_gvcf, out_vcf, ref, interval,  dbsnp):

    dbsnp_option = '' if ref != default_ref else '--dbsnp ' + dbsnp 
    return r"""
    gatk GenotypeGVCFs -R {ref}  {dbsnp_option}\
    -V {in_gvcf} \
    -L {interval} \
    -O {out_vcf} \
    -OVI -OVM
    """.format(ref=ref, in_gvcf=in_gvcf, out_vcf=out_vcf, interval=interval,
                dbsnp_option=dbsnp_option)



def generateContaminationTable(in_tumor_bam, in_common_biallelic_vcf, out_prefix, interval):
    return r"""
    sed 's/\(:\|-\)/\t/g' {interval} > {interval}.bed; \
    gatk GetPileupSummaries \
    -I {in_bam} \
    -L {interval}.bed \
    -V {biallelic_vcf} \
    -O {out_prefix}_pileups.table ; \
    gatk CalculateContamination \
    -I {out_prefix}_pileups.table \
    -O {out_prefix}_contamination.table ; \
    rm -f {interval}.bed
    """.format(in_bam=in_tumor_bam, biallelic_vcf=in_common_biallelic_vcf,
                out_prefix=out_prefix,interval=interval)



def filterMutectCalls(ref, in_vcf, out_vcf):
    return f"gatk FilterMutectCalls -R {ref} -V {in_vcf} -O {out_vcf}"


def variantRecalibrator(ref, in_vcf, out_prefix):

 #   hapmap=     "/ref/truthSets/GM12878/gatk/hapmap_3.3.hg38.vcf.gz"
 #   g1000_omni= "/ref/truthSets/GM12878/gatk/1000G_omni2.5.hg38.vcf.gz"
 #   g1000=      "/ref/truthSets/GM12878/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

#    --resource:hapmap,known=false,training=true,truth=true,prior=15.0   {hapmap} \
#    --resource:omni,known=false,training=true,truth=false,prior=12.0    {1000G_omni} \
#    --resource:1000G,known=false,training=true,truth=false,prior=10.0   {1000G} \


    return r"""
    gatk VariantRecalibrator \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0    {dbsnp} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  \
    -mode BOTH \
    -R {ref} \
    -V {inVCF} \ 
    -O {out_prefix}.recal \
    --tranches-file {out_prefix}.tranches \
    -rscript-file  {out_prefix}.plots.R
    """.format(ref=ref, inVCF=in_vcf, out_prefix=out_prefix)


def applyVqsrFilter(ref, in_vcf, out_vcf, out_prefix):
    outVCF = out_vcf if out_vcf.endswith('.gz') else out_vcf+'.gz'

    return r"""
    gatk ApplyVQSR \ 
    -V {inVCF} \
    -O {outVCF} \
    -truth-sensitivity-filter-level 99.0 \
    --tranches-file {out_prefix}.tranches \
    --recal-file  {out_prefix}.recal \
    -mode BOTH
    """.format(ref=ref, inVCF=in_vcf, out_prefix=out_prefix, **locals())


'''
    ----------------
    VCF manipulation
    ----------------
'''
def mergeVCFs(in_vcf_list, out_vcf, memory):
    out_vcf = out_vcf if out_vcf.endswith('.gz') else out_vcf+'.gz'

    tmp_file = os.path.join(os.path.dirname(out_vcf),'%s_vcf_%s.list' % (base_no_ext(out_vcf, 'vcf'), hex_uid()) )
    with open(tmp_file, 'w') as fp:
        for split_vcf in in_vcf_list:
            print(split_vcf, file=fp)

    return f"picard  MergeVcfs I={tmp_file} O={out_vcf}"


def mergeStats(in_vcf_list, out_stats):
    out_dir = os.path.dirname(out_stats)
    in_stats = ' '.join(['-stats '+ x.strip()+'.stats' for x in in_vcf_list])
    return r"""gatk MergeMutectStats {in_stats} -O {out_stats}""".format(in_stats=in_stats, out_stats=out_stats, out_dir=out_dir)


def gvcf2vcf(ref, in_gvcf, out_vcf):
    out_vcf = out_vcf if out_vcf.endswith('.gz') else out_vcf+'.gz'
    return r"""
    bcftools convert --gvcf2vcf --fasta-ref {ref} -O z -o {out_vcf} {in_gvcf}
    """.format(in_gvcf=in_gvcf, out_vcf=out_vcf, ref=ref)


def remove_split_files(split_dir):
    os.system("rm -rf %s" % split_dir)



'''
    ---------------------------------
    Evaluation of SNP and INDEL calls 
    ---------------------------------
'''
def concordance(ref, in_truthVCF, in_sampleVCF, out_dir, intervals, threads):
    tuth_label      = base_no_ext(in_truthVCF,'vcf')
    sample_label    = base_no_ext(in_sampleVCF,'vcf')
    out_prefix = '%s/%s_%s' % (out_dir, sample_label, tuth_label)
    summary = "%s_summary.tsv" %( out_prefix)
    outVCF = "%s_tpfp.vcf" % (out_prefix)
    outVCF2 = "%s_tpfn.vcf" % (out_prefix)
    log = "%s.log" % (out_prefix) 
    intervals_opt = "-L %s -L %s -isr INTERSECTION " %(ConfidentRegions, intervals) if ref == default_ref else ""  
    return r"""
    gatk Concordance -R {ref}  {intervals_opt} \
    -eval {inVCF} \
    -truth {truthVCF} \
    --summary {summary}  \
    -tpfp {outVCF} \
    -tpfn {outVCF2}  &> {log}; \
    bgzip -f -@{threads} {outVCF}; \
    tabix -f -p vcf {outVCF}.gz ; \
    bgzip -f -@{threads} {outVCF2}; \
    tabix -f -p vcf {outVCF2}.gz; \
    """.format(ref=ref, inVCF=in_sampleVCF, truthVCF=in_truthVCF, summary=summary, outVCF=outVCF, outVCF2=outVCF2, intervals_opt=intervals_opt, threads=threads, log=log)




'''
----------------------------------------
        *** WORKFLOW STEPS ***
----------------------------------------
'''
''' STEP 0 ''' 
def split_to_intervals2(workflow, ref, interval_list, intervals=None, window=0, dependencies=[]):
    if intervals is not None:
        split_task = workflow.add_task(func=get_intervals,
                                        params={'in_interval_list':intervals,
                                                'out_interval_list':interval_list},
                                        core_req=1,
                                        uid="get_intervals")
        return [split_task]
 
    window = abs(window)
    if window == 0:
        split_task = workflow.add_task(func=split_intervals_byNs,
                                        params={'ref': ref, 
                                                'out_interval_list':interval_list},
                                        core_req=1,
                                        uid="split_to_intervals_" +os.path.basename(ref))
    else:
        split_task = workflow.add_task(func=split_intervals_byWindow,
                                        params={'ref': ref, 
                                                'out_interval_list':interval_list,
                                                'window_size':window},
                                        core_req=1,
                                        uid="split_to_intervals_" +os.path.basename(ref))
    return [split_task]

''' STEP 1 '''
def base_recalibration(workflow, ref, bam, out_prefix, interval, dependencies=[]):
    recal_task = workflow.add_task(func=baseRecalibrator,
                                params={'ref':ref,
                                        'in_bam':bam,
                                        'out_prefix':out_prefix,
                                        'interval':interval,
                                        'memory':'4096m'},
                                    parents=dependencies,
                                    core_req=2,
                                    max_attempts=2,
                                    mem_req=4096+256,
                                    uid='baseRecalibrator_'+os.path.basename(out_prefix))
    apply_bqsr_task = workflow.add_task(func=applyBqsrFilter,
                                    params={'ref':ref,
                                            'in_bam':bam,
                                            'out_prefix':out_prefix,
                                            'interval':interval,
                                            'memory':'4096m'},
                                    parents=[recal_task],
                                    core_req=2,
                                    max_attempts=2,
                                    mem_req=4096+256,\
                                    uid='ApplyBQSR_'+os.path.basename(out_prefix))
    return apply_bqsr_task

    

''' STEP 2 '''  
''' Alternative A '''
def germline_snp_calling(workflow, dependencies, ref, bam, sample, interval, out_prefix, *ignore):
    snp_calling_task = workflow.add_task(func=haplotypeCaller,
                                        params={'ref':ref,
                                                'in_bam':bam,
                                                'sample':sample,
                                                'interval':interval,
                                                'out_prefix':out_prefix,
                                                'memory':'4096m'},
                                        parents=dependencies,
                                        core_req=2,
                                        mem_req=4096+256,
                                        uid='HaplotypeCaller_'+os.path.basename(out_prefix))
    
    gvcf = out_prefix+'.g.vcf.gz'
    vcf  = out_prefix+'.vcf.gz'
    genotype_task =  workflow.add_task(func=genotypeGvcf,
                                       params={ 'ref':ref,
                                                'in_gvcf':gvcf,
                                                'out_vcf':vcf,
                                                'interval':interval,
                                                'dbsnp':dbsnp},
                                        parents=[snp_calling_task],
                                        core_req=2,
                                        uid='Genotype_'+os.path.basename(out_prefix))                                               
    return genotype_task

''' Alternative B '''   
def somatic_snp_calling(workflow, dependencies, ref, bam, sample, interval, out_prefix, *ignore):
    snp_calling_task = workflow.add_task(func=mutect2_tumor,
                                        params={'ref':ref,
                                                'in_bam':bam,
                                                'sample':sample,
                                                'interval':interval,
                                                'out_prefix':out_prefix,
                                                'memory':'4096m'},
                                        parents=dependencies,
                                        core_req=2,
                                        mem_req=4096+256,
                                        uid='Mutect2_tumor_'+os.path.basename(out_prefix))

    return snp_calling_task
   
 

''' Alternative C '''
def somatic_snp_calling2(workflow, dependencies, ref, bam, sample, interval, out_prefix, n_bam, n_sample):
    n_bam=os.path.abspath(n_bam)
    snp_calling_task = workflow.add_task(func=mutect2_tumor_n_normal,
                                        params={'ref':ref,
                                                'in_tbam':bam,
                                                't_sample':sample,
                                                'in_nbam':n_bam,
                                                'n_sample':n_sample,
                                                'interval':interval,
                                                'out_prefix':out_prefix,
                                                'memory':'4096m'},
                                        parents=dependencies,
                                        core_req=2,
                                        mem_req=4096+256,
                                        uid='Mutect2_tumor_N_normal_'+os.path.basename(out_prefix))
    return snp_calling_task




''' STEP 3 '''
def get_contamination_table(workflow, dependencies, in_bam, in_vcf, in_common_biallelic_vcf, intervals):
    out_prefix = in_vcf.split('.vcf')[0]
    
    pileup_task = workflow.add_task(func=generateContaminationTable,
                                    params={"in_tumor_bam":in_bam,
                                            "in_common_biallelic_vcf":in_common_biallelic_vcf,
                                            "out_prefix":out_prefix,
                                            "interval":intervals},
                                    parents=dependencies,
                                    core_req=1,
                                    uid='generateContaminationTable_'+os.path.basename(in_bam))
    return [pileup_task]


''' STEP 4 '''
def merge_vcfs(workflow,  dependencies, in_vcf_list, out_vcf):
    merge_task = workflow.add_task(func=mergeVCFs,
                                   params={'in_vcf_list':in_vcf_list,
                                           'out_vcf':out_vcf,
                                            'memory':'4096m'},
                                    parents=dependencies,
                                    core_req=1,
                                    uid='Merge_VCFs_to_'+os.path.basename(out_vcf))
    return [merge_task]


def merge_stats(workflow,  dependencies, in_vcf_list, out_stats):
    merge_task = workflow.add_task(func=mergeStats,
                                   params={'in_vcf_list':in_vcf_list,
                                           'out_stats':out_stats},
                                    parents=dependencies,
                                    core_req=4,
                                    uid='Merge_stats_to_'+os.path.basename(out_stats))
    return [merge_task]


def filter_mutect_calls(workflow, dependencies, ref, in_vcf):
    out_prefix = in_vcf.split('.vcf')[0]
    out_vcf = out_prefix +'_flt.vcf.gz'
    filter_task = workflow.add_task(func=filterMutectCalls,
                                    params={"ref":ref,
                                            "in_vcf":in_vcf,
                                            'out_vcf':out_vcf},
                                    parents=dependencies,
                                    core_req=1,
                                    uid='filterMutectCalls_'+ os.path.basename(in_vcf))

    return [filter_task]


def clean_split_vcfs(workflow,  dependencies, split_dir):
    clean_task = workflow.add_task(func=remove_split_files,
                                    params={'split_dir':split_dir},
                                    parents=dependencies,
                                    core_req=1,
                                    uid='remove_split_vcfs')
    return [clean_task]
    


''' STEP 5 '''
''' Optional, when truth set is provided,
    calculate concordance with raw calls 
'''
def eval_calls(workflow, dependencies, ref, in_truthVCF, in_vcf, out_dir, intervals):
    tuth_label =os.path.basename(in_truthVCF).split('.vcf')[0]
    sample_label =os.path.basename(in_vcf).split('.vcf')[0] 
    out_prefix = '%s/%s_%s' % (out_dir, sample_label, tuth_label)
    calls_eval_task = workflow.add_task(func=concordance,
                                        params={'ref':ref,
                                                'in_truthVCF':in_truthVCF,
                                                'in_sampleVCF':in_vcf,
                                                'out_dir':out_dir,
                                                'intervals':intervals,
                                                'threads':4},
                                        parents=dependencies,
                                        core_req=4,
                                        uid='_'.join([sample_label,tuth_label]))
    return calls_eval_task



'''
----------------------------------------
        *** WORKFLOW ***
----------------------------------------
'''

def snp_calling(workflow, ref, bam, out_dir, merged_dir, intervals, truthset=None, sample=None, somatic=False,  normal_bam=None, normal_sample=None, parents=[]):

    gatk_tasks  = []
    vcf_list    = []

    _gatk_workflow  = germline_snp_calling if not somatic else somatic_snp_calling if normal_bam is None else somatic_snp_calling2
    interval_list=os.path.join(out_dir,'ACGT_intervals.list') 

    for interval in intervals:
        # base recalibration
        bqsr_tasks  = []
        bqsr_bam    = "{}/{}_bqsr_{}".format(out_dir, os.path.basename(bam).split('.bam')[0], interval)
        bqsr_tasks.append( base_recalibration(workflow, ref, bam, bqsr_bam, interval, parents))
        bqsr_bam += '.bam'
        bqsr_normal_bam = None
        if normal_bam is not None:
            bqsr_normal_bam = "{}/{}_bqsr_{}".format(out_dir, os.path.basename(normal_bam).split('.bam')[0], interval)
            bqsr_tasks.append( base_recalibration(workflow, ref, normal_bam, bqsr_normal_bam, interval, parents))
            bqsr_normal_bam +='.bam'

        # gatk
        out_prefix =  os.path.join(out_dir, os.path.basename(bam).split('.bam')[0])
        out_prefix += '_'+interval if normal_bam is None else '_%s_%s' % (os.path.basename(normal_bam).split('.bam')[0], interval)
        gatk_task =  _gatk_workflow(workflow, bqsr_tasks, ref, bqsr_bam, sample, \
                                    interval, out_prefix, bqsr_normal_bam, normal_sample)

        vcf_file = out_prefix +'.vcf.gz'
        vcf_list.append(vcf_file)
        gatk_tasks.append(gatk_task)

    # merge and filter
    out_vcf = os.path.join(merged_dir, os.path.basename(bam).split('.bam')[0])
    out_vcf += '.vcf.gz' if normal_bam is None else '_%s' % (os.path.basename(normal_bam).split('.bam')[0] +'.vcf.gz')
    out_stats = out_vcf + '.stats'
    merge_vcfs_task = merge_vcfs(workflow, gatk_tasks, vcf_list, out_vcf)

    if somatic:
        merge_stats_task = merge_stats(workflow, gatk_tasks, vcf_list, out_stats)
        filter_task =   filter_mutect_calls(workflow, merge_vcfs_task+merge_stats_task, ref, out_vcf)
        merge_vcfs_task = filter_task
        out_vcf = out_vcf.split('.vcf')[0] +'_flt.vcf.gz'    
    
    if truthset is not None:
        if not os.path.exists(truthset):
            logging.error("%s does not exist." % truthset)
        concordance_task = eval_calls(workflow, merge_vcfs_task, ref, truthset, out_vcf, merged_dir, interval_list)


    return merge_vcfs_task


'''
----------------------------------------
    *** Run workflow per input  ***
----------------------------------------
'''

def snps_n_indels_workflow(workflow, ref, input_config, out_dir, intervals=None, window=0, clean=False, parents=[]):    
    # setup log config
    logging.basicConfig(
                handlers=[ logging.FileHandler(os.path.join(out_dir, 'workflow.log')),
                           logging.StreamHandler()],
                format='%(levelname)s: %(asctime)s  %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]',
                level=logging.DEBUG)

    # check reference
    if not os.path.exists(ref):
        logging.error("%s does not exist." % ref)
        sys.exit(1)
    ref=os.path.abspath(ref)

    # create output directory
    split_dir = out_dir+'/split_vcfs'
    merged_dir= out_dir+'/small_variants'
    os.system('mkdir -p %s' % split_dir)
    os.system('mkdir -p %s' % merged_dir)

    # create / get intervals
    interval_list=os.path.join(split_dir,'ACGT_intervals.list')
    split_to_intervals(ref, interval_list, intervals, window)

    #interval_task = split_to_intervals(workflow, ref, interval_list, intervals, window, parents)
    with open( interval_list ) as fp:
        intervals = fp.readlines()

    intervals = [interval.strip() for interval in intervals] 

    # RUN SNP calling workflow per input
    value   = lambda d, key, x:  x if key not in d else d[key]
    inputs  = readJSON(input_config)
    germline_list   = value(inputs, 'germline', [])
    somatic_list    = value(inputs, 'somatic', [])
    truthset=None

    snp_calling_tasks = []
    for record in germline_list:
        truthset = value(record,'truthSet',None)
        if truthset and truthset.upper() in default_truth_sets:
            truthset =default_truth_sets[truthset.upper()]
        if truthset is not None: truthset = os.path.abspath(truthset)
        for bam in record['bams']:
            snp_calling_task = snp_calling(workflow, ref, bam, split_dir, merged_dir, intervals, truthset, parents=parents)
            snp_calling_tasks.append(snp_calling_task)

    for record in somatic_list:
        bam     = value(record, 'bam',      None)
        sample  = value(record, 'sample',   None)
        nbam    = value(record, 'nbam',     None)
        nsample = value(record, 'nsample',  None)
        truthset= value(record, 'truthSet', None)

        if bam == None or sample == None:
            logging.warning("A bam file and a sample name are required. Skipping this input $s." %record)
            continue
        somatic = True
        if nbam != nsample and (nbam == None or nsample == None):
            logging.warning("Normal mode requires a bam file and a normal sample name. Skipping this input $s." %record)
            continue
        normal = False if nbam is None else True
        if truthset and truthset.upper() in default_truth_sets:
            truthset =default_truth_sets[truthset.upper()]

        snp_calling_task = snp_calling(workflow, ref, bam, split_dir, merged_dir, intervals, truthset, sample, somatic, nbam, nsample, parents=parents)
        snp_calling_tasks.append(snp_calling_task)
        
    if clean:
        clean_split_vcfs(workflow, snp_calling_tasks, split_dir)
 
    return


'''
----------------------------------------
            *** MAIN ***
----------------------------------------
'''

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--in_config', required=True)
    parser.add_argument('-O', '--out_dir', required=True)
    parser.add_argument('-R', '--ref', required=False, default=default_ref)
    parser.add_argument('-L', '--intervals', required=False, default=None)
    parser.add_argument('-W', '--window', required=False, default='0', type=int, 
                        help='Split reference  genome to intervals using the provided window size. \
                              Default behaviour is to split to ACGT intervals')
    parser.add_argument('-clean','--clean',  required=False, action='store_true', help='Removes split vcf files')
    parser.add_argument('-vqsr', required=False, action='store_true')
    parser.add_argument('-t','--max_cores', required=False, default=16, type=int)
    parser.add_argument('-dbsnp',required=True,"DB SNP known variants for base recalibration.")
    parser.add_argument('-drm', required=False, default="slurm")
    args = parser.parse_args()

    global dbsnp # Sorry... retrofitting someone else's code, a bigger refactor is needed. 
    dbsnp=args.dbsnp
    
    # create output directory
    os.system('mkdir -p %s' % args.out_dir)
    
    # init cosmos
    cosmos = Cosmos('sqlite:///%s/sqlite.db' % args.out_dir,
                    default_drm=args.drm, default_max_attempts=2)
    cosmos.initdb() 

    # workflow
    workflow = cosmos.start('snp_calls', restart=False, skip_confirm=True)
    snps_n_indels_workflow(workflow, args.ref, args.in_config, args.out_dir, args.intervals, args.window, args.clean) 
    workflow.run(set_successful=True, dry=False, max_cores=args.max_cores)
    sys.exit(0 if workflow.successful else 1) 
