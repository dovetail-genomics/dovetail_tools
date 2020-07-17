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

#dbsnp=None

'''
    -----------------------------------
    Files and directories Manipulations
    -----------------------------------
'''

hex_uid = lambda: uuid.uuid4().hex[:6].upper()
base_no_ext = lambda f, e: os.path.basename(f).split('.'+e)[0]


def tmpDir(dir=os.getcwd):
    tmp_dir = 'tmp'+ hex_uid()
    while os.path.isdir(os.path.join(dir,tmp_dir)):
        tmpDir = 'tmp'+hex_uid()
    tmp_dir= os.path.join(dir,tmp_dir)
    os.system('mkdir -p %s' % tmp_dir)
    return tmp_dir



'''
    ----------------------    
    Pre-processing modules
    ---------------------- 
'''
def baseRecalibrator(ref, in_bam, out_prefix, interval, memory):

    return r"""
    gatk  --java-options -Xmx{memory} \
    BaseRecalibrator -R {ref} \
    -I {in_bam} \
    -L {interval} \
    -O {out_prefix}.table \
    --known-sites  {dbsnp} \
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, interval=interval, memory=memory, dbsnp=dbsnp)


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
    print(f"tmp_dir = {tmp_dir}")
    tmp_prefix = os.path.join(tmp_dir, os.path.basename(out_prefix))
    print(f"tmp_prefix = {tmp_prefix}")
    user = os.environ["USER"]
    os.system(f"mkdir -p {tmp_dir}/tmp")    

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
    --tmp-dir {tmp_dir}/tmp \
    --QUIET ; \
    gatk  SelectVariants -R {ref} \
    -V {tmp_prefix}.g.vcf.gz \
    -L {interval}\
    -O {out_prefix}.g.vcf ;  \
    bgzip -f  {out_prefix}.g.vcf ; \
    tabix -f -p vcf {out_prefix}.g.vcf.gz; \
    rm -rf {tmp_dir}
    """.format(ref=ref, in_bam=in_bam, out_prefix=out_prefix, tmp_dir=tmp_dir, tmp_prefix=tmp_prefix, interval=interval, user=user, memory=memory)


'''
    -------------------------------------------------------------------------
     Post-precessing modules, includs genotyping and calls filtering modules 
    -------------------------------------------------------------------------
'''
def genotypeGvcf(in_gvcf, out_vcf, ref, interval,  dbsnp):

    dbsnp_option = f' --dbsnp {dbsnp}'  
    return r"""
    gatk GenotypeGVCFs -R {ref}  {dbsnp_option}\
    -V {in_gvcf} \
    -L {interval} \
    -O {out_vcf} \
    -OVI -OVM
    """.format(ref=ref, in_gvcf=in_gvcf, out_vcf=out_vcf, interval=interval,
                dbsnp_option=dbsnp_option)




def variantRecalibrator(ref, in_vcf, out_prefix):

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

def remove_split_files(split_dir):
    os.system("rm -rf %s" % split_dir)




'''
----------------------------------------
        *** WORKFLOW STEPS ***
----------------------------------------
'''

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




def clean_split_vcfs(workflow,  dependencies, split_dir):
    clean_task = workflow.add_task(func=remove_split_files,
                                    params={'split_dir':split_dir},
                                    parents=dependencies,
                                    core_req=1,
                                    uid='remove_split_vcfs')
    return [clean_task]
    


'''
----------------------------------------
        *** WORKFLOW ***
----------------------------------------
'''

def snp_calling(workflow, ref, bam, out_dir, merged_dir, intervals, parents=[]):

    gatk_tasks  = []
    vcf_list    = []

    normal_bam = None
    sample = None
    normal_sample = None
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
        gatk_task =  germline_snp_calling(workflow, bqsr_tasks, ref, bqsr_bam, sample, \
                                    interval, out_prefix, bqsr_normal_bam, normal_sample)

        vcf_file = out_prefix +'.vcf.gz'
        vcf_list.append(vcf_file)
        gatk_tasks.append(gatk_task)

    # merge and filter
    out_vcf = os.path.join(merged_dir, os.path.basename(bam).split('.bam')[0])
    out_vcf += '.vcf.gz' if normal_bam is None else '_%s' % (os.path.basename(normal_bam).split('.bam')[0] +'.vcf.gz')
    out_stats = out_vcf + '.stats'
    merge_vcfs_task = merge_vcfs(workflow, gatk_tasks, vcf_list, out_vcf)

    return merge_vcfs_task


'''
----------------------------------------
    *** Run workflow per input  ***
----------------------------------------
'''

def snps_n_indels_workflow(workflow, ref, bam, out_dir, intervals_file=None, clean=False, parents=[]):    
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


    #interval_task = split_to_intervals(workflow, ref, interval_list, intervals, window, parents)
    intervals = []
    with open(intervals_file,'r') as f:
        for line in f:
            attrs  = line.split()
            token = attrs[0]+':'+attrs[1]+'-'+attrs[2]
            intervals.append(token)


    snp_calling_tasks = []
    snp_calling_task = snp_calling(workflow, ref, bam, split_dir, merged_dir, intervals, parents=parents)
    #snp_calling_tasks.append(snp_calling_task)
       
    if clean:
        clean_split_vcfs(workflow, snp_calling_task, split_dir)
 
    return


'''
----------------------------------------
            *** MAIN ***
----------------------------------------
'''

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-B', '--bam', required=True)
    parser.add_argument('-O', '--out_dir', required=True)
    parser.add_argument('-R', '--ref', required=True)
    parser.add_argument('-L', '--intervals', required=True, default=None, \
            help="Regions to call intervals on. Helpful for parallelization as \
            Each interval is processed independently")
    parser.add_argument('-dbsnp','--dbsnp',  required=True)
    parser.add_argument('-clean','--clean',  required=False, action='store_true', help='Removes split vcf files')
    parser.add_argument('-t','--max_cores', required=False, default=16, type=int)
    parser.add_argument('-drm', required=False, default="local")
    args = parser.parse_args()

    global dbsnp  
    dbsnp=args.dbsnp
    
    # create output directory
    os.system('mkdir -p %s' % args.out_dir)
    
    # init cosmos
    cosmos = Cosmos('sqlite:///%s/sqlite.db' % args.out_dir,
                    default_drm=args.drm, default_max_attempts=2)
    cosmos.initdb() 

    # workflow
    workflow = cosmos.start('snp_calls', restart=False, skip_confirm=True)
    snps_n_indels_workflow(workflow, args.ref, args.bam, args.out_dir, args.intervals,  args.clean) 
    workflow.run(set_successful=True, dry=False, max_cores=args.max_cores)
    sys.exit(0 if workflow.successful else 1) 
