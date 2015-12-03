import PipelineCaptureC as Capt_C

if __name__ == "__main__":
    
    Capt_C.sites2fragments('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/digest.dir/digest.bed.gz',
                           '/shared/sudlab1/General/annotations/hg19_ensembl75/contigs.tsv',
                           '/fastdata/mbp15ja/capturec-pilot/out_fragments.bed.gz')
#     Capt_C.countInteractions('/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/dedupped.dir/IGF0003828.bwa.by_name.bam', 
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/digest.dir/fragments.bed.gz',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/probe_fragments.bed.gz',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/test_output1.txt',
#                               '/fastdata/mbp15ja/capturec-pilot/capture_c_sudlab_pe_raw/test_output2.txt')
#     
