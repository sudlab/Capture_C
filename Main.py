import pydevd
import PipelineCaptureC as Capt_C

if __name__ == "__main__":
    print("Before settrace")
    pydevd.settrace("86.177.59.46")
    print("After settrace")
    Capt_C.countInteractions('/fastdata/mbp15ja/capturec-pilot/capture_c_mb1imsbams/dedupped.dir/Jurkat-MERGE-R1.by_name.bam', 
                              '/fastdata/mbp15ja/capturec-pilot/capture_c_mb1imsbams/digest.dir/fragments.bed.gz',
                              '/fastdata/mbp15ja/capturec-pilot/capture_c_mb1imsbams/probe_fragments.bed.gz',
                              '/fastdata/mbp15ja/capturec-pilot/capture_c_mb1imsbams/test_output1.txt',
                              '/fastdata/mbp15ja/capturec-pilot/capture_c_mb1imsbams/test_output2.txt')
    
    pydevd.settrace("86.177.59.46")
