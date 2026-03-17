### A example pipeline to identify DGRs in 29 known DGRs from human gut virome
    ./MetaCSSTmain -build config.json -in example/hv29.fa -out example/out_tmp1 -thread 4
    python3 src/call_vr.py example/out_tmp1/raw.gtf example/hv29.fa example/out-DGR.gtf
    
