#########################################
# Change below directories accordingly.
#########################################
CONFIG={}

# Make sure those directories are exist.
CONFIG["result-dir"] = "./result/" # all the running result will be stored in this directory.
CONFIG["report-dir"] = "./report/" # all the reports will be stored in this directory.
CONFIG["log-dir"] = "./log/"

CONFIG["ref-dir"] = "/srv/public/kjk/noname/data/refs/" # directory for genome files in .fa format
CONFIG["read-dir"] = "/srv/public/kjk/noname/paper/simulated/" # directory for read files
CONFIG["p-dir"] = "/srv/public/kjk/noname/data/refs/tp/" # directory for true positive (vcf)

CONFIG["read-mappers"] = {}
CONFIG["read-mappers"]["bwa"] = "/srv/public/kjk/noname/bin/bwa-mem/bwa" # bwa
CONFIG["read-mappers"]["stellar"] = "/srv/public/kjk/noname/bin/stellar/stellar" # stellar

CONFIG["variant-callers"] = {}
#CONFIG["variant-callers"]["delly"] = "/srv/public/kjk/noname/bin/delly/delly" # delly binary
#CONFIG["variant-callers"]["pindel"] = "/srv/public/kjk/noname/bin/pindel/pindel" # pindel binary 
#CONFIG["variant-callers"]["lumpy"] = "/srv/public/kjk/noname/bin/lumpy-sv/bin/lumpyexpress" # pindel binary 
#CONFIG["variant-callers"]["gasvpro"] = "/srv/public/kjk/noname/bin/gasv/bin/GASVPro.sh" # pindel binary 
#CONFIG["variant-callers"]["crest"] = "/srv/public/kjk/noname/bin/crest/CRESTRun.sh" # crest
#CONFIG["variant-callers"]["gustaf"] = "/srv/public/kjk/noname/bin/gustaf/gustaf" # gustaf binary
CONFIG["variant-callers"]["v_seqan"] = "/srv/public/kjk/vaquita-build-develop/bin/vaquita" # v binary

CONFIG["tools"] = {}
CONFIG["tools"]["samtools"] = "/srv/public/kjk/noname/bin/samtools/samtools" # samtools binary
CONFIG["tools"]["bcftools"] = "/srv/public/kjk/noname/bin/bcftools/bcftools/bcftools" #bcftools binary 
CONFIG["tools"]["sam2stellargff"] = "/srv/public/kjk/noname/bin/sam2stellargff.pl" # do not change this

CONFIG["genome-simulators"] = ["mason","varsim", "dgv"]

#CONFIG["read-white-list"] = [".5.readsL", ".30.readsL"]
#CONFIG["read-white-list"] = [".5.readsL"]
#CONFIG["read-white-list"] = [".30.readsL"]
#CONFIG["read-white-list"] = [".30"]
#CONFIG["read-white-list"] = ["_30x_"]
#CONFIG["read-white-list"] = ["5x","50x"]
#CONFIG["read-white-list"] = ["10x", "50x", ]
#CONFIG["read-white-list"] = ["10x", "50x", "30x","5x"]
#CONFIG["read-white-list"] = ["10x", "50x"]
#CONFIG["read-white-list"] = ["50x","30x", "10x", "5x"]
#CONFIG["read-white-list"] = ["ERR194147"]
CONFIG["read-white-list"] = ["ERR1941"]
#CONFIG["read-white-list"] = ["5x", "10x", "30x", "50x"]
#CONFIG["read-white-list"] = ["5x", "10x","30x"]
#CONFIG["read-white-list"] = ["5x","10x","30x","50x"]
#CONFIG["read-white-list"] = ["30x", "5x"]
#CONFIG["read-white-list"] = ["5x_1.", "10x_1.", "30x_1.", "50x_1."]
#CONFIG["read-white-list"] = ["r1.", "30x_1."]
#CONFIG["read-white-list"] = ["i1.fq", "30x_1."]
#CONFIG["read-white-list"] = ["ERR19416"]
#CONFIG["read-white-list"] = ["145"]
#CONFIG["read-white-list"] = ["50x_1."]
#CONFIG["read-white-list"] = ["30x_1."]
#CONFIG["read-white-list"] = ["5x_1.","30x_1."]
#CONFIG["read-white-list"] = ["30x_1."]
#CONFIG["read-white-list"] = ["30x_1.","50x_1."]
#CONFIG["read-white-list"] = ["i1.fq", "30x_1."]
CONFIG["read-black-list"] = []
CONFIG["vc-white-list"] = [""]
#CONFIG["vc-white-list"] = ["-pd:1.5"] 
#CONFIG["vc-white-list"] = ["pi:8.0"] 
#CONFIG["vc-white-list"] = ["v:3,-m:50,-a:50,-ps:501,-pi:5.0"] 
#CONFIG["vc-white-list"] = ["-rw:10#-c:4#-v:3"] 
#CONFIG["vc-white-list"] = ["use-global-th"] 
#CONFIG["vc-white-list"] = ["v_seqan#-m:30#-a:10#-ps:500#-pi:5.0#-cs:20#-ce:0.1#-cg:_srv_public_kjk_noname_data_refs_hg19.fa#-rw:50#-sr:0.5#-pr:1.0#-cr:0.5#-rm:30#-rx:1.0#-rs:0.2#-v:1"] 
#CONFIG["vc-white-list"] = ["v_seqan#-m:30#-a:20#-ps:500#-pi:5.0#-cs:20#-ce:0.1#-cg:_srv_public_kjk_noname_data_refs_hg19.fa#-rw:500#-sr:0.5#-pr:1.0#-cr:0.5#-rm:30#-rx:1.0#-rs:0.2#-v:1"] 
CONFIG["vc-black-list"] = [] 
CONFIG["mapper-white-list"] = ["bwa"] 
CONFIG["mapper-black-list"] = [] 
