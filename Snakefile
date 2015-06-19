shell.prefix("set -euo pipefail;")

sample= "test"

### Executables
longorfs    = "./src/TransDecoder-2.0.1/TransDecoder.LongOrfs"
predict     = "./src/TransDecoder-2.0.1/TransDecoder.Predict"
hmmscan     = "hmmscan"
signalp     = "/usr/local/src/signalp-4.1/signalp"
tmhmm       = "/usr/local/src/tmhmm-2.0c/bin/tmhmm"
rnammer_transcriptome= "./src/Trinotate-2.0.2/util/rnammer_support/RnammerTranscriptome.pl"
rnammer     = "/usr/local/src/rnammer-1.2/rnammer"
trinotate   = "./src/Trinotate-2.0.2/Trinotate"
###


### URLs
# TODO
###

### Folders
# TODO
###



rule all:
    input:
#        "data/transdecoder/test.pep",
#        "data/trinotate/test_pep_uniref90.tsv.gz",
#        "data/trinotate/test_rna_uniref90.tsv.gz",
#        "data/trinotate/test_pep_sprot.tsv.gz",
#        "data/trinotate/test_rna_sprot.tsv.gz",
#        "data/trinotate/test_pep_pfam.tsv.gz"
#        "data/trinotate/test_pep_loaded.txt",
#        "data/trinotate/test_rna_loaded.txt",
        "data/trinotate/test.tsv"



rule clean:
    shell:
        """
        rm -rf data/trinotate
        rm -rf data/transdecoder
        """



################################################################################
# Transdecoder                                                                 #
################################################################################

rule transdecoder_longorfs:
    input:
        assembly= "data/assembly/{sample}.fasta"
    output:
        orfs= "{sample}.fasta.transdecoder_dir/longest_orfs.pep"
    params:
        folder= "data/transdecoder"
    threads:
        1
    log:
        "data/transdecoder/{sample}_longorfs.log"
    shell:
        """
        mkdir -p {params.folder}
        
        {longorfs} -t {input.assembly}  \
        > {log} 2>&1
        """



rule transdecoder_pep_sprot:
    input:
        orfs=   "{sample}.fasta.transdecoder_dir/longest_orfs.pep",
        db=     "data/db/uniprot_sprot.trinotate.pep"
    output:
        table=  "data/transdecoder/{sample}_sprot.tsv.gz"
    threads:
        24
    params:
    log:
        "data/transdecoder/{sample}_sprot.log"
    shell:
        """
        mkdir -p data/transdecoder
        
        blastp                      \
            -query {input.orfs}     \
            -db    {input.db}       \
            -max_target_seqs 1      \
            -outfmt 6               \
            -evalue 1e-5            \
            -num_threads {threads}  \
        2> {log}                    |
        gzip -9 > {output.table}
        """



rule transdecoder_pep_uniref90:
    input:
        orfs=   "{sample}.fasta.transdecoder_dir/longest_orfs.pep",
        db=     "data/db/uniprot_uniref90.trinotate.pep"
    output:
        table=  "data/transdecoder/{sample}_uniref90.tsv.gz"
    params:
        folder= "data/transdecoder"
    threads:
        24
    params:
    log:
        "data/transdecoder/{sample}_uniref90.log"
    shell:
        """
        mkdir -p {params.folder}
        
        blastp                                  \
            -query              {input.orfs}    \
            -db                 {input.db}      \
            -max_target_seqs    1               \
            -outfmt             6               \
            -evalue             1e-5            \
            -num_threads        {threads}       \
        2> {log}                                |
        gzip -9 > {output.table}
        """



rule transdecoder_pep_pfam:
    input:
        orfs=   "{sample}.fasta.transdecoder_dir/longest_orfs.pep",
        db=     "data/db/Pfam-A.hmm"
    output:
        table=  "data/transdecoder/{sample}_pep_pfam.tsv.gz"
    params:
        folder= "data/transdecoder"
    log:
        "data/transdecoder/{sample}_pfam.log"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {hmmscan}                                   \
            --cpu       {threads}                   \
            --domtblout >(gzip -9 > {output.table}) \
            {input.db}                              \
            {input.orfs}                            \
        > {log} 2>&1
        """



rule transdecoder_predict:
    input:
        assembly=       "data/assembly/{sample}.fasta",
        uniref_table=   "data/transdecoder/{sample}_uniref90.tsv.gz",
        pfam_table=     "data/transdecoder/{sample}_pep_pfam.tsv.gz"
    output:
        bed=    "data/transdecoder/{sample}.bed",
        cds=    "data/transdecoder/{sample}.cds",
        gff3=   "data/transdecoder/{sample}.gff3",
        mRNA=   "data/transdecoder/{sample}.mRNA",
        pep=    "data/transdecoder/{sample}.pep"
    params:
        folder= "data/transdecoder",
        bed=    "{sample}.fasta.transdecoder.bed",
        cds=    "{sample}.fasta.transdecoder.cds",
        gff3=   "{sample}.fasta.transdecoder.gff3",
        mRNA=   "{sample}.fasta.transdecoder.mRNA",
        pep=    "{sample}.fasta.transdecoder.pep"
    threads:
        1
    log:
        "data/transdecoder/{sample}_predict.log"
    shell:
        """
        {predict}                                           \
            -t                      {input.assembly}        \
            --retain_pfam_hits      {input.pfam_table}      \
            --retain_blastp_hits    {input.uniref_table}    \
        > {log} 2>&1
        
        mv  {params.bed}  {output.bed}
        mv  {params.cds}  {output.cds}
        mv  {params.gff3} {output.gff3}
        mv  {params.mRNA} {output.mRNA}
        mv  {params.pep}  {output.pep}
        
        rm -rf {sample}.fasta.transdecoder_dir
        """



################################################################################
# Trinotate                                                                    #
################################################################################

rule trinotate_pep_sprot:
    input:
        pep= "data/transdecoder/{sample}.pep",
        db=  "data/db/uniprot_sprot.trinotate.pep"
    output:
        table= "data/trinotate/{sample}_pep_sprot.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        "data/trinotate/{sample}_pep_sprot.log"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        blastp                              \
            -query              {input.pep} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2> {log}                            |
        gzip -9 > {output.table}
        """



rule trinotate_rna_sprot:
    input:
        rna= "data/assembly/{sample}.fasta",
        db=  "data/db/uniprot_sprot.trinotate.pep"
    output:
        table= "data/trinotate/{sample}_rna_sprot.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        "data/trinotate/{sample}_rna_sprot.log"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        blastx                              \
            -query              {input.rna} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2> {log}                            |
        gzip -9 > {output.table}
        """



rule trinotate_pep_uniref90:
    input:
        pep= "data/transdecoder/{sample}.pep",
        db=  "data/db/uniprot_uniref90.trinotate.pep"
    output:
        table= "data/trinotate/{sample}_pep_uniref90.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        "data/trinotate/{sample}_pep_uniref90.log"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        blastp                              \
            -query              {input.pep} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           |
        gzip -9 > {output.table}            \
        2> {log}
        """



rule trinotate_rna_uniref90:
    input:
        rna= "data/assembly/{sample}.fasta",
        db=  "data/db/uniprot_uniref90.trinotate.pep"
    output:
        table= "data/trinotate/{sample}_rna_uniref90.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        "data/trinotate/{sample}_rna_uniref90.log"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        blastx                              \
            -query              {input.rna} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2> {log}                            |
        gzip -9 > {output.table}
        """



rule trinotate_pep_pfam:
    input:
        pep= "data/transdecoder/{sample}.pep",
        db=  "data/db/Pfam-A.hmm"
    output:
        table= "data/trinotate/{sample}_pep_pfam.tsv.gz"
    params:
    log:
        "data/trinotate/{sample}_pep_pfam.log"
    threads:
        24
    shell:
        """
        {hmmscan}                                   \
            --cpu       {threads}                   \
            --domtblout >(gzip -9 > {output.table}) \
            {input.db}                              \
            {input.pep}                             \
        > {log}
        """



rule trinotate_pep_signalp:
    input:
        pep= "data/transdecoder/{sample}.pep",
    output:
        table= "data/trinotate/{sample}_pep_signalp.tsv"
    log:
        "data/trinotate/{sample}_pep_signalp.log"
    threads:
        1
    shell:
        """
        {signalp}               \
            -f  short           \
            -n  {output.table}  \
            {input.pep}         \
        > {log} 2>&1
        """



rule trinotate_pep_tmhmm:
    input:
        pep= "data/transdecoder/{sample}.pep"
    output:
        table= "data/trinotate/{sample}_pep_tmhmm.tsv"
    log:
        "data/trinotate/{sample}_pep_tmhmm.log"
    threads:
        1
    shell:
        """
        {tmhmm}             \
            --short         \
        < {input.pep}       \
        > {output.table}    \
        2> {log}
        """



rule trinotate_rna_rnammer:
    input:
        fasta= "data/assembly/{sample}.fasta"
    output:
        gff= "data/trinotate/{sample}_rna_rnammer.gff"
    threads:
        24 # Because it uses nonspecific temporary filenames
    params:
        tmp=   "tmp.superscaff.rnammer.gff",
        fasta= "transcriptSuperScaffold.fasta",
        bed=   "transcriptSuperScaffold.bed"
    log:
        "data/trinotate/{sample}_rna_rnammer.log"
    shell:
        """
        {rnammer_transcriptome}             \
            --transcriptome {input.fasta}   \
            --path_to_rnammer {rnammer}     \
        > {log} 2>&1
        
        mv {sample}.fasta.rnammer.gff {output.gff}
        
        rm {params}
        """



rule trinotate_unzip_db:
    input:
        db_gz=  "data/db/Trinotate.sqlite.gz"
    output:
        db=     "data/trinotate/{sample}.sqlite"
    params:
        folder= "data/trinotate"
    log:
                "data/trinotate/{sample}.sqlite.log"
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        
        gzip -dc {input.db_gz}  \
        > {output.db}           \
        2> {log}
        """



rule trinotate_init_database:
    input:
        db=         "data/trinotate/{sample}.sqlite",
        assembly=   "data/assembly/{sample}.fasta",
        pep=        "data/transdecoder/{sample}.pep",
        mapfile=    "data/assembly/{sample}_gene_to_trans_map.tsv"
    output:
        mock=       "data/trinotate/{sample}_db_init.txt"
    params:
        folder=     "data/trinotate"
    log:
                    "data/trinotate/{sample}_db_init.txt"
    threads:
        24 # Avoid other threads
    shell:
        """
        {trinotate}                                 \
            {input.db}                              \
            init                                    \
            --gene_trans_map    {input.mapfile}     \
            --transcript_fasta  {input.assembly}    \
            --transdecoder_pep  {input.pep}         \
        > {log} 2>&1
        
        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_results:
    input:
        mock=       "data/trinotate/{sample}_db_init.txt",
        db=         "data/trinotate/{sample}.sqlite",
        sprot=      "data/trinotate/{sample}_pep_sprot.tsv.gz",
        trembl=     "data/trinotate/{sample}_pep_uniref90.tsv.gz",
        pfam=       "data/trinotate/{sample}_pep_pfam.tsv.gz",
        tmhmm=      "data/trinotate/{sample}_pep_tmhmm.tsv",
        signalp=    "data/trinotate/{sample}_pep_signalp.tsv"
    output:
        mock=       "data/trinotate/{sample}_pep_loaded.txt"
    threads:
        24 # Avoid other accesses
    params:
        folder= "data/trinotate"
    log:
        "data/trinotate/{sample}_pep_loaded.log"
    shell:
        """
        mkdir -p {params.folder}
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_swissprot_blastp   <(gzip -dc {input.sprot})   \
        > {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_trembl_blastp      <(gzip -dc {input.trembl})  \
        >> {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_pfam               <(gzip -dc {input.pfam})    \
        >> {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_tmhmm              {input.tmhmm}               \
        >> {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_signalp            {input.signalp}             \
        >> {log} 2>&1
        
        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_transcript_results:
    input:
        mock=       "data/trinotate/{sample}_db_init.txt",
        db=         "data/trinotate/{sample}.sqlite",
        sprot=      "data/trinotate/{sample}_rna_sprot.tsv.gz",
        trembl=     "data/trinotate/{sample}_rna_uniref90.tsv.gz",
        rnammer=    "data/trinotate/{sample}_rna_rnammer.gff"
    output:
        mock=   "data/trinotate/{sample}_rna_loaded.txt"
    params:
        folder= "data/trinotate"
    threads:
        24 # Avoid other threads
    log:
        "data/trinotate/{sample}_rna_loaded.txt"
    shell:
        """
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_swissprot_blastp   <(gzip -dc {input.sprot})   \
        > {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_trembl_blastx      <(gzip -dc {input.trembl})  \
        >> {log} 2>&1
        
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_rnammer            {input.rnammer}             \
        >> {log} 2>&1
        
        printf "This is a mock file" > {output.mock}
        """



rule trinotate_generate_report:
    input:
        db=         "data/trinotate/{sample}.sqlite",
        mock_pep=   "data/trinotate/{sample}_pep_loaded.txt",
        mock_rna=   "data/trinotate/{sample}_rna_loaded.txt"
    output:
        report=     "data/trinotate/{sample}.tsv"
    params:
        evalue=         "1e-5",
        pfam_cutoff=    "DNC" # DNC DGC DTC SNC SGC STC
    log:
        "data/trinotate/{sample}_report.log"
    threads:
        24 # Avoid other threads
    shell:
        """
        {trinotate}                                 \
            {input.db}                              \
            report                                  \
            -E              {params.evalue}         \
            --pfam_cutoff   {params.pfam_cutoff}    \
        > {output.report}                           \
        2> {log}
        """
