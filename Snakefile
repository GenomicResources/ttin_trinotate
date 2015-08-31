shell.prefix("set -euo pipefail;")

SAMPLES, = glob_wildcards("./data/assembly/{sample}.fasta")


### Executables
longorfs    = "./src/TransDecoder-2.0.1/TransDecoder.LongOrfs"
blastp      = "blastp"
blastx      = "blastx"
predict     = "./src/TransDecoder-2.0.1/TransDecoder.Predict"
hmmscan     = "hmmscan"
signalp     = "./src/signalp-4.1/signalp"
tmhmm       = "./src/tmhmm-2.0c/bin/tmhmm"
rnammer_transcriptome = "./src/Trinotate-2.0.2/util/rnammer_support/RnammerTranscriptome.pl"
rnammer     = "./src/rnammer-1.2/rnammer"
trinotate   = "./src/Trinotate-2.0.2/Trinotate"
gzip        = "pigz"
###


### URLs
# TODO
###

### Folders
# TODO
###



rule all:
    input:
        expand("data/trinotate/{sample}.tsv", sample= SAMPLES)



rule only_transdecoder:
    input:
        expand("data/transdecoder/{sample}.pep", sample= SAMPLES)



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
        assembly=   "data/assembly/{sample}.fasta"
    output:
        orfs=       "{sample}.fasta.transdecoder_dir/longest_orfs.pep"
    params:
        folder=     "data/transdecoder"
    threads:
        1
    log:
        out=        "data/transdecoder/{sample}_longorfs.out",
        err=        "data/transdecoder/{sample}_longorfs.err"
    shell:
        """
        mkdir -p {params.folder}
        
        {longorfs} -t {input.assembly}  \
        >   {log.out}                   \
        2>  {log.err}
        """



#rule transdecoder_pep_sprot:
#	"""
#	Written just in case you do not want to use the Blastp against Uniref
#	"""
#    input:
#        orfs=   "{sample}.fasta.transdecoder_dir/longest_orfs.pep",
#        db=     "data/db/uniprot_sprot.trinotate.pep"
#    output:
#        table=  "data/transdecoder/{sample}_sprot.tsv.gz"
#    threads:
#        24
#    log:
#        "data/transdecoder/{sample}_sprot.log"
#    shell:
#        """
#        mkdir -p data/transdecoder
#         
#        {blastp}                    \
#            -query {input.orfs}     \
#            -db    {input.db}       \
#            -max_target_seqs 1      \
#            -outfmt 6               \
#            -evalue 1e-5            \
#            -num_threads {threads}  \
#        2> {log}                    |
#        {gzip} -9 > {output.table}
#        """



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
    log:
        out=    "data/transdecoder/{sample}_uniref90.out",
        err=    "data/transdecoder/{sample}_uniref90.err"
    shell:
        """
        mkdir -p {params.folder}
        
        {blastp}                                \
            -query              {input.orfs}    \
            -db                 {input.db}      \
            -max_target_seqs    1               \
            -outfmt             6               \
            -evalue             1e-5            \
            -num_threads        {threads}       \
        2>  {log.err}                           |
        {gzip} -9 > {output.table}

		touch {log.out}
        """



rule transdecoder_pep_pfam:
    input:
        orfs=   "{sample}.fasta.transdecoder_dir/longest_orfs.pep",
        db=     "data/db/Pfam-A.hmm"
    output:
        table=  "data/transdecoder/{sample}_pfam.tsv.gz"
    params:
        folder= "data/transdecoder"
    log:
        out=    "data/transdecoder/{sample}_pfam.out",
        err=    "data/transdecoder/{sample}_pfam.err"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {hmmscan}                                       \
            --cpu       {threads}                       \
            --domtblout >({gzip} -9 > {output.table})   \
            {input.db}                                  \
            {input.orfs}                                \
        >   {log.out}                                   \
        2>  {log.err}
        """



rule transdecoder_predict:
    input:
        assembly=       "data/assembly/{sample}.fasta",
        uniref_table=   "data/transdecoder/{sample}_uniref90.tsv.gz",
        pfam_table=     "data/transdecoder/{sample}_pfam.tsv.gz"
    output:
        bed=    "data/transdecoder/{sample}.bed",
        cds=    "data/transdecoder/{sample}.cds",
        gff3=   "data/transdecoder/{sample}.gff3",
        mRNA=   "data/transdecoder/{sample}.mRNA",
        pep=    "data/transdecoder/{sample}.pep"
    params:
        folder= "data/transdecoder",
        tmp_dir="{sample}.fasta.transdecoder_dir",
        bed=    "{sample}.fasta.transdecoder.bed",
        cds=    "{sample}.fasta.transdecoder.cds",
        gff3=   "{sample}.fasta.transdecoder.gff3",
        mRNA=   "{sample}.fasta.transdecoder.mRNA",
        pep=    "{sample}.fasta.transdecoder.pep"
    threads:
        1
    log:
        out=    "data/transdecoder/{sample}_predict.out",
        err=    "data/transdecoder/{sample}_predict.err"
    shell:
        """
        {predict}                                           		    \
            -t                      {input.assembly}        		    \
            --retain_pfam_hits      <({gzip} -dc {input.pfam_table})    \
            --retain_blastp_hits    <({gzip} -dc {input.uniref_table})	\
        >   {log.out}                                                   \
        2>  {log.err}
        
        mv  {params.bed}  {output.bed}
        mv  {params.cds}  {output.cds}
        mv  {params.gff3} {output.gff3}
        mv  {params.mRNA} {output.mRNA}
        mv  {params.pep}  {output.pep}
        
        rm -rf {params.tmp_dir}
        """






################################################################################
# Trinotate                                                                    #
################################################################################

rule trinotate_pep_sprot:
    input:
        pep=    "data/transdecoder/{sample}.pep",
        db=     "data/db/uniprot_sprot.trinotate.pep"
    output:
        table=  "data/trinotate/{sample}_pep_sprot.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        out=    "data/trinotate/{sample}_pep_sprot.out",
        err=    "data/trinotate/{sample}_pep_sprot.err"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {blastp}                            \
            -query              {input.pep} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2>  {log.err}                       |
        {gzip} -9 > {output.table}

		touch {log.out}
        """



rule trinotate_rna_sprot:
    input:
        rna=    "data/assembly/{sample}.fasta",
        db=     "data/db/uniprot_sprot.trinotate.pep"
    output:
        table=  "data/trinotate/{sample}_rna_sprot.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        out=    "data/trinotate/{sample}_rna_sprot.out",
        err=    "data/trinotate/{sample}_rna_sprot.err"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {blastx}                            \
            -query              {input.rna} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2>  {log.err}                       |
        {gzip} -9 > {output.table}

		touch {log.out}
        """



rule trinotate_pep_uniref90:
    input:
        pep=    "data/transdecoder/{sample}.pep",
        db=     "data/db/uniprot_uniref90.trinotate.pep"
    output:
        table=  "data/trinotate/{sample}_pep_uniref90.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        out=    "data/trinotate/{sample}_pep_uniref90.out",
        err=    "data/trinotate/{sample}_pep_uniref90.err"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {blastp}                            \
            -query              {input.pep} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2>  {log.err}                       |
        {gzip} -9 > {output.table}

		touch {log.out}
        """



rule trinotate_rna_uniref90:
    input:
        rna=    "data/assembly/{sample}.fasta",
        db=     "data/db/uniprot_uniref90.trinotate.pep"
    output:
        table=  "data/trinotate/{sample}_rna_uniref90.tsv.gz"
    params:
        folder= "data/trinotate"
    log:
        out=    "data/trinotate/{sample}_rna_uniref90.out",
        err=    "data/trinotate/{sample}_rna_uniref90.err"
    threads:
        24
    shell:
        """
        mkdir -p {params.folder}
        
        {blastx}                            \
            -query              {input.rna} \
            -db                 {input.db}  \
            -num_threads        {threads}   \
            -max_target_seqs    1           \
            -outfmt             6           \
        2>  {log.err}                       |
        {gzip} -9 > {output.table}

		touch   {log.out}
        """



rule trinotate_pep_pfam:
    input:
        pep=    "data/transdecoder/{sample}.pep",
        db=     "data/db/Pfam-A.hmm"
    output:
        table=  "data/trinotate/{sample}_pep_pfam.tsv.gz"
    params:
    log:
        out=    "data/trinotate/{sample}_pep_pfam.out",
        err=    "data/trinotate/{sample}_pep_pfam.err"
    threads:
        24
    shell:
        """
        {hmmscan}                                       \
            --cpu       {threads}                       \
            --domtblout >({gzip} -9 > {output.table})   \
            {input.db}                                  \
            {input.pep}                                 \
        >   {log.out}                                   \
        2>  {log.err}
        """



rule trinotate_pep_signalp:
    input:
        pep=    "data/transdecoder/{sample}.pep",
    output:
        table=  "data/trinotate/{sample}_pep_signalp.tsv"
    log:
        out=    "data/trinotate/{sample}_pep_signalp.out",
        err=    "data/trinotate/{sample}_pep_signalp.err"
    threads:
        1
    shell:
        """
        {signalp}               \
            -f  short           \
            -n  {output.table}  \
            {input.pep}         \
        >   {log.out}           \
        2>  {log.err}
        """



rule trinotate_pep_tmhmm:
    input:
        pep=    "data/transdecoder/{sample}.pep"
    output:
        table=  "data/trinotate/{sample}_pep_tmhmm.tsv"
    log:
        out=    "data/trinotate/{sample}_pep_tmhmm.out",
        err=    "data/trinotate/{sample}_pep_tmhmm.err"
    threads:
        1
    shell:
        """
        {tmhmm}             \
            --short         \
        < {input.pep}       \
        > {output.table}    \
        2> {log.err}
        
        rm -rf TMHMM_*
        """



rule trinotate_rna_rnammer:
    input:
        fasta=  "data/assembly/{sample}.fasta"
    output:
        gff=    "data/trinotate/{sample}_rna_rnammer.gff"
    threads:
        24 # Because it uses nonspecific temporary filenames
    params:
        gff=    "{sample}.fasta.rnammer.gff",
        tmp=    "tmp.superscaff.rnammer.gff",
        fasta=  "transcriptSuperScaffold.fasta",
        bed=    "transcriptSuperScaffold.bed"
    log:
        out=    "data/trinotate/{sample}_rna_rnammer.out",
        err=    "data/trinotate/{sample}_rna_rnammer.err"
    shell:
        """
        {rnammer_transcriptome}                 \
            --transcriptome     {input.fasta}   \
            --path_to_rnammer   {rnammer}       \
        >   {log.out}                           \
        2>  {log.err}
        
        mv {params.gff} {output.gff}
        
        rm {params.tmp} {params.fasta} {params.bed}
        """



rule trinotate_unzip_db:
    input:
        db_gz=  "data/db/Trinotate.sqlite.gz"
    output:
        db=     "data/trinotate/{sample}.sqlite"
    params:
        folder= "data/trinotate"
    log:
        out=    "data/trinotate/{sample}_sqlite.out",
        err=    "data/trinotate/{sample}_sqlite.err"
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        
        {gzip} -dc {input.db_gz}    \
        >   {output.db}             \
        2>  {log.err}
        """



rule trinotate_init_database:
    input:
        db=         "data/trinotate/{sample}.sqlite",
        assembly=   "data/assembly/{sample}.fasta",
        pep=        "data/transdecoder/{sample}.pep",
        mapfile=    "data/assembly/{sample}_gene_to_trans_map.tsv"
    output:
        mock=       "data/trinotate/{sample}_db_init.txt"
    log:
        out=        "data/trinotate/{sample}_db_init.out",
        err=        "data/trinotate/{sample}_db_init.err"
    threads:
        24 # Lock the database
    shell:
        """
        {trinotate}                                 \
            {input.db}                              \
            init                                    \
            --gene_trans_map    {input.mapfile}     \
            --transcript_fasta  {input.assembly}    \
            --transdecoder_pep  {input.pep}         \
        >   {log.out}                               \
        2>  {log.err}
        
        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_sprot:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        sprot=  "data/trinotate/{sample}_pep_sprot.tsv.gz"
    output:
        mock=   "data/trinotate/{sample}_pep_sprot_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_pep_sprot_loaded.out",
        err=    "data/trinotate/{sample}_pep_sprot_loaded.err"
    shell:
        """
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_swissprot_blastp   <({gzip} -dc {input.sprot}) \
        >   {log.out}                                           \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_uniref90:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        trembl= "data/trinotate/{sample}_pep_uniref90.tsv.gz"
    output:
        mock=   "data/trinotate/{sample}_pep_uniref90_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_pep_uniref90_loaded.out",
        err=    "data/trinotate/{sample}_pep_uniref90_loaded.err"
    shell:
        """
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_trembl_blastp   <({gzip} -dc {input.trembl})   \
        >   {log.out}                                           \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_pfam:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        pfam=   "data/trinotate/{sample}_pep_pfam.tsv.gz"
    output:
        mock=   "data/trinotate/{sample}_pep_pfam_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_pep_pfam_loaded.out",
        err=    "data/trinotate/{sample}_pep_pfam_loaded.err"
    shell:
        """
        {trinotate}                                 \
            {input.db}                              \
            LOAD_pfam   <({gzip} -dc {input.pfam})  \
        >   {log.out}                               \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_signalp:
    input:
        mock=       "data/trinotate/{sample}_db_init.txt",
        db=         "data/trinotate/{sample}.sqlite",
        signalp=    "data/trinotate/{sample}_pep_signalp.tsv"
    output:
        mock=       "data/trinotate/{sample}_pep_signalp_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=        "data/trinotate/{sample}_pep_signalp_loaded.out",
        err=        "data/trinotate/{sample}_pep_signalp_loaded.err"
    shell:
        """
        {trinotate}                         \
            {input.db}                      \
            LOAD_signalp    {input.signalp} \
        >   {log.out}                       \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_pep_tmhmm:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        tmhmm=  "data/trinotate/{sample}_pep_tmhmm.tsv"
    output:
        mock=   "data/trinotate/{sample}_pep_tmhmm_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_pep_tmhmm_loaded.out",
        err=    "data/trinotate/{sample}_pep_tmhmm_loaded.err"
    shell:
        """
        {trinotate}                     \
            {input.db}                  \
            LOAD_tmhmm    {input.tmhmm} \
        >   {log.out}                   \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_rna_sprot:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        sprot=  "data/trinotate/{sample}_rna_sprot.tsv.gz",
    output:
        mock=   "data/trinotate/{sample}_rna_sprot_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_rna_sprot_loaded.out",
        err=    "data/trinotate/{sample}_rna_sprot_loaded.err"
    shell:
        """
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_swissprot_blastx   <({gzip} -dc {input.sprot}) \
        >   {log.out}                                           \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_rna_uniref90:
    input:
        mock=   "data/trinotate/{sample}_db_init.txt",
        db=     "data/trinotate/{sample}.sqlite",
        trembl= "data/trinotate/{sample}_rna_uniref90.tsv.gz",
    output:
        mock=   "data/trinotate/{sample}_rna_uniref90_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=    "data/trinotate/{sample}_rna_uniref90_loaded.out",
        err=    "data/trinotate/{sample}_rna_uniref90_loaded.err"
    shell:
        """
        {trinotate}                                             \
            {input.db}                                          \
            LOAD_trembl_blastx  <({gzip} -dc {input.trembl})    \
        >   {log.out}                                           \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_load_rna_rnammer:
    input:
        mock=       "data/trinotate/{sample}_db_init.txt",
        db=         "data/trinotate/{sample}.sqlite",
        rnammer=    "data/trinotate/{sample}_rna_rnammer.gff",
    output:
        mock=       "data/trinotate/{sample}_rna_rnammer_loaded.txt"
    threads:
        24 # Lock the database
    log:
        out=        "data/trinotate/{sample}_rna_rnammer_loaded.out",
        err=        "data/trinotate/{sample}_rna_rnammer_loaded.err"
    shell:
        """
        {trinotate}                         \
            {input.db}                      \
            LOAD_rnammer {input.rnammer}    \
        >   {log.out}                       \
        2>  {log.err}

        printf "This is a mock file" > {output.mock}
        """



rule trinotate_generate_report:
    input:
        db=                 "data/trinotate/{sample}.sqlite",
        mock_pep_sprot=     "data/trinotate/{sample}_pep_sprot_loaded.txt",
        mock_pep_uniref90=  "data/trinotate/{sample}_pep_uniref90_loaded.txt",
        mock_pep_pfam=      "data/trinotate/{sample}_pep_pfam_loaded.txt",
        mock_pep_signalp=   "data/trinotate/{sample}_pep_signalp_loaded.txt",
        mock_pep_tmhmm=     "data/trinotate/{sample}_pep_tmhmm_loaded.txt",
        mock_rna_sprot=     "data/trinotate/{sample}_rna_sprot_loaded.txt",
        mock_rna_uniref90=  "data/trinotate/{sample}_rna_uniref90_loaded.txt",
        mock_rna_rnammer=   "data/trinotate/{sample}_rna_rnammer_loaded.txt",
    output:
        report=             "data/trinotate/{sample}.tsv"
    params:
        evalue=             "1e-5",
        pfam_cutoff=        "DNC" # DNC DGC DTC SNC SGC STC
    log:
        out=                "data/trinotate/{sample}_report.out",
        err=                "data/trinotate/{sample}_report.err"
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
        2> {log.err}
        """

