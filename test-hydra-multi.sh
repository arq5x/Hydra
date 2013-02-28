echo "Downloading 3 sample files from 1000 Genomes (~1.5Gb total)...\c"
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00419/alignment/HG00419.chrom11.ILLUMINA.bwa.CHS.low_coverage.20121211.bam
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG01615/alignment/HG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam
echo "done"


echo "creating a basic configuration file from the downloaded 100G files...\c"
echo "HG00096\tHG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
HG00419\tHG00419.chrom11.ILLUMINA.bwa.CHS.low_coverage.20121211.bam
HG01615\tHG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam" > config.stub.txt
echo "done"


echo "creating a complete configuration file by sampling BAM to create library stats...\c"
python scripts/make_hydra_config.py -i config.stub.txt > config.hydra.txt
echo "done"


echo "extracting discordant alignments from BAM files...\c"
python scripts/extract_discordants.py -i config.hydra.txt
echo "done"


echo "running hydra router on the discordant alignments...\c"
hydra-router -config config.hydra.txt -routedList routed-files.txt
echo "done"


echo "running hydra-assembler on the routed files of discordant alignments...\c"
sh scripts/assemble-routed-files.sh config.hydra.txt routed-files.txt
echo "done"


echo "re-combining the individual assembled files...\c"
sh scripts/combine-assembled-files.sh   ./   all.1000G.assembled
echo "done"


echo "finalizing SV breakpoint calls...\c"
python scripts/finalizeBreakpoints.py -i all.1000G.assembled -o all.1000G.sv
echo "done"

# final
#11      1914926 1915263 11      1936982 1937295 19      6       +       -       0       0       1       1       37      37      22370   6       6.0     6       6.0     6       0       0


# detail
#11      1914926 1914998 11      1937045 1937121 ERR016031.2210220       1       +       -       2       1       1       1       37      37      1       Y       19
#11      1915040 1915142 11      1937006 1937084 ERR023221.8596282       2       +       -       1       0       1       1       37      37      1       Y       19
#11      1915051 1915151 11      1936982 1937074 ERR023222.25347361      1       +       -       0       1       1       1       37      37      1       Y       19
#11      1915088 1915158 11      1937219 1937295 ERR016031.9514135       2       +       -       1       0       1       1       37      37      1       Y       19
#11      1915089 1915191 11      1936998 1937089 ERR023222.29574320      2       +       -       1       0       1       1       37      37      1       Y       19
#11      1915202 1915263 11      1937103 1937205 ERR023222.17735071      1       +       -       0       2       1       1       37      37      1       Y       19

# 1000G call
http://genome.ucsc.edu/cgi-bin/hgc?hgsid=327887259&c=chr11&o=1915273&t=1936956&g=tgpPhase1&i=G%2F%3CDEL%3E