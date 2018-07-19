#NEED TO CHECK SEQUENCE OUTPUT FILES AGAINST ACTUAL RESULTS

# install ncbi tools
# pull nr db if local blasting
# install clustal-omega (install argstable for this)

# $ export PATH=$PATH:/Users/jason/code/find-interruption-elements-side/NCBI_Resources/ncbi-blast-2.4.0+/bin/ ----- how to get this permanent?

# first requires the db to go through makeblastdb - can i do that here in python? (make it happen only if it isn't already made)
# makeblastdb -in genome -dbtype nucl -parse_seqids
# makeblastdb -in known_interrupted_genes_protein.faa -dbtype prot -parse_seqids
# and for nr
# add an error message so it prints out a more human-friendly message and the command required to create the blastdb

# make input either 1 genome, 1 nickname OR 1 table with genomes & nicknames

# could include blast options - xis_candidate_e_value, etc

from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast import NCBIXML
from difflib import SequenceMatcher
import argparse

EPILOG = '''
    Write something here maybe
    '''

# CHANGE - need to align read headers in full_xis_set.faa w/ these keys
recombinase_class = {
    '131C_1156': 'tyrosine',
    '131C_1121': 'tyrosine',
    '310F_2885': 'serine',
    '310F_2836': 'tyrosine',
    '310F_1894': 'tyrosine',
    '310F_2868': 'tyrosine',
    'Anacy_2120': 'serine',
    'Anacy_2144': 'serine',
    'Anacy_6026': 'serine',
    'Anacy_2214': 'tyrosine',
    'Anacy_2118': 'tyrosine',
    'Anacy_2116': 'tyrosine',
    'Anacy_1762': 'tyrosine',
    'ANA_C13501': 'serine',
    'ANA_C13532': 'tyrosine',
    'ANA_C13510': 'tyrosine',
    'ANA_C11423': 'tyrosine',
    'Ana7108_2785': 'tyrosine',
    'Ana7108_2847': 'tyrosine',
    'Ana7108_2850': 'tyrosine',
    'Ana7108_2755': 'tyrosine',
    'Ava_3928': 'tyrosine',
    'Cal7102DRAFT_03510': 'serine',
    'Cal7102DRAFT_03503': 'serine',
    'Cal7102DRAFT_00467': 'serine',
    'Cal7102DRAFT_07832': 'serine',
    'Cal7102DRAFT_03532': 'serine',
    'Cal7102DRAFT_02973': 'serine',
    'Cal7102DRAFT_00422': 'tyrosine',
    'Cal7102DRAFT_03508': 'tyrosine',
    'CalSC01_3760': 'serine',
    'CalSC01_4808': 'serine',
    'CalSC01_1389': 'serine',
    'CalSC01_2817': 'tyrosine',
    'Cal6303_5631': 'serine',
    'Cal6303_1564': 'serine',
    'Cal6303_2021': 'serine',
    'Cal6303_1693': 'serine',
    'Cal6303_0394': 'serine',
    'Cal6303_0220': 'serine',
    'Cal6303_0224': 'tyrosine',
    'Cal6303_0402': 'tyrosine',
    'Cal7103DRAFT_00003030': 'serine',
    'Cal7103DRAFT_00042020': 'tyrosine',
    'Cal7103DRAFT_00042040': 'tyrosine',
    'Cal7103DRAFT_00044010': 'serine',
    'Cal7103DRAFT_00042750': 'serine',
    'Cal7103DRAFT_00041060': 'serine',
    'Cal7103DRAFT_00055580': 'serine',
    'Cal7103DRAFT_00041480': 'serine',
    'Cal7103DRAFT_00006510': 'serine',
    'Cal7103DRAFT_00041270': 'serine',
    'Cal7103DRAFT_00041710': 'tyrosine',
    'Cal7507_0111': 'tyrosine',
    'UYCDRAFT_06247': 'tyrosine',
    'UYEDRAFT_00081': 'tyrosine',
    'CylstDRAFT_5431': 'tyrosine',
    'CylstDRAFT_5410': 'tyrosine',
    'CylstDRAFT_5455': 'tyrosine',
    'PCC9339DRAFT_03891': 'serine',
    'PCC9339DRAFT_03888': 'serine',
    'Fis9431DRAFT_4159': 'tyrosine',
    'Mas10914DRAFT_0641': 'tyrosine',
    'Mic7126DRAFT_0254': 'serine',
    'Mic7126DRAFT_0221': 'tyrosine',
    'N9414_15040': 'tyrosine',
    'N9414_15065': 'tyrosine',
    'N9414_14930': 'tyrosine',
    'Npun_R2637': 'serine',
    'Npun_F0392': 'tyrosine',
    'Nos7107_3361': 'serine',
    'Nos7107_3373': 'tyrosine',
    'alr1459': 'serine',
    'alr1442': 'tyrosine',
    'alr0677': 'tyrosine',
    'Nos7524_1241': 'serine',
    'Nos7524_1212': 'tyrosine',
    'Nos7524_1209': 'tyrosine',
    'Nos7524_1231': 'tyrosine',
    'RINTHH_3080': 'serine',
    'RintRC_2216': 'tyrosine',
    'Riv7116_6196': 'serine',
    'Riv7116_6219': 'serine',
    'Riv7116_6354': 'serine',
    'Riv7116_6308': 'serine',
    'Riv7116_6281': 'tyrosine',
    'Riv7116_6305': 'tyrosine',
    'Riv7116_6198': 'tyrosine',
    'Riv7116_6362': 'tyrosine',
    'Riv7116_6312': 'tyrosine',
    'WA1DRAFT_10544': 'serine',
    'WA1DRAFT_02505': 'serine',
    'WA1DRAFT_03849': 'serine',
    'WA1DRAFT_03874': 'serine',
    'WA1DRAFT_03870': 'tyrosine',
    'Tol9009DRAFT_00005880': 'serine',
    'Tol9009DRAFT_00061200': 'serine',
    'Tol9009DRAFT_00060970': 'serine',
    'Tol9009DRAFT_00060810': 'tyrosine',
    'Tol9009DRAFT_00004090': 'serine'
}

reference_genes = {
    '646565097': 'transposase (Ava_B0242)',
    '648049578': 'nifE',
    '648049580': 'nifK',
    '648049581': 'nifD',
    '648049582': 'nifH',
    '648049585': 'fdxN',
    '648049586': 'NifB',
    '648050508': 'coxA',
    '648050542': 'integrase family protein (Aazo_2682)',
    '648051374': 'hupL',
    '648051375': 'hupS',
    '648051418': 'hglE',
    '648051575': 'flv3B',
    '648052264': 'NADPH-dependent FMN reductase (Aazo_5221)',
    '2505798824': 'phospholipase D/transphosphatidylase (Cal7507_0570)',
    '2505803755': 'nifJ',
    '2505803979': 'FAD-dependent oxidoreductase (Cal7507_5656)',
    '2506492697': 'primase P4 (Ana7108_2845)',
    '2507085412': 'ATP-dependent DNA helicase (Pse6802_0098)',
    '2507088806': 'predicted integral membrane protein (Pse6802_3453)',
    '2507477837': 'caspase domain-containing protein (Cal7103DRAFT_00047390)',
    '2509768028': 'hypothetical protein (CylstDRAFT_1988)',
    '2509784323': 'arabinose efflux permease (Mic7126DRAFT_5075)',
    '2517243815': 'transposase (ISSoc8) (Mas10914DRAFT_5058)'
}


def getArgs():
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    parser.add_argument('--genome',
                        help="The fasta file of the genome to search for interruption elements")
    parser.add_argument('--name',
                        help="Short nickname for genome, to be used in output files")
    parser.add_argument('--xis_set',
                        help="The protein sequences of known xis genes")
    parser.add_argument('--protein_set',
                        help="The protein sequence database used to identify interrupted genes")
    parser.add_argument('--flank',
                        help="The distance away from candidate xis genes to search for interrupted genes. Deafault is 3000",
                        default=3000)
    args = parser.parse_args()
    return args


def find_xis_candidates(name):
    '''
    Given a genome and protein sequences of xis, returns xis candidates and their flanking regions from that genome.

    Does a BLAST of known xis genes against the genome.
    For each BLAST hit, walks the BLAST in each direction on the genome to identify the ORF containing the hit based on start/stop codons.
    Outputs the xis candidate ORFs to *xis_candidates.fna.
    Outputs the xis candidate ORFs with flanking regions to *xis_flank.fna.
    Passes the xis candidate ORFs with flanking regions to find_interrupted_gene() to look for interrupted genes in the flanking regions.
    '''

    # BLAST xis protein gene sequences into genome
    tblastn_cline = NcbitblastnCommandline(query=args.xis_set, db=args.genome, evalue=1e-20, outfmt=5, out=name + 'xis_tn_genome.xml')
    print(str(tblastn_cline))
    stdout, stderr = tblastn_cline()
    print(name + ':xis_set-genome BLAST complete')

    # from BLAST results, pull hits
    blast_results = open(name + 'xis_tn_genome.xml', 'r')
    blast_records = NCBIXML.parse(blast_results)
    xis_hits_dict = {}
    count = 0
    off_limits_ranges = {}
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            if alignment.title not in off_limits_ranges.keys():
                off_limits_ranges[alignment.title] = []
            for hsp in alignment.hsps:
                middle_of_hit = round((hsp.query_start+hsp.query_end)/2, 0)
                if middle_of_hit not in off_limits_ranges[alignment.title]:
                    count += 1
                    xis_hits_dict['hit' + str(count)] = {}
                    xis_hits_dict['hit' + str(count)]['title'] = alignment.title
                    xis_hits_dict['hit' + str(count)]['accession'] = alignment.accession
                    xis_hits_dict['hit' + str(count)]['start'] = hsp.sbjct_start
                    xis_hits_dict['hit' + str(count)]['end'] = hsp.sbjct_end
                    xis_hits_dict['hit' + str(count)]['top_xis_hit'] = blast_record.query
                    for i in range(hsp.query_start, hsp.query_end):
                        off_limits_ranges[alignment.title].append(i)
    blast_results.close()
    print(name + ':' + str(count) + ' xis candidates identified')

    # if no hits were found, then quite the program
    if count == 0:
        return

    # cut the xis_candidates from the contigs based on BLAST results
    xis_count = 0
    for hit in xis_hits_dict.values():
        xis_count += 1
        xis_flank = open(name + '-' + str(xis_count) + '_xis_flank.fna', 'w')
        xis_candidate_sequence = open(name + '-' + str(xis_count) + '_xis_candidate.fna', 'w')
        contig_title = hit['title']
        genome = open(args.genome, 'r')
        contig_sequence = sequence_cutter(genome, contig_title, 'all', 'FALSE')
        xis_coordinates = orf_finder(contig_sequence, hit['start'], hit['end'])
        genome.close()
        genome = open(args.genome, 'r')
        xis_fasta = sequence_cutter(genome, contig_title, xis_coordinates)
        xis_candidate_sequence.write(xis_fasta + '\n')
        genome.close()

        # log characteristics of the xis candidate to be used in later processing
        xis_dict = {}
        xis_dict['xis coordinates'] = xis_coordinates
        xis_dict['count'] = xis_count
        xis_dict['contig'] = hit['title']
        xis_dict['contig accession'] = hit['accession']
        xis_dict['class'] = recombinase_class[hit['top_xis_hit']]
        
        # cut the xis_candidates and flanking regions from the contigs based on BLAST results
        contig_sequence_length = len(contig_sequence)
        xis_start = xis_coordinates[0]
        xis_stop = xis_coordinates[1]

        if args.flank < xis_start:
            flank_start = xis_start - args.flank
        else:
            flank_start = 0

        if args.flank < contig_sequence_length:
            flank_end = xis_stop + args.flank
        else:
            flank_end = contig_sequence_length

        genome = open(args.genome, 'r')
        flank_coordinates = (flank_start, flank_end)
        flank_fasta = sequence_cutter(genome, contig_title, flank_coordinates)
        xis_flank.write(flank_fasta + '\n')
        genome.close()
        xis_flank.close()
        xis_candidate_sequence.close()

        # pass on the xis_candidates and flanking regions to search for interrupted genes in the flanking regions
        find_interrupted_gene(name, xis_dict, name + '-' + str(xis_count) + '_xis_flank.fna')


def sequence_cutter(genome, sequence_id, cut_coordinates='all', header='TRUE'):
    '''
    Given a fasta file, sequence identifier of a sequence in that file, and coordinate range, returns the given sequence cut to those coordinates.

    header option allows the user to note if they would like only the sequence without a header. Default is to return the sequence with a header that specifies the cut coordinates.
    '''

    count = 0
    for line in genome:
        if line.startswith('>' + sequence_id):
            sequence = []
            count += 1
            header_line = line.strip('\n')
            sequence_flag = 'TRUE'
            while sequence_flag == 'TRUE':
                try:
                    line = next(genome)
                    if line.startswith('>'):
                        sequence_flag = 'FALSE'
                    else:
                        for base in line.strip('\n'):
                            sequence.append(base)
                except Exception:
                    sequence_flag = 'FALSE'
    if count == 0:
        print('sequence name:' + sequence_id + ' not found in genome:' + str(genome))
    elif count == 1:
        if header == 'TRUE':
            if cut_coordinates == 'all':
                return(header_line + '\n' + ''.join(sequence[:]))
            else:
                start, stop = cut_coordinates
                return(header_line + ':' + str(cut_coordinates[0]) + ':' + str(cut_coordinates[1]) + '\n' + ''.join(sequence[start:stop + 1]))
        elif header == 'FALSE':
            if cut_coordinates == 'all':
                return(''.join(sequence[:]))
            else:
                start, stop = cut_coordinates
                return(''.join(sequence[start:stop + 1]))
    else:
        print('sequence name not unique')

def orf_finder(contig_sequence, blast_hit_start, blast_hit_end):
    '''
    Given a sequence and the start and stop of a BLAST hit, returns the ORF that the sequence range is a part of.

    Takes the start coordinate of the BLAST hit and traces the sequence back each codon at a time until it finds a start codon.
    Take the end coordinate of the BLAST hit and traces the sequence forward each codon at a time until it finds a stop codon.
    Returns the predicted ORF sequence.
    '''

    # do i need to switch things up if the frame is negative?
    start_codon = 'ATG'
    start_flag = 'FALSE'
    start_position = blast_hit_start
    while start_flag == 'FALSE':
        codon_to_check = contig_sequence[start_position:start_position+3]
        if codon_to_check == start_codon:
            start_flag = 'TRUE'
            orf_start_position = start_position
        start_position -= 3

    stop_codons = ['TAA', 'TAG', 'TGA']
    stop_flag = 'FALSE'
    stop_position = blast_hit_end
    while stop_flag == 'FALSE':
        codon_to_check = contig_sequence[stop_position:stop_position+3]
        if codon_to_check in stop_codons:
            stop_flag = 'TRUE'
            orf_stop_position = stop_position + 2
        stop_position += 3

    return orf_start_position, orf_stop_position


def find_interrupted_gene(name, xis_dict, xis_plus_flank):
    xis_count = xis_dict['count']

    # BLAST the xis + flanking region against the known reference genes for previously found interrupted genes
    blastx_cline = NcbiblastxCommandline(query=xis_plus_flank, db='interrupted_genes/known_interrupted_genes_protein.faa', culling_limit=15, evalue=1e-15, outfmt=5, out=name + '-' + str(xis_count) + '_xisflank_x_KnownIntGenes.xml')
    print(str(blastx_cline))
    stdout, stderr = blastx_cline()
    print(name + '-' + str(xis_count) + ':xis w/ flank - known Interrupted Genes BLAST complete')

    hitsdW = {}
    sortme = {}
    off_limits_range = []

    blast_results = open(name + '-' + str(xis_count) + '_xisflank_x_KnownIntGenes.xml', 'r')
    blast_records = NCBIXML.parse(blast_results)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                seq_flank_length = blast_record.query_length
                distance_from_edge = min(hsp.query_start, seq_flank_length - hsp.query_end)
                middle_of_hit = round((hsp.query_start+hsp.query_end)/2, 0)
                if middle_of_hit not in off_limits_range and hsp.query_start not in off_limits_range and hsp.query_end not in off_limits_range and distance_from_edge > 15:
                    sortme[alignment.accession] = round(hsp.align_length/alignment.length*100, 2)
                    hitsdW[alignment.accession] = {}
                    hitsdW[alignment.accession]['title'] = alignment.hit_def
                    hitsdW[alignment.accession]['start'] = hsp.query_start
                    hitsdW[alignment.accession]['end'] = hsp.query_end
                    hitsdW[alignment.accession]['cover'] = round(hsp.align_length/alignment.length*100, 2)
                    for i in range(hsp.query_start, hsp.query_end):
                        off_limits_range.append(i)
    blast_results.close()
  
    # if interrupted gene isn't a previously known gene, BLAST against a larger protein database
    if sortme == {}: # and some other criteria about the BLAST results?
        # for web blast, use...
        if args.protein_set == 'web':
            blastx_cline = NcbiblastxCommandline(remote='remote', query='xis_flank.fasta', db='refseq_protein', culling_limit=15, evalue=0.00000000000001, outfmt=5, out='xisflank_x_refseqprot.xml')
        # for local blast against a protein db, use...
        else:
            blastx_cline = NcbiblastxCommandline(query=xis_plus_flank, db=args.protein_set, culling_limit=15, evalue=1e-15, outfmt=5, out=name + '-' + str(xis_count) + '_xisflank_x_nr.xml')
        print(str(blastx_cline))
        stdout, stderr = blastx_cline()
        print(name + '-' + str(xis_count) + ':xis w/ flank - nr BLAST complete')
        # pull the best BLAST hit per region
        blast_results = open(name + '-' + str(xis_count) + '_xisflank_x_nr.xml', 'r')
        blast_records = NCBIXML.parse(blast_results)
        

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    seq_flank_length = blast_record.query_length
                    distance_from_edge = min(hsp.query_start, seq_flank_length - hsp.query_end)
                    middle_of_hit = round((hsp.query_start+hsp.query_end)/2, 0)
                    if middle_of_hit not in off_limits_range and hsp.query_start not in off_limits_range and hsp.query_end not in off_limits_range and distance_from_edge > 15:
                        sortme[alignment.accession] = round(hsp.align_length/alignment.length*100, 2)
                        hitsdW[alignment.accession] = {}
                        hitsdW[alignment.accession]['title'] = alignment.hit_def
                        hitsdW[alignment.accession]['start'] = hsp.query_start
                        hitsdW[alignment.accession]['end'] = hsp.query_end
                        hitsdW[alignment.accession]['cover'] = round(hsp.align_length/alignment.length*100, 2)
                        for i in range(hsp.query_start, hsp.query_end):
                            off_limits_range.append(i)
        blast_results.close()

        if sortme == {}:
            print(name + '-' + str(xis_count) + ':No reference genes found in flanking region')
            return

        # pull the hit from above with the lowest % coverage
        lowest_coverage = sorted(sortme.values())[0]
        for key in sortme.keys():
            if sortme[key] == lowest_coverage:
                lowest_coverage_accession = key

        # use the accession of the hit to pull the full gene sequence
        # CHANGE - need protein sequence to id contig regions, but nucleotide sequence to align
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=lowest_coverage_accession)
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()

        reference_gene_aa_file = open(name + '-' + str(xis_count) + '_reference_gene.faa', 'w')
        reference_gene_aa_file.write(">" + str(seq_record.id) + '\n' + str(seq_record.seq))
        reference_gene_aa_file.close()

        reference_gene_dna_file = open(name + '-' + str(xis_count) + '_reference_gene.fna', 'w')
        reference_gene_dna_file.write(">" + str(seq_record.id) + '\n' + str(seq_record.seq))
        reference_gene_dna_file.close()

    # if a previously known interrutped gene was found in the flanking regions, pull the DNA sequence of that gene
    else:
        lowest_coverage = sorted(sortme.values())[0]
        for key in sortme.keys():
            if sortme[key] == lowest_coverage:
                lowest_coverage_accession = key

        reference_gene_aa_file = open(name + '-' + str(xis_count) + '_reference_gene.faa', 'w')
        reference_gene_aa_file.write(sequence_cutter(open('interrupted_genes/known_interrupted_genes_protein.faa', 'r'), lowest_coverage_accession, 'all'))
        reference_gene_aa_file.close()

        reference_gene_dna_file = open(name + '-' + str(xis_count) + '_reference_gene.fna', 'w')
        reference_gene_dna_file.write(sequence_cutter(open('interrupted_genes/known_interrupted_genes_protein.fna', 'r'), lowest_coverage_accession, 'all'))
        reference_gene_dna_file.close()

    if lowest_coverage_accession in reference_genes.keys():
        xis_dict['reference gene'] = reference_genes[lowest_coverage_accession]
        print(name + '-' + str(xis_count) + ':reference gene identified: ' + reference_genes[lowest_coverage_accession])
    else:
        xis_dict['reference gene'] = lowest_coverage_accession
        print(name + '-' + str(xis_count) + ':reference gene identified: ' + lowest_coverage_accession)
    
    # pass on the reference gene to search the whole genome for all gene regions
    find_all_gene_regions(name, xis_dict, name + '-' + str(xis_count) + '_reference_gene.faa', name + '-' + str(xis_count) + '_reference_gene.fna')


def find_all_gene_regions(name, xis_dict, reference_gene_aa, reference_gene_dna):
    xis_count = xis_dict['count']
    xis_coordinates = xis_dict['xis coordinates']
    xis_contig = xis_dict['contig']

    # BLAST reference genes into genome
    tblastn_cline = NcbitblastnCommandline(query=reference_gene_aa, db=args.genome, evalue=1e-10, outfmt=5, out=name + '-' + str(xis_count) + "_referencegene_tn_genome.xml")
    print(str(tblastn_cline))
    stdout, stderr = tblastn_cline()
    print(name + '-' + str(xis_count) + ':reference gene - genome BLAST complete')

    # pull gene regions from blast hits
    hitsd = {}
    count = 0
    blast_results = open(name + '-' + str(xis_count) + '_referencegene_tn_genome.xml', 'r')
    blast_records = NCBIXML.read(blast_results)
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            count += 1
            hitsd['hit' + str(count)] = {}
            hitsd['hit' + str(count)]['contig_title'] = alignment.title
            hitsd['hit' + str(count)]['contig_start'] = hsp.sbjct_start
            hitsd['hit' + str(count)]['contig_end'] = hsp.sbjct_end
            if hsp.sbjct_start < xis_coordinates[0]:
                hitsd['hit' + str(count)]['gene_location'] = 'before'
            else:
                hitsd['hit' + str(count)]['gene_location'] = 'after'

    blast_results.close()

    # if only 1 gene section found in the genome, quit the program
    if count == 1:
        print(name + '-' + str(xis_count) + ':only 1 gene section located in genome')
        return

    # determine the gene region that is closest to the xis candidate gene
    dist_from_xis_list = []
    quick_dict = {}
    for hit in hitsd.keys():
        if hitsd[hit]['contig_title'] == xis_contig:
            quick_dict[hit] = {}
            dist_from_xis = min(abs(hitsd[hit]['contig_start'] - xis_coordinates[0]),
                    abs(hitsd[hit]['contig_start'] - xis_coordinates[1]),
                    abs(hitsd[hit]['contig_end'] - xis_coordinates[0]),
                    abs(hitsd[hit]['contig_end'] - xis_coordinates[1])
                    )
            dist_from_xis_list.append(dist_from_xis)
            quick_dict[hit]['dist_from_xis'] = dist_from_xis
    for hit in quick_dict.keys():
        if quick_dict[hit]['dist_from_xis'] == min(dist_from_xis_list):
            nearest_gene_section = hit

    # determine the gene region that is on the other side of the xis candidate gene (relative to the above gene region) 
    dist_from_xis_list = []
    quick_dict = {}
    for hit in hitsd.keys():
        if hitsd[hit]['contig_title'] == xis_contig and hitsd[hit]['gene_location'] != hitsd[nearest_gene_section]['gene_location']:
            quick_dict[hit] = {}
            dist_from_xis = min(abs(hitsd[hit]['contig_start'] - xis_coordinates[0]),
                    abs(hitsd[hit]['contig_start'] - xis_coordinates[1]),
                    abs(hitsd[hit]['contig_end'] - xis_coordinates[0]),
                    abs(hitsd[hit]['contig_end'] - xis_coordinates[1])
                    )
            dist_from_xis_list.append(dist_from_xis)
            quick_dict[hit]['dist_from_xis'] = dist_from_xis

    # if a second gene regions can't be found on the same contig as the xis candidate gene and the first gene region, take the best BLAST hti that isn't the first gene region
    if quick_dict == {}:
        if nearest_gene_section == 'hit1':
            second_gene_section = 'hit2'
        else:
            second_gene_section = 'hit1'
    else:
        for hit in quick_dict.keys():
            if quick_dict[hit]['dist_from_xis'] == min(dist_from_xis_list):
                second_gene_section = hit

    # CHANGE - do we want to adjust this in some way so that we capture the completed gene? do we need to rope ORF finder back in???
    # write the two gene regions to a file for alignment
    sequences_to_align = open(name + '-' + str(xis_count) + '_gene_sections_w_reference.fna', 'w')
    for hit in [nearest_gene_section, second_gene_section]:
        genome = open(args.genome, 'r')
        # extend the regions 20 bp beyond the blast hit to ensure overlap
        start_coord = hitsd[hit]['contig_start'] - 20
        end_coord = hitsd[hit]['contig_end'] + 20
        cut_coordinates = (start_coord, end_coord)
        contig_title = hitsd[hit]['contig_title']
        sequences_to_align.write(sequence_cutter(genome, contig_title, cut_coordinates) + '\n')
        genome.close()

    # write the reference gene for the interrupted gene to the file for alignment with teh gene regions
    for line in open(reference_gene_dna, 'r'):
        sequences_to_align.write(line)
    sequences_to_align.close()

    print(name + '-' + str(xis_count) + ':gene sections extracted')

    # pass on the alignment file to be aligned
    align_gene_sections(name, xis_dict, name + '-' + str(xis_count) + '_gene_sections_w_reference.fna')


def align_gene_sections(name, xis_dict, sequences_to_align):
    xis_count = xis_dict['count']

    out_file = name + '-' + str(xis_count) + '_aligned.fasta'
    clustalomega_cline = ClustalOmegaCommandline(infile=sequences_to_align, outfile=out_file, verbose=True, auto=True)
    print(clustalomega_cline)
    stdout, stderr = clustalomega_cline()
    print(name + '-' + str(xis_count) + ':reference gene & gene sections- alignment complete')

    dic = {}
    alignment = AlignIO.read(out_file, "fasta")
    for i in range(1, alignment.get_alignment_length() + 1):
        dic[i] = {}
    seq_count = 0
    for record in alignment:
        seq_count += 1
        base_count = 0
        if seq_count == 3:
            ref_base_count = 0
            for base in record.seq:
                base_count += 1
                dic[base_count][seq_count] = base
                if base != '-':
                    ref_base_count += 1
                dic[base_count]['ref_position'] = ref_base_count
                max_ref_position = ref_base_count
        else:
            seq_base_count = int(record.description.split(':')[1]) - 1
            for base in record.seq:
                base_count += 1
                dic[base_count][seq_count] = base
                if base != '-':
                    seq_base_count += 1
                dic[base_count]['seq' + str(seq_count) + '_position'] = seq_base_count

    #create a list of absolute positions of alignment where 2 gene sections match
    equal_positions = []
    for position in sorted(dic):
        if dic[position][1] == dic[position][2] and dic[position][1] != '-':
            #print(str(position) + '\t' + str(dic[position]))
            equal_positions.append(position)

    #get each string of matching sequences
    #OPTION1
    flag = 'match'
    matching_dict = {}
    mismatching_dict = {}
    count = 1
    start_position = equal_positions[0]
    matching_dict[count] = {}
    matching_dict[count]['sequence'] = []
    matching_dict[count]['start_position'] = start_position
    mismatching_REF = []
    for position in range(equal_positions[0], equal_positions[-1]+1):
        if dic[position][1] == dic[position][2]:
            test = 'match'
        else:
            test = 'mismatch'
            mismatching_REF.append(dic[position][3])
        if test == flag:
            if test == 'match':
                matching_dict[count]['sequence'].append(dic[position][1])
            else:
                mismatching_dict['seq' + str(1) + '_' + 'stretch' + str(count)].append(dic[position][1])
                mismatching_dict['seq' + str(2) + '_' + 'stretch' + str(count)].append(dic[position][2])
        else:
            if test == 'match':
                count += 1
                #matching_dict[count] = []
                #matching_dict[count].append(dic[position][1])
                matching_dict[count] = {}
                matching_dict[count]['start_position'] = position
                matching_dict[count]['sequence'] = []
                matching_dict[count]['sequence'].append(dic[position][1])
            else:
                mismatching_dict['seq' + str(1) + '_' + 'stretch' + str(count)] = []
                mismatching_dict['seq' + str(2) + '_' + 'stretch' + str(count)] = []
                mismatching_dict['seq' + str(1) + '_' + 'stretch' + str(count)].append(dic[position][1])
                mismatching_dict['seq' + str(2) + '_' + 'stretch' + str(count)].append(dic[position][2])
            flag = test

    #for each possible direct repeat, get the gene sequence before & after, compare to reference
    uber_dict = {}
    matching_ratios = []
    for i in range(1,count + 1):
        uber_dict[i] = {}
        uber_dict[i]['direct_repeat'] = ''.join(matching_dict[i]['sequence'])
        uber_dict[i]['start_position'] = matching_dict[i]['start_position']
        temp_seq_list = []
        for n in range(1, count):
            if n < i:
                for item in mismatching_dict['seq1_stretch' + str(n)]:
                    temp_seq_list.append(item)
            else:
                for item in mismatching_dict['seq2_stretch' + str(n)]:
                    temp_seq_list.append(item)
        uber_dict[i]['gene_seq'] = ''.join(temp_seq_list)
        uber_dict[i]['ref_match_ratio'] = SequenceMatcher(None, uber_dict[i]['gene_seq'], ''.join(mismatching_REF)).ratio()
        matching_ratios.append(uber_dict[i]['ref_match_ratio'])

    #pick the possible direct repeat that has the best corresponding gene match to reference
    #NO PLAN IF 2 DIR_REPEATS HAVE MAX RATIOS & SAME LENGTH
    second_count = 0
    narrowed_dict = {}
    direct_repeat_lengths = []
    for i in range(1,count + 1):
        if uber_dict[i]['ref_match_ratio'] == max(matching_ratios):
            second_count += 1
            absolute_start_position = uber_dict[i]['start_position']
            ref_gene_start_position = dic[uber_dict[i]['start_position']]['ref_position']
            direct_repeat = uber_dict[i]['direct_repeat']
            direct_repeat_length = len(uber_dict[i]['direct_repeat'])
            direct_repeat_lengths.append(direct_repeat_length)
            #UNSURE IF I NEED TO ADD ONE TO THE 1ST KEY
            #ALSO NEED TO ADD THIS TO THE 2ND CALCULATION AFTER IF SECOND_COUNT...
            element_start_positon = dic[uber_dict[i]['start_position']]['seq1_position'] + 1
            element_end_positon = dic[uber_dict[i]['start_position']]['seq2_position']
            narrowed_dict[i] = uber_dict[i]
    if second_count > 1:
        for i in narrowed_dict.keys():
            if len(narrowed_dict[i]['direct_repeat']) == max(direct_repeat_lengths):
                absolute_start_position = narrowed_dict[i]['start_position']
                ref_gene_start_position = dic[narrowed_dict[i]['start_position']]['ref_position']
                direct_repeat = narrowed_dict[i]['direct_repeat']

    xis_dict['direct repeat'] = direct_repeat
    xis_dict['interruption position'] = ref_gene_start_position
    xis_dict['element start position'] = element_start_positon
    xis_dict['element end position'] = element_end_positon
    xis_dict['element length'] = element_end_positon - element_start_positon + 1

    start_to_start = abs(element_start_positon - xis_dict['xis coordinates'][0])
    end_to_end = abs(element_end_positon - xis_dict['xis coordinates'][1])
    xis_dict['xis to edge'] = (min(start_to_start, end_to_end))
    if start_to_start == xis_dict['xis to edge']:
        xis_dict['xis location'] = 'begin'
    else:
        xis_dict['xis location'] = 'end'

    genome = open(args.genome, 'r')
    element_cut_coordinates = (element_start_positon, element_end_positon)
    element_file = open(name + '-' + str(xis_count) + '_interruption_element.fna', 'w')
    element_sequence = sequence_cutter(genome, xis_dict['contig accession'], element_cut_coordinates)
    element_file.write(element_sequence + '\n')
    element_file.close()
    genome.close()

    genome = open(args.genome, 'r')
    sequence = sequence_cutter(genome, xis_dict['contig accession'], element_cut_coordinates, 'FALSE')
    xis_dict['element GC'] = (sequence.count('G') + sequence.count('C'))/len(sequence)
    genome.close()

    '''
    #OPTION2
    #get the two overlapping gene sections & find the longest matching segment
    seq_1 = []
    seq_2 = []
    for position in range(equal_positions[0], equal_positions[-1]+1):
        seq_1.append(dic[position][1])
        seq_2.append(dic[position][2])
    length = equal_positions[-1]+1 - equal_positions[0]
    matching = SequenceMatcher(None, ''.join(seq_1), ''.join(seq_2)).find_longest_match(0, length, 0, length)
    start_position = matching[1] + equal_positions[0]
    length = matching[2]
    dir_rep_list = []
    for position in range(start_position, start_position + length):
        dir_rep_list.append(dic[position][1])
    print('direct repeat:' + ''.join(dir_rep_list))
    print('start position:' + str(start_position))
    '''

    #collect bases from each gene section before/after direct repeats to complete the gene
    full_gene_list = []
    for position in range(1,absolute_start_position):
        if dic[position]['ref_position'] != 0:
            full_gene_list.append(dic[position][1])
    for position in range(absolute_start_position, alignment.get_alignment_length() + 1):
        full_gene_list.append(dic[position][2])
        if dic[position]['ref_position'] == max_ref_position:
            break
    full_gene = ''.join(full_gene_list)
    full_gene_file = open(name + '-' + str(xis_count) + '_complete_gene.fna', 'w')
    full_gene_file.write('>' + name + '-' + str(xis_count) + ' full gene sequence' + '\n' + full_gene)
    full_gene_file.close()

    print('-----------------')
    print(args.name)
    print(xis_dict['count'])
    print(xis_dict['reference gene'])
    print(xis_dict['direct repeat'])
    print(xis_dict['interruption position'])
    print(xis_dict['class'])
    print(xis_dict['xis coordinates'][0]) #NEED TO CHECK THE XIS START
    print(xis_dict['xis coordinates'][1]) #NEED TO CHECK THE XIS END
    print(xis_dict['xis to edge']) # THIS IS OFF, NEED TO CHECK THE XIS COORDINATES VS WHAT I FOUND IN PUB.
    print(xis_dict['xis location'])
    #print xis direction
    print(xis_dict['element length'])
    print(xis_dict['element GC']) # THIS IS OFF, NEED TO CHECK ELEMENT SEQUENCE
    print(xis_dict['contig accession'])
    print(xis_dict['element start position'])
    print(xis_dict['element end position'])

    # need to print results overview to a file:

    # ADD TO PRINTED UPDATES 'COMPLETE' OR SOMETHING

if __name__ == '__main__':
    args = getArgs()
    find_xis_candidates(args.name)
