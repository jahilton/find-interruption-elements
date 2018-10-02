# I THINK COMPLETED GENE HAS 3 TOO MANY BP AT FRONT - BUT THAT IS HOW IT IS ALIGNED
# MAYBE INCLUDE REFERENCE GENE IN THE HEADER OF THE COMPLETED GENE
# NEED TO FIGURE OUT WHAT IS GOING ON IN ALIGNMENT ETC WITH XIS GENES THAT ARENT ON INT ELE & HOW THEY ARE REPORTED
# IF MAKEBLASTDB NEEDS TO BE RUN, CREATE A HUMAN-FRIENDLY ERROR MESSAGE

# install ncbi tools
# pull nr db if local blasting
# install clustal-omega (install argstable for this)
# $ export PATH=$PATH:/Users/jason/code/find-interruption-elements-side/NCBI_Resources/ncbi-blast-2.4.0+/bin/ ----- how to get this permanent?
# first requires the db to go through makeblastdb - can i do that here in python? (make it happen only if it isn't already made)
# makeblastdb -in genome -dbtype nucl -parse_seqids
# makeblastdb -in known_interrupted_genes_protein.faa -dbtype prot -parse_seqids
# and for nr

# how to get concat_cyano.genes.fna/faa (AND NEED TO CHANGE THAT NAME, maybe UDPATE DB)
# Pulled all IMG JGI cyano genome assemblies, got the IMG_nnnnn ID (https://img.jgi.doe.gov/cgi-bin/w/main.cgi)
# Crossed that with all available downloadable projects at https://genome.jgi.doe.gov/portal/
# curl 'https://genome.jgi.doe.gov/portal/IMG_2675903261/download/download_bundle.tar.gz' -b cookies > download_2675903261.tar.gz ; gunzip download_2675903261.tar.gz ; tar -xvf download_2675903261.tar ; mv 2675903261/2675903261.genes* cyano_genes/ ; rm -r 2675903261/ ; rm download_2675903261.tar
# above command for 754 projects: downloads tar bundle, unzips it, extracts the genes faa & fna file, deletes the rest
# then i concatenated 754 fna files to 1 file & 754 faa files to 1 file


from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast import NCBIXML
from difflib import SequenceMatcher
import argparse
from operator import itemgetter

EPILOG = '''
    Write something here maybe
    '''

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
    'Ga0055315_143016': 'tyrosine',
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
    'Ava_B0242': 'transposase (Ava_B0242)',
    'Aazo_1350': 'nifE',
    'Aazo_1352': 'nifK',
    'Aazo_1353': 'nifD',
    'Aazo_1354': 'nifH',
    'Aazo_1357': 'fdxN',
    'Aazo_1358': 'NifB',
    'Aazo_2640': 'coxA',
    'Aazo_2682': 'integrase family protein (Aazo_2682)',
    'Aazo_3865': 'hupL',
    'Aazo_3866': 'hupS',
    'Aazo_3917': 'hglE',
    'Aazo_4140': 'flv3B',
    'Aazo_5221': 'NADPH-dependent FMN reductase (Aazo_5221)',
    'Cal7507_0570': 'phospholipase D/transphosphatidylase (Cal7507_0570)',
    'Cal7507_5433': 'nifJ',
    'Cal7507_5656': 'FAD-dependent oxidoreductase (Cal7507_5656)',
    'Ana7108_2845': 'primase P4 (Ana7108_2845)',
    'Pse6802_0098': 'ATP-dependent DNA helicase (Pse6802_0098)',
    'Pse6802_3453': 'predicted integral membrane protein (Pse6802_3453)',
    'Cal7103DRAFT_00047390': 'caspase domain-containing protein (Cal7103DRAFT_00047390)',
    'CylstDRAFT_1988': 'hypothetical protein (CylstDRAFT_1988)',
    'Mic7126DRAFT_5075': 'arabinose efflux permease (Mic7126DRAFT_5075)',
    'Mas10914DRAFT_5058': 'transposase (ISSoc8) (Mas10914DRAFT_5058)'
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


def get_xis_candidate_structure():
    return {
        'xis orientation': '',
        'class': '',
        'element GC': '',
        'start': '',
        'reference gene': '',
        'interruption position': '',
        'element length': '',
        'xis coordinates': ['', ''],
        'element coordinates': ['', ''],
        'xis_orientation_on_contig': '',
        'score': '',
        'end': '',
        'xis to edge': '',
        'evalue_to_known_xis': '',
        'name': '',
        'xis location': '',
        'contig accession': '',
        'distal_gene_orientation': '',
        'proximal_gene_orientation': '',
        'count': '',
        'top_xis_hit': '',
        'direct repeat': '',
        'notes': ''
    }


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
    tblastn_cline = NcbitblastnCommandline(query=args.xis_set, db=args.genome, evalue=1e-20, outfmt=5, out=name + '_xis_tn_genome.xml')
    print(str(tblastn_cline))
    stdout, stderr = tblastn_cline()
    print(name + ':xis_set-genome BLAST complete')

    # from BLAST results, pull hits
    blast_results = open(name + '_xis_tn_genome.xml', 'r')
    blast_records = NCBIXML.parse(blast_results)
    list_of_hits_dicts = []
    # build a dictionary of all the hits to be sorted
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                xis_dict = get_xis_candidate_structure()
                xis_dict['contig accession'] = alignment.accession
                xis_dict['start'] = hsp.sbjct_start
                xis_dict['end'] = hsp.sbjct_end
                xis_dict['evalue_to_known_xis'] = hsp.expect
                xis_dict['top_xis_hit'] = blast_record.query
                xis_dict['score'] = hsp.score
                if (hsp.frame[0] < 0 and hsp.frame[1] < 0) or (hsp.frame[0] > -1 and hsp.frame[1] > -1):
                    xis_dict['xis_orientation_on_contig'] = 'plus'
                else:
                    xis_dict['xis_orientation_on_contig'] = 'minus'
                list_of_hits_dicts.append(xis_dict)
    blast_results.close()

    # sort blast hits and pull one per genome region
    off_limits_ranges = {}
    list_of_xis_dicts = []
    xis_count = 0
    sorted_by_score = sorted(list_of_hits_dicts, key=itemgetter('score'), reverse=True)
    sorted_by_evalue = sorted(sorted_by_score, key=itemgetter('evalue_to_known_xis'), reverse=False)
    for xis_dict in sorted_by_evalue:
        if xis_dict['contig accession'] not in off_limits_ranges.keys():
            off_limits_ranges[alignment.accession] = []
        middle_of_hit = round((xis_dict['start']+xis_dict['end'])/2, 0)
        if middle_of_hit not in off_limits_ranges[xis_dict['contig accession']]:
            xis_count += 1
            xis_dict['count'] = xis_count
            for i in range(xis_dict['start'], xis_dict['end']):
                off_limits_ranges[xis_dict['contig accession']].append(i)

            # cut the xis_candidates from the contigs based on BLAST results
            contig_title = xis_dict['contig accession']
            genome = open(args.genome, 'r')
            contig_sequence = sequence_cutter(genome, contig_title, 'all', 'FALSE')
            xis_coordinates = orf_finder(contig_sequence, xis_dict['start'], xis_dict['end'])
            genome.close()

            # log characteristics of the xis candidate to be used in later processing and reporting
            xis_dict['name'] = name
            xis_dict['xis coordinates'] = xis_coordinates
            locus_tag = xis_dict['top_xis_hit'].split()[1]
            xis_dict['top_xis_hit'] = locus_tag
            xis_dict['class'] = recombinase_class[locus_tag]
            list_of_xis_dicts.append(xis_dict)

            genome = open(args.genome, 'r')
            xis_fasta = sequence_cutter(genome, contig_title, xis_coordinates)
            xis_candidate_sequence = open(name + '-' + str(xis_count) + '_xis_candidate.fna', 'w')
            xis_candidate_sequence.write(xis_fasta + '\n')
            xis_candidate_sequence.close()
            genome.close()

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
            xis_flank = open(name + '-' + str(xis_count) + '_xis_flank.fna', 'w')
            xis_flank.write(flank_fasta + '\n')
            xis_flank.close()
            genome.close()

    print(name + ':' + str(xis_count) + ' xis candidates identified')

    # if no hits were found, then quit the program
    if xis_count == 0:
        return

    # pass on the xis_candidates and flanking regions to search for interrupted genes in the flanking regions
    for xis_dict in list_of_xis_dicts:
        find_interrupted_gene(xis_dict, name + '-' + str(xis_dict['count']) + '_xis_flank.fna')


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
                return(header_line + ':' + str(cut_coordinates[0]) + ':' + str(cut_coordinates[1]) + '\n' + ''.join(sequence[start - 1:stop]))
        elif header == 'FALSE':
            if cut_coordinates == 'all':
                return(''.join(sequence[:]))
            else:
                start, stop = cut_coordinates
                return(''.join(sequence[start - 1:stop]))
    else:
        print('sequence name not unique')


def orf_finder(contig_sequence, blast_hit_start, blast_hit_end):
    '''
    Given a sequence and the start and stop of a BLAST hit, returns the ORF that the sequence range is a part of.

    Takes the start coordinate of the BLAST hit and traces the sequence back each codon at a time until it finds a start codon.
    Take the end coordinate of the BLAST hit and traces the sequence forward each codon at a time until it finds a stop codon.
    Returns the predicted ORF sequence.
    '''

    # HOW DOES THIS WORK IF FRAME IS NEGATIVE?
    start_codon = 'ATG'
    start_flag = 'FALSE'
    start_position = blast_hit_start - 1
    while start_flag == 'FALSE':
        codon_to_check = contig_sequence[start_position:start_position+3]
        if codon_to_check == start_codon:
            start_flag = 'TRUE'
            orf_start_position = start_position + 1
        start_position -= 3

    stop_codons = ['TAA', 'TAG', 'TGA']
    stop_flag = 'FALSE'
    stop_position = blast_hit_end
    while stop_flag == 'FALSE':
        codon_to_check = contig_sequence[stop_position:stop_position+3]
        if codon_to_check in stop_codons:
            stop_flag = 'TRUE'
            orf_stop_position = stop_position + 3
        stop_position += 3

    return orf_start_position, orf_stop_position


def find_interrupted_gene(xis_dict, xis_plus_flank):
    xis_count = xis_dict['count']
    name = xis_dict['name']

    # BLAST the xis + flanking region against the known reference genes for previously found interrupted genes
    blastx_cline = NcbiblastxCommandline(query=xis_plus_flank, db='interrupted_genes/known_interrupted_genes_protein.faa', culling_limit=15, evalue=1e-15, outfmt=5, out=name + '-' + str(xis_count) + '_xisflank_x_KnownIntGenes.xml')
    print(str(blastx_cline))
    stdout, stderr = blastx_cline()
    print(name + '-' + str(xis_count) + ':xis w/ flank - known Interrupted Genes BLAST complete')

    strand = {}
    sortme = {}
    off_limits_range = []
    xis_locus_tag = {}

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
                    xis_locus_tag[alignment.accession] = alignment.hit_def.split()[0]
                    if (hsp.frame[0] < 0 and hsp.frame[1] < 0) or (hsp.frame[0] > -1 and hsp.frame[1] > -1):
                        strand[alignment.accession] = 'plus'
                    else:
                        strand[alignment.accession] = 'minus'
                    for i in range(hsp.query_start, hsp.query_end):
                        off_limits_range.append(i)
    blast_results.close()

    # if interrupted gene isn't a previously known gene, BLAST against a larger protein database
    if sortme == {}:  # AND SOME OTHER CRITERIA ABOUT THE BLAST RESULTS?
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
                        if (hsp.frame[0] < 0 and hsp.frame[1] < 0) or (hsp.frame[0] > -1 and hsp.frame[1] > -1):
                            strand[alignment.accession] = 'plus'
                        else:
                            strand[alignment.accession] = 'minus'
                        for i in range(hsp.query_start, hsp.query_end):
                            off_limits_range.append(i)
        blast_results.close()

        if sortme == {}:
            print(name + '-' + str(xis_count) + ':No reference genes found in flanking region')
            xis_dict['notes'] = 'No reference genes found in flanking region'
            reporter(xis_dict)
            return

        # pull the hit from above with the lowest % coverage
        lowest_coverage = sorted(sortme.values())[0]
        for key in sortme.keys():
            if sortme[key] == lowest_coverage:
                lowest_coverage_accession = key

        # use the accession of the hit to pull the full gene sequence
        reference_gene_aa_file = open(name + '-' + str(xis_count) + '_reference_gene.faa', 'w')
        reference_gene_aa_file.write(sequence_cutter(open('/Users/jason/code/find-interruption-elements-side/cyano_genes/concat_cyano_genes.faa', 'r'), lowest_coverage_accession, 'all'))
        reference_gene_aa_file.close()

        reference_gene_dna_file = open(name + '-' + str(xis_count) + '_reference_gene.fna', 'w')
        reference_gene_dna_file.write(sequence_cutter(open('/Users/jason/code/find-interruption-elements-side/cyano_genes/concat_cyano_genes.fna', 'r'), lowest_coverage_accession, 'all'))
        reference_gene_dna_file.close()

        # try this
        xis_dict['reference gene'] = lowest_coverage_accession
        print(name + '-' + str(xis_count) + ':reference gene identified: ' + lowest_coverage_accession)

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

        xis_dict['reference gene'] = reference_genes[xis_locus_tag[lowest_coverage_accession]]
        print(name + '-' + str(xis_count) + ':reference gene identified: ' + reference_genes[xis_locus_tag[lowest_coverage_accession]])

    xis_dict['proximal_gene_orientation'] = strand[lowest_coverage_accession]

    if xis_dict['proximal_gene_orientation'] == xis_dict['xis_orientation_on_contig']:
        xis_dict['xis orientation'] = 'plus'
    else:
        xis_dict['xis orientation'] = 'minus'

    # pass on the reference gene to search the whole genome for all gene regions
    find_all_gene_regions(xis_dict, name + '-' + str(xis_count) + '_reference_gene.faa', name + '-' + str(xis_count) + '_reference_gene.fna')


def find_all_gene_regions(xis_dict, reference_gene_aa, reference_gene_dna):
    name = xis_dict['name']
    xis_count = xis_dict['count']
    xis_coordinates = xis_dict['xis coordinates']
    xis_contig = xis_dict['contig accession']

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
            hitsd['hit' + str(count)]['contig_title'] = alignment.accession
            hitsd['hit' + str(count)]['contig_start'] = hsp.sbjct_start
            hitsd['hit' + str(count)]['contig_end'] = hsp.sbjct_end
            if (hsp.frame[0] < 0 and hsp.frame[1] < 0) or (hsp.frame[0] > -1 and hsp.frame[1] > -1):
                hitsd['hit' + str(count)]['strand'] = 'plus'
            else:
                hitsd['hit' + str(count)]['strand'] = 'minus'
            if hsp.sbjct_start < xis_coordinates[0]:
                hitsd['hit' + str(count)]['gene_location'] = 'before'
            else:
                hitsd['hit' + str(count)]['gene_location'] = 'after'

    blast_results.close()

    # if only 1 gene section found in the genome, quit the program
    if count == 1:
        print(name + '-' + str(xis_count) + ':only 1 gene section located in genome')
        xis_dict['notes'] = 'only 1 gene section located in genome'
        reporter(xis_dict)
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

    if (hitsd[nearest_gene_section]['contig_start'] in range(hitsd[second_gene_section]['contig_start'], hitsd[second_gene_section]['contig_end'])
        or hitsd[nearest_gene_section]['contig_end'] in range(hitsd[second_gene_section]['contig_start'], hitsd[second_gene_section]['contig_end'])):
        print(name + '-' + str(xis_count) + ':2 gene sections were found to overlap')
        xis_dict['notes'] = '2 gene sections were found to overlap'
        reporter(xis_dict)
        return

    xis_dict['distal_gene_orientation'] = hitsd[second_gene_section]['strand']

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
    align_gene_sections(xis_dict, name + '-' + str(xis_count) + '_gene_sections_w_reference.fna')


def align_gene_sections(xis_dict, sequences_to_align):
    name = xis_dict['name']
    xis_count = xis_dict['count']

    # FOR CLUSTO, --residuenumber PRINTS THE RESIDUE NUMBER, MAY HELP
    # VERBOSE OPTION MAY NOT BE NEEDED
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

    # create a list of absolute positions of alignment where 2 gene sections match
    equal_positions = []
    for position in sorted(dic):
        if dic[position][1] == dic[position][2] and dic[position][1] != '-':
            equal_positions.append(position)

    # get each string of matching sequences
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

    # for each possible direct repeat, get the gene sequence before & after, compare to reference
    uber_dict = {}
    matching_ratios = []
    for i in range(1, count + 1):
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

    # pick the possible direct repeat that has the best corresponding gene match to reference
    # NO PLAN IF 2 DIR_REPEATS HAVE MAX RATIOS & SAME LENGTH
    second_count = 0
    narrowed_dict = {}
    direct_repeat_lengths = []
    for i in range(1, count + 1):
        if uber_dict[i]['ref_match_ratio'] == max(matching_ratios):
            second_count += 1
            absolute_start_position = uber_dict[i]['start_position']
            ref_gene_start_position = dic[uber_dict[i]['start_position']]['ref_position']
            direct_repeat = uber_dict[i]['direct_repeat']
            direct_repeat_length = len(uber_dict[i]['direct_repeat'])
            direct_repeat_lengths.append(direct_repeat_length)
            # UNSURE IF I NEED TO ADD ONE TO THE 1ST KEY
            # ALSO NEED TO ADD THIS TO THE 2ND CALCULATION AFTER IF SECOND_COUNT...
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
    xis_dict['element coordinates'] = (element_start_positon, element_end_positon)
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
    element_file.write(sequence_cutter(genome, xis_dict['contig accession'], element_cut_coordinates) + '\n')
    element_file.close()
    genome.close()

    genome = open(args.genome, 'r')
    sequence = sequence_cutter(genome, xis_dict['contig accession'], element_cut_coordinates, 'FALSE')
    if len(sequence) == 0:
        xis_dict['element GC'] = 0
    else:
        xis_dict['element GC'] = round(((sequence.count('G') + sequence.count('C'))/len(sequence))*100,2)
    genome.close()

    # collect bases from each gene section before/after direct repeats to complete the gene
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

    print(name + '-' + str(xis_count) + ':direct repeat identified & interrupted gene completed')

    reporter(xis_dict)


def reporter(xis_dict):
    genome_file = args.genome

    try:
        open('xis_summary.txt', 'r')
    except FileNotFoundError:
        outfile = open('xis_summary.txt', 'w')
        outfile.write('Genome Name' + '\t' + 'Genome File' + '\t' + 'xis Number' + '\t' + 'Interrupted Gene' + '\t' + 'Contig Accession' + '\t' + 'Top xis Hit' + '\t' + 'Top xis Hit evalue' + '\t' + 'xis Class' + '\t' + 'xis Start' + '\t' + 'xis Stop' + '\t' + 'xis Orientation' + '\t' + 'Direct Repeat' + '\t' + 'Position Interrupted' + '\t' + 'xis Distance to edge' + '\t' + 'xis Location' + '\t' + 'Element Length' + '\t' + 'Element GC' + '\t' + 'Element Start' + '\t' + 'Element Stop' + '\t' + 'Notes')
        outfile.write('\n')
        outfile.close()
    outfile = open('xis_summary.txt', 'a')
    outfile.write(xis_dict['name'] + '\t' + genome_file + '\t' + '#' + str(xis_dict['count']) + '\t' + xis_dict['reference gene'] + '\t' + xis_dict['contig accession'] + '\t' + xis_dict['top_xis_hit'] + '\t' + str(xis_dict['evalue_to_known_xis']) + '\t' + xis_dict['class'] + '\t' + str(xis_dict['xis coordinates'][0]) + '\t' + str(xis_dict['xis coordinates'][1]) + '\t' + xis_dict['xis orientation'] + '\t' + xis_dict['direct repeat'] + '\t' + str(xis_dict['interruption position']) + '\t' + str(xis_dict['xis to edge']) + '\t' + xis_dict['xis location'] + '\t' + str(xis_dict['element length']) + '\t' + str(xis_dict['element GC']) + '\t' + str(xis_dict['element coordinates'][0]) + '\t' + str(xis_dict['element coordinates'][1]) + '\t' + xis_dict['notes'])
    outfile.write('\n')
    outfile.close()

    print(xis_dict['name'] + '-' + str(xis_dict['count']) + ':resuls writen to xis_summary.txt')


if __name__ == '__main__':
    args = getArgs()
    find_xis_candidates(args.name)
