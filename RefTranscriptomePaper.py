import CGAT.GTF as GTF
import CGAT.IOTools as IOTools

import CGATPipelines.Pipeline as P

def tabulateGENCODETranscriptGTF(infile, info_out, tag_out):
    '''
    extract the required transcript information from the GENCODE GTF
    '''

    info_outf = IOTools.openFile(info_out, 'w')
    tag_outf = IOTools.openFile(tag_out, 'w')

    info_columns = ['transcript_name', 'transcript_status',
                    'transcript_support_level', 'transcript_type',
                    'level', 'protein_id']

    info_outf.write('transcript_id\tgene_id\t%s\n' % '\t'.join(info_columns))

    tag_outf.write('%s\n' % '\t'.join(('transcript_id', 'gene_id', 'tag')))

    gtfs = GTF.track_iterator(IOTools.openFile(infile), repeated_attributes=True)

    for gtf in gtfs:

        attr_dict = gtf.asDict()

        if 'tag' in attr_dict:
            for tag in attr_dict['tag']:
                tag_outf.write('%s\n' % "\t".join(
                    (gtf.transcript_id, gtf.gene_id, tag)))

        values = [gtf.transcript_id, gtf.gene_id]
        for column in info_columns:
            if column in attr_dict:
                values.append(attr_dict[column][0])
            else:
                values.append('')

        info_outf.write('%s\n' % '\t'.join(map(str, values)))


    info_outf.close()
    tag_outf.close()


def tabulateGENCODEGeneGTF(infile, info_out):
    '''
    extract the required gene information from the GENCODE GTF
    '''

    info_outf = IOTools.openFile(info_out, 'w')

    info_columns = ['gene_name', 'gene_status',
                    'gene_type', 'level']
    info_key = 'gene_id'
    info_outf.write('%s\t%s\n' % (info_key, '\t'.join(info_columns)))


    gtfs = GTF.track_iterator(
        IOTools.openFile(infile), repeated_attributes=True, strict=False)

    for gtf in gtfs:

        attr_dict = gtf.asDict()

        values = [gtf.gene_id]
        for column in info_columns:
            if column in attr_dict:
                values.append(attr_dict[column][0])
            else:
                values.append('')

        info_outf.write('%s\n' % '\t'.join(map(str, values)))


    info_outf.close()


def getTranscript2GeneMap(infile, outfile):
    ''' parse a gtf to output a transript to gene map'''

    with IOTools.openFile(outfile, 'w') as outf:

        outf.write("%s\t%s\n" % ("transcript_id", "gene_id"))

        gtfs = GTF.iterator(IOTools.openFile(infile, 'r'))

        observed = set()

        for gtf in gtfs:
            if gtf.transcript_id in observed:
                continue
            else:
                observed.add(gtf.transcript_id)
                outf.write("%s\t%s\n" % (gtf.transcript_id, gtf.gene_id))


def buildJunctions(infile, outfile):
    '''code copied from CGATPipelines.pipeline_mapping.py'''

    outf = IOTools.openFile(outfile, "w")
    njunctions = 0
    for gffs in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile, "r"))):

        gffs.sort(key=lambda x: x.start)
        end = gffs[0].end
        for gff in gffs[1:]:
            # subtract one: these are not open/closed coordinates but
            # the 0-based coordinates
            # of first and last residue that are to be kept (i.e., within the
            # exon).
            outf.write("%s\t%i\t%i\t%s\n" %
                       (gff.contig, end - 1, gff.start, gff.strand))
            end = gff.end
            njunctions += 1

    outf.close()

    # make unique
    statement = '''mv %(outfile)s %(outfile)s.tmp;
                   cat < %(outfile)s.tmp | sort | uniq > %(outfile)s;
                   rm -f %(outfile)s.tmp; '''
    P.run()
