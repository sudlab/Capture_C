import pysam
import CGAT.Bed as Bed
import collections
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
from CGATPipelines.Pipeline import cluster_runnable


@cluster_runnable
def quantifyAnomolies(bam, probes, outfile):

    bamfile = pysam.AlignmentFile(bam)
    results = dict()
    mapped = bamfile.mapped
    total = 0
    seeks = 0
    for probe in Bed.iterator(IOTools.openFile(probes)):
        
        c = collections.Counter()
 
        for read in bamfile.fetch(probe.contig, probe.start, probe.end,
                                  multiple_iterators=True):

            if read.is_unmapped:
                continue

            c["total"] += 1
    
            if not (read.is_secondary or read.is_supplementary):
                c["primary"] += 1
            else:
                continue

            if read.pos < (probe.start-4) or read.aend > (probe.end +4):
                c["undigested"] += 1

            if read.is_read1:
                c["read1"] += 1

            if not read.mate_is_unmapped:
                c["paired"] += 1

                if read.is_read1 and (
                        read.mpos >= probe.start and read.mpos <= probe.end):
                        c["self_lig"] += 1

            if (total + c["total"] % 10000) == 0:
                E.debug("%s/%s done" % (total + c["total"], mapped))

        E.debug("%s processed, %i found" % (probe.name, c["total"]))
        results[probe.name] = c
        total += c["total"]

    headers = ["Probe", "total", "primary", "undigested",
               "read1", "paired", "self_lig"]

    with IOTools.openFile(outfile, "w") as outf:

        outf.write("\t".join(headers) + "\n")
        
        for probe in results:
            outf.write("\t".join(
                [probe]+[str(results[probe][col])
                         for col in headers[1:]]) + "\n")


@cluster_runnable
def fetchProbeFragments(probe_bed, digest_bed, outfile,
                        lookup_out):

    digest_fragments = pysam.TabixFile(digest_bed)
    bed = Bed.Bed()
    with IOTools.openFile(outfile, "w") as outf, \
         IOTools.openFile(lookup_out,"w") as lookup:

        lookup.write("probe\tfragment\n")
        for probe in Bed.iterator(IOTools.openFile(probe_bed)):
            
            frag = digest_fragments.fetch(probe.contig,
                                          probe.start,
                                          probe.end,
                                          parser=pysam.asBed())
            frag = list(frag)
            if not len(frag) == 1:
                E.warn("%i fragments found for probe %s, skipping" %
                       (len(frag), probe.name))
                continue

            frag = frag[0]
            bed.start = frag.start
            bed.end = frag.end
            bed.contig = frag.contig
            bed["name"] = probe.name
            bed["score"] = "."
            bed["strand"] = "+"

            lookup.write("%s\t%s\n" % (probe.name, frag.name))
            outf.write(str(bed) + "\n")


@cluster_runnable
def sites2fragments(infile, genomefile, outfile):
    '''Convert bedfile of deigestion sites into bedfile of fragments'''

    contig_lengths = {line.split()[0]: int(line.split()[1][:-1])
                      for line in IOTools.openFile(genomefile)}

    last_end = 0
    last_contig = None
    name = 0
    new_bed = Bed.Bed()
    new_bed["strand"] = "+"
    new_bed["score"] = "."
    with IOTools.openFile(outfile, "w") as outf:
        for bed in Bed.iterator(IOTools.openFile(infile)):
       
            if last_contig is not None and not bed.contig == last_contig:
                name += 1
                new_bed.start = last_end
                new_bed.contig = last_contig
                new_bed.end = contig_lengths[bed.contig]
                new_bed["name"] = str(name)

                outf.write(str(new_bed) + "\n")
        
                last_end = 0
                
            last_contig = bed.contig
            new_bed.contig = last_contig
            new_bed.start = last_end
            new_bed.end = bed.start
            name += 1
            new_bed["name"] = str(name)
            outf.write(str(new_bed) + "\n")
            last_end = bed.end

        name += 1
        new_bed.start = last_end
        new_bed.contig = last_contig
        new_bed.end = contig_lengths[bed.contig]
        new_bed["name"] = str(name)

        outf.write(str(new_bed) + "\n")

    pysam.tabix_index(outfile, force=True, preset="bed")


def readBundle(bamfile):

    last_name = None
    bundle = []

    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        if last_name is not None and not last_name == read.query_name:
            yield bundle
            bundle = []

        bundle.append(read)
        last_name = read.query_name


@cluster_runnable
def countInteractions(reads, digest, probe_fragments, outfile,
                      metrics_file):

    reads = pysam.AlignmentFile(reads)
    digest_intervals = Bed.readAndIndex(
        IOTools.openFile(digest), with_values="name")
    probe_fragments = Bed.readAndIndex(
        IOTools.openFile(probe_fragments), with_values="name")
    c = collections.Counter()
    results = collections.defaultdict(int)
    contigs = digest_intervals.keys()
    rejects=pysam.AlignmentFile(outfile + ".rejects.bam", "wb", template=reads)

    E.debug("Starting Counting Run...")
    for fragment in readBundle(reads):

        if c["fragments"] % 1000 == 0:
            E.debug("\t".join(["%s:\t%s" % x for x in c.iteritems()]))

        c["fragments"] += 1
        c["reads"] += len(fragment)

        def _get_read_position(aln):

            if aln.is_reverse:
                c = aln.cigartuples[::-1]
            else:
                c = aln.cigartuples

            pos = 0
            for op, length in c:
                if op == 0:
                    return pos
                else:
                    pos += length

            return pos

        def _get_first(alignments):
            '''find alignment that maps earliest part of read'''
            slist = sorted(alignments, key=lambda x: _get_read_position(x))
            return (slist[0], slist[1:])
            
        def _get_probes(read):
            try:
                m = list(
                    probe_fragments[reads.getrname(read.reference_id)].find(
                        read.pos, read.aend))
            except KeyError:
                return []
            return m

        first_read = [aln for aln in fragment if aln.is_read1]
        second_read = [aln for aln in fragment if aln.is_read2]

        if len(first_read) == 0 or len(second_read) == 0:
            c["Unpaired"] += 1
            continue

        assert len(first_read) + len(second_read) == len(fragment)

        primary_aln, other_aln = zip(_get_first(first_read),
                                     _get_first(second_read))
        other_aln = sum(other_aln, [])
        
        probe_matches = [_get_probes(read) for read in primary_aln]

        if len(sum(probe_matches, [])) == 0:
            c["no_probe"] += 1
            for read in fragment:
                rejects.write(read)
            continue
    
        other_matches = set(sum([
            list(digest_intervals[reads.getrname(read.reference_id)].find(
                read.pos, read.aend))
            for read in other_aln
            if reads.getrname(read.reference_id) in contigs],[]))

        primary_matches = set(sum([
            list(digest_intervals[reads.getrname(read.reference_id)].find(
                read.pos, read.aend))
            for read in primary_aln
            if reads.getrname(read.reference_id) in contigs],[]))

        if len(primary_matches) > 2:
            c["multi-hit"] += 1
            continue
            
        if not all([match in primary_matches for match in other_matches]):
            c["multi-hit"] += 1
            continue
        
        primary_matches = list(primary_matches)
        if len(primary_matches) == 2:
            results[(primary_matches[0][2], primary_matches[1][2])] += 1
            results[(primary_matches[1][2], primary_matches[0][2])] += 1
        elif len(primary_matches) == 1:
            results[(primary_matches[0][2], primary_matches[0][2])] += 1
        else:
            raise ValueError("Matches not found")

    with IOTools.openFile(outfile, "w") as outf:

        outf.write("Frag1\tFrag2\tCount\n")
        for pair in results:
            outf.write("%s\t%s\t%s\n" % (pair+(results[pair],)))

    with IOTools.openFile(metrics_file, "w") as outf:
        
        outf.write("\t".join(c.keys())+"\n")
        outf.write("\t".join(map(str,c.values())) + "\n")

    rejects.close()


def interactions2BedGraph(track, probe, db, outfile, chrom=None, start=None, end=None,
                          region_size=1000000, window=2000, step=200):

    import CGATPipelines.PipelineUtilities as PUtils
    import pandas
    import numpy as np

    E.debug("Get Probe fragment")
    #get probe fragment
    statement = '''SELECT chr, start, end, fragment
                   FROM probe_fragments 
                   INNER JOIN probe_fragments_lookup as lu
                   ON probe_fragments.name = lu.probe
                   WHERE name = '%(probe)s' '''

    probes = PUtils.fetch_DataFrame(statement % locals(), database=db)
    probe_frag = probes.fragment.iloc[0]

    if chrom is None:
        probe_chr = probes.chr.iloc[0]

    probe_centre = 0.5*(probes.start.iloc[0] + probes.end.iloc[0])

    if start is None:
        start = probe_centre - region_size/2

    if end is None:
        end = probe_centre + region_size/2

    E.debug("Fetch data")
    statement = '''SELECT (start+end)/2 as centre,
                          count
                   FROM interaction_counts as ic
                   INNER JOIN fragments 
                     ON fragments.name = ic.Frag2
                   WHERE abs(%(probe_frag)s - ic.Frag2) > 1
                    AND Frag1 = '%(probe_frag)s' AND track = '%(track)s'
                    AND chr='%(probe_chr)s'
                    AND centre > %(start)s AND centre <= %(end)s
                    AND NOT ABS(centre - %(probe_start)) < 2000
                   ORDER BY centre '''

    E.debug(statement % locals())
    interactions = PUtils.fetch_DataFrame(statement % locals(), db)
    E.debug("Got %i interacting fragments" % interactions.shape[0])

    interactions["centre"] = interactions["centre"].astype("float64")
    interactions.set_index("centre", inplace=True)
    interactions = interactions.sort_index()

    window_starts = pandas.Series(np.arange(start, end, step))
    window_starts.index = window_starts.values
    window_starts.index.name = "start"

    def _applyToWindow(val):

        #E.debug(val)
        return interactions.loc[(val-window/2):(val+window/2)]["Count"].sum()

    E.debug("rolling")
    windowed = window_starts.apply(_applyToWindow)
    assert windowed.sum() >= interactions["Count"].sum(), \
        "windowed sum %s, total count %s" % (windowed.sum(), interactions["Count"].sum())

    E.debug("format output")
    windowed = windowed.reset_index()

    windowed["end"] = windowed["start"] + step
    windowed["chr"] = probe_chr
    windowed = windowed[["chr", "start", "end", 0]]
    windowed.to_csv(IOTools.openFile(outfile, "w"), sep = "\t",
                    index=False, header=False)

    

