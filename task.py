#!/usr/bin/env python
from __future__ import division, with_statement

"""
task: wraps a GMTK subtask to reduce size of output
"""

__version__ = "$Revision: 5777 $"

# Copyright 2009, 2011 Michael M. Hoffman <mmh1@washington.edu>

from errno import ENOENT
from os import extsep, fdopen
import re
import sys
from tempfile import gettempdir, mkstemp

from genomedata import Genome
from numpy import array, zeros, empty, where, diff, r_
from path import path

from ._util import (BED_SCORE, BED_STRAND, ceildiv, DTYPE_IDENTIFY, EXT_FLOAT,
                    EXT_INT, EXT_LIST, find_segment_starts, get_label_color,
                    _make_continuous_cells, POSTERIOR_PROG,
                    POSTERIOR_SCALE_FACTOR, read_posterior,
                    _save_observations_window, VITERBI_PROG)

MSG_SUCCESS = "____ PROGRAM ENDED SUCCESSFULLY WITH STATUS 0 AT"


SCORE_MIN = 100
SCORE_MAX = 1000

TEMP_DIRPATH = path(gettempdir())

EXT_OPTIONS = {}
EXT_OPTIONS[EXT_INT] = "-of1"
EXT_OPTIONS[EXT_FLOAT] = "-of2"

def make_track_indexes(text):
    return array(map(int, text.split(",")))

re_seg = re.compile(r"^seg\((\d+)\)=(\d+)$")
def parse_viterbi(lines):
    """
    returns: numpy.ndarray of size (num_frames,), type DTYPE_IDENTIFY
    """
    lines = iter(lines)

    # Segment 0, after Island[...]
    assert lines.next().startswith("Segment ")

    # ========
    assert lines.next().startswith("========")

    # Segment 0, number of frames = 1001, viteri-score = -6998.363710
    line = lines.next()
    assert line.startswith("Segment ")

    num_frames_text = line.split(", ")[1].partition(" = ")
    assert num_frames_text[0] == "number of frames"
    assert num_frames_text[1] == " = "

    num_frames = int(num_frames_text[2])

    # Printing random variables from (P,C,E)=(1,999,0) partitions
    line = lines.next()
    assert line.startswith("Printing random variables from (P,C,E)")

    res = zeros(num_frames, DTYPE_IDENTIFY)

    for line in lines:
        # Ptn-0 P': seg(0)=24,seg(1)=24
        if line.startswith(MSG_SUCCESS):
            return res

        assert line.startswith("Ptn-")

        values = line.rpartition(": ")[2]

        for pair in values.split(","):
            match = re_seg.match(pair)
            if not match:
                continue

            index = int(match.group(1))
            val = int(match.group(2))

            res[index] = val

    # shouldn't get to this point
    raise ValueError("%s did not complete successfully" % VITERBI_PROG.prog)

# num_cols is for use by genomedata_count_condition
# XXX: should move this function somewhere else
def write_bed(outfile, start_pos, labels, coord, resolution, num_labels,
              num_cols=None):
    """
    start_pos is an array
    """
    (chrom, region_start, region_end) = coord

    start_pos = region_start + (start_pos * resolution)

    # correct last position which may not be a full sample beyond
    if __debug__:
        remainder = (region_end - region_start) % resolution
        if remainder == 0:
            remainder = resolution
        assert region_end == start_pos[-1] - resolution + remainder
    start_pos[-1] = region_end

    # XXX: why is score_step never used? should this line be removed?
    score_step = (SCORE_MAX - SCORE_MIN) / (num_labels - 1)

    zipper = zip(start_pos[:-1], start_pos[1:], labels)

    # this is easily concatenated since it has no context
    for seg_start, seg_end, seg_label in zipper:
        name = str(seg_label)

        chrom_start = str(seg_start)
        chrom_end = str(seg_end)
        item_rgb = get_label_color(seg_label)

        row = [chrom, chrom_start, chrom_end, name, BED_SCORE, BED_STRAND,
               chrom_start, chrom_end, item_rgb][:num_cols]

        print >>outfile, "\t".join(row)

    # assert that the whole region is mapped
    # seg_end here means the last seg_end in the loop
    assert seg_end == region_end

def save_bed(outfilename, *args, **kwargs):
    with open(outfilename, "w") as outfile:
        write_bed(outfile, *args, **kwargs)

def read_posterior_save_bed(coord, resolution,  outfilename_tmpl, num_labels,
                            infile):
    (chrom, start, end) = coord
    num_frames = ceildiv(end - start, resolution)
    probs = read_posterior(infile, num_frames, num_labels)
    probs_rounded = empty(probs.shape, int)

    outfilenames = []
    for label_index in xrange(num_labels):
        outfilenames.append(outfilename_tmpl % label_index)

    # scale, round, and cast to int
    (probs * POSTERIOR_SCALE_FACTOR).round(out = probs_rounded)

    # print array columns as text to each outfile
    zipper = zip(outfilenames, probs_rounded.T, xrange(num_labels))
    for outfilename, probs_rounded_label, label_index in zipper:
        # run-length encoding on the probs_rounded_label

        outfile = open(outfilename, "w")
        pos, = where(diff(probs_rounded_label) != 0)
        pos = r_[start, pos[:]+start+1, end]

        for bed_start, bed_end in zip(pos[:-1], pos[1:]):
            chrom_start = str(bed_start)
            chrom_end = str(bed_end)
            value = str(probs_rounded_label[bed_start-start])

            row = [chrom, chrom_start, chrom_end, value]
            print >>outfile, "\t".join(row)

def load_posterior_save_bed(coord, resolution, outfilename, num_labels,
                            infilename):
    with open(infilename) as infile:
        read_posterior_save_bed(coord, resolution, outfilename,
                                int(num_labels), infile)

def parse_viterbi_save_bed(coord, resolution, viterbi_lines, bed_filename, num_labels):
    data = parse_viterbi(viterbi_lines)

    start_pos, labels = find_segment_starts(data)

    save_bed(bed_filename, start_pos, labels, coord, resolution, int(num_labels))

def load_viterbi_save_bed(coord, resolution, outfilename, num_labels, infilename):
    with open(infilename) as infile:
        lines = infile.readlines()

    return parse_viterbi_save_bed(coord, resolution, lines, outfilename,
                                  num_labels)

def replace_args_filelistname(args, temp_filepaths, ext):
    option = EXT_OPTIONS[ext]
    filelistname_index = args.index(option) + 1
    fd, filelistname = mkstemp(EXT_LIST + extsep, extsep + ext)
    filelistpath = path(filelistname)

    # side-effect on args, temp_filepaths
    args[filelistname_index] = filelistname
    temp_filepaths.append(filelistpath)

    return fd

def print_to_fd(fd, line):
    with fdopen(fd, "w") as outfile:
        print >>outfile, line

def run_posterior_save_bed(coord, resolution, outfilename, num_labels,
                           genomedata_dirname, float_filename, int_filename,
                           distribution, track_indexes_text, *args):
    # XXX: this whole function is duplicative of run_viterbi_save_bed and needs to be reduced
    # convert from tuple
    args = list(args)

    # a 2,000,000-frame output file is only 84 MiB so it is okay to
    # read the whole thing into memory
    (chrom, start, end) = coord
    track_indexes = make_track_indexes(track_indexes_text)

    float_filepath = TEMP_DIRPATH / float_filename
    int_filepath = TEMP_DIRPATH / int_filename
    temp_filepaths = [float_filepath, int_filepath]

    # XXX: should do something to ensure of1 matches with int, of2 with float
    int_filelistfd = replace_args_filelistname(args, temp_filepaths, EXT_INT)
    float_filelistfd = replace_args_filelistname(args, temp_filepaths,
                                                 EXT_FLOAT)
    with Genome(genomedata_dirname) as genome:
        supercontigs = genome[chrom].supercontigs[start:end]
        assert len(supercontigs) == 1
        supercontig = supercontigs[0]

        continuous_cells = _make_continuous_cells(supercontig, start, end,
                                                  track_indexes)
    try:
        print_to_fd(float_filelistfd, float_filename)
        print_to_fd(int_filelistfd, int_filename)

        _save_observations_window(float_filename, int_filename,
                                  continuous_cells, resolution, distribution)

        # XXXopt: does this actually free the memory? or do we need to
        # use a subprocess to do the loading?

        # remove from memory
        del continuous_cells

        output = POSTERIOR_PROG.getoutput(*args)
    finally:
        for filepath in temp_filepaths:
            # don't raise a nested exception if the file was never created
            try:
                filepath.remove()
            except OSError, err:
                if err.errno == ENOENT:
                    pass

    lines = output.splitlines()
    return read_posterior_save_bed(coord, resolution,  outfilename,
                                   int(num_labels), lines)


def run_viterbi_save_bed(coord, resolution, outfilename, num_labels,
                         genomedata_dirname, float_filename, int_filename,
                         distribution, track_indexes_text, *args):
    # convert from tuple
    args = list(args)
    # a 2,000,000-frame output file is only 84 MiB so it is okay to
    # read the whole thing into memory

    (chrom, start, end) = coord

    track_indexes = make_track_indexes(track_indexes_text)

    float_filepath = TEMP_DIRPATH / float_filename
    int_filepath = TEMP_DIRPATH / int_filename

    temp_filepaths = [float_filepath, int_filepath]

    # XXX: should do something to ensure of1 matches with int, of2 with float
    int_filelistfd = replace_args_filelistname(args, temp_filepaths, EXT_INT)
    float_filelistfd = replace_args_filelistname(args, temp_filepaths,
                                                   EXT_FLOAT)

    with Genome(genomedata_dirname) as genome:
        supercontigs = genome[chrom].supercontigs[start:end]
        assert len(supercontigs) == 1
        supercontig = supercontigs[0]

        continuous_cells = _make_continuous_cells(supercontig, start, end,
                                                  track_indexes)

    try:
        print_to_fd(float_filelistfd, float_filename)
        print_to_fd(int_filelistfd, int_filename)

        _save_observations_window(float_filename, int_filename,
                                  continuous_cells, resolution, distribution)

        # XXXopt: does this work? or do we need to use a subprocess to
        # do the loading?
        # remove from memory
        del continuous_cells

        output = VITERBI_PROG.getoutput(*args)
    finally:
        for filepath in temp_filepaths:
            # don't raise a nested exception if the file was never created
            try:
                filepath.remove()
            except OSError, err:
                if err.errno == ENOENT:
                    pass

    lines = output.splitlines()

    return parse_viterbi_save_bed(coord, resolution, lines, outfilename,
                                  num_labels)

TASKS = {("run", "viterbi"): run_viterbi_save_bed,
         ("load", "viterbi"): load_viterbi_save_bed,
         ("run", "posterior"): run_posterior_save_bed,
         ("load", "posterior"): load_posterior_save_bed}

def task(verb, kind, outfilename, chrom, start, end, resolution, *args):
    start = int(start)
    end = int(end)
    resolution = int(resolution)

    TASKS[verb, kind]((chrom, start, end), resolution, outfilename, *args)

def main(args=sys.argv[1:]):
    if len(args) < 7:
        print >>sys.stderr, \
            "args: VERB KIND OUTFILE CHROM START END RESOLUTION [ARGS...]"
        sys.exit(2)

    return task(*args)

if __name__ == "__main__":
    sys.exit(main())
