#!/usr/bin/env python3

"""
License

Copyright (c) 2014, Ole Lund, Technical University of Denmark All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
"""


from __future__ import division
from argparse import ArgumentParser
from tabulate import tabulate
from distutils.spawn import find_executable
import sys
import os
import time
import subprocess
import json
import gzip
import pprint
import pandas as pd


##########################################################################
# FUNCTIONS
##########################################################################


def text_table(headers, rows, empty_replace='-'):
    ''' Create text table

    USAGE:
        >>> from tabulate import tabulate
        >>> headers = ['A','B']
        >>> rows = [[1,2],[3,4]]
        >>> print(text_table(headers, rows))
        **********
          A     B
        **********
          1     2
          3     4
        ==========
    '''
    # Replace empty cells with placeholder
    rows = map(lambda row: map(lambda x: x if x else empty_replace, row), rows)
    # Create table
    table = tabulate(rows, headers, tablefmt='simple').split('\n')
    # Prepare title injection
    width = len(table[0])
    # Switch horisontal line
    table[1] = '*' * (width + 2)
    # Update table with title
    table = (("%s\n" * 3)
             % ('*' * (width + 2), '\n'.join(table), '=' * (width + 2)))
    return table


def is_gzipped(file_path):
    ''' Returns True if file is gzipped and False otherwise.
         The result is inferred from the first two bits in the file read
         from the input path.
         On unix systems this should be: 1f 8b
         Theoretically there could be exceptions to this test but it is
         unlikely and impossible if the input files are otherwise expected
         to be encoded in utf-8.
    '''
    with open(file_path, mode='rb') as fh:
        bit_start = fh.read(2)
    if(bit_start == b'\x1f\x8b'):
        return True
    else:
        return False


def get_file_format(input_files):
    """
    Takes all input files and checks their first character to assess
    the file format. Returns one of the following strings; fasta, fastq,
    other or mixed. fasta and fastq indicates that all input files are
    of the same format, either fasta or fastq. other indiates that all
    files are not fasta nor fastq files. mixed indicates that the inputfiles
    are a mix of different file formats.
    """

    # Open all input files and get the first character
    file_format = []
    invalid_files = []
    for infile in input_files:
        if is_gzipped(infile):
            f = gzip.open(infile, "rb")
            fst_char = f.read(1)
        else:
            f = open(infile, "rb")
            fst_char = f.read(1)
        f.close()
        # Assess the first character
        if fst_char == b"@":
            file_format.append("fastq")
        elif fst_char == b">":
            file_format.append("fasta")
        else:
            invalid_files.append("other")
    if len(set(file_format)) != 1:
        return "mixed"
    return ",".join(set(file_format))

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################


parser = ArgumentParser()

parser.add_argument("-i", "--infile",
                    help="FASTA or FASTQ input files.",
                    nargs="+",
                    required=True)
parser.add_argument("-o", "--outputPath",
                    dest="outdir",
                    help="Path to blast output",
                    default='.')
parser.add_argument("-tmp", "--tmp_dir",
                    dest="tmpdir",
                    help=("Temporary directory for storage of the results "
                          "from the external software."))
parser.add_argument("-kp", "--kmaPath",
                    dest="method_path",
                    help="Path to kma")
parser.add_argument("-p", "--databasePath",
                    dest="db_path",
                    help="Path to the database",
                    default='/database')
parser.add_argument("-q", "--quiet",
                    action="store_true")


args = parser.parse_args()

##########################################################################
# MAIN
##########################################################################

if args.quiet:
    f = open('/dev/null', 'w')
    sys.stdout = f

# Defining varibales
method_path = args.method_path

# Check if valid database is provided
if args.db_path is None:
    sys.exit("Input Error: No database directory was provided!\n")
elif not os.path.exists(args.db_path):
    sys.exit("Input Error: The specified database directory does not exist!\n")
else:
    # Check existence of config file
    db_config_file = '%s/config' % (args.db_path)
    if not os.path.exists(db_config_file):
        sys.exit("Input Error: The database config file could not be found!")
    db_path = args.db_path

# Check if valid input files are provided
if args.infile is None:
    sys.exit("Input Error: No input file was provided!\n")
elif not os.path.exists(args.infile[0]):
    sys.exit("Input Error: Input file does not exist!\n")
elif len(args.infile) > 1:
    if not os.path.exists(args.infile[1]):
        sys.exit("Input Error: Input file does not exist!\n")
    infile = args.infile
else:
    infile = args.infile

# Check if valid output directory is provided
if not os.path.exists(args.outdir):
    sys.exit("Input Error: Output dirctory does not exist!\n")
outdir = os.path.abspath(args.outdir)

# Check if valid tmp directory is provided
if args.tmpdir:
    if not os.path.exists(args.tmpdir):
        sys.exit("Input Error: Tmp dirctory, {}, does not exist!\n"
                 .format(args.tmpdir))
    else:
        tmpdir = os.path.abspath(args.tmpdir)
else:
    tmpdir = outdir

# Check if databases and config file are correct/correponds
dbs = {}
extensions = []
db_description = {}
with open(db_config_file) as f:
    for i in f:
        i = i.strip()
        if i == '':
            continue
        if i[0] == '#':
            if 'extensions:' in i:
                extensions = [
                    s.strip() for s in i.split('extensions:')[-1].split(',')]
            continue
        tmp = i.split('\t')
        if len(tmp) != 3:
            sys.exit(("Input Error: Invalid line in the database"
                      " config file!\nA proper entry requires 3 tab "
                      "separated columns!\n%s") % (i))
        db_prefix = tmp[0].strip()
        name = tmp[1].split('#')[0].strip()
        description = tmp[2]

        # Check if all db files are present
        for ext in extensions:
            db = "%s/%s.%s" % (db_path, db_prefix, ext)
            if not os.path.exists(db):
                sys.exit(("Input Error: The database file (%s) "
                          "could not be found!") % (db))
        if db_prefix not in dbs:
            dbs[db_prefix] = []
        dbs[db_prefix].append(name)
        db_description[db_prefix] = description

if len(dbs) == 0:
    sys.exit("Input Error: No databases were found in the "
             "database config file!")

database = str(db_path)+str(list(dbs.keys())[0])
file_format = get_file_format(infile)

# Call appropriate method (kma or blastn) based on file format
if file_format == "fastq":

    if not method_path:
        method_path = "kma"
    if find_executable(method_path) is None:
        sys.exit("No valid path to a kma program was provided. Use the -kp "
                 "flag to provide the path.")

    # Check the number of files
    if len(infile) == 1:
        infile_1 = infile[0]
        infile_2 = ""
    elif len(infile) == 2:
        infile_1 = infile[0]
        infile_2 = infile[1]
    else:
        sys.exit("Only 2 input file accepted for raw read data,\
                     if data from more runs is avaliable for the same\
                     sample, please concatinate the reads into two files")

    sample_name = os.path.basename(sorted(args.infile)[0])

    cmd = str(method_path)+" -i "+str(infile_1)+" "+str(infile_2)+" -o " + \
        str(args.tmpdir)+str(sample_name)+"_16srna"+" -t_db "+str(database) + \
        " -1t1 -mem_mode -Sparse"
    print(cmd)

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = process.communicate()

elif file_format == "fasta":

    if not method_path:
        method_path = "kma"
    if find_executable(method_path) is None:
        sys.exit("No valid path to a kma program was provided. Use the "
                 "-mp flag to provide the path.")

    # Assert that only one fasta file is inputted
    assert len(infile) == 1, "Only one input file accepted for assembled data"
    infile = infile[0]
    sample_name = os.path.basename(sorted(args.infile)[0])
    # Call BLASTn
    cmd = str(method_path)+" -i "+str(infile)+" -o "+str(args.tmpdir) + \
        str(sample_name)+"_16srna"+" -t_db "+str(database) + \
        " -mem_mode -Sparse"
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = process.communicate()

else:
    sys.exit("Input file must be fastq or fasta format, not " + file_format)

filename = str(args.tmpdir)+str(sample_name)+"_16srna"".spa"
json_results = dict()

if pd.read_csv(filename, sep="\t").empty:
    json_results['Message'] = 'No Results Found!'
    json_results['Species'] = ''
    json_results["Confidence of result"] = ''
    json_results['Match'] = ''
else:
    alignment_data = pd.read_csv(filename, sep="\t").iloc[0]

    json_results['Template'] = alignment_data['#Template']
    json_results['Species'] = alignment_data['#Template'].split(';')[-1]
    json_results['Match'] = alignment_data['#Template'].split(';')[0].split('.')[0]
    json_results['Database'] = database
    if file_format == "fasta":
        if alignment_data['Template_Coverage'] >= 98.0:
            json_results["Confidence of result"] = "PASS"
        else:
            json_results["Confidence of result"] = "FAIL"
    else:
        if alignment_data['Template_Coverage'] >= 99.5:
            json_results["Confidence of result"] = "PASS"
        else:
            json_results["Confidence of result"] = "FAIL"

# Get run info for JSON file
service = os.path.basename(__file__).replace(".py", "")
date = time.strftime("%d.%m.%Y")
time = time.strftime("%H:%M:%S")

# Make JSON output file
data = {service: {}}

userinput = {"filename(s)": args.infile,
             "method": "kma",
             "file_format": file_format}
run_info = {"date": date, "time": time}

data[service]["user_input"] = userinput
data[service]["run_info"] = run_info
data[service]["results"] = json_results

pprint.pprint(data)

# Save json output
result_file = "{}/data.txt".format(outdir)
with open(result_file, "w") as outfile:
    json.dump(data, outfile)

if args.quiet:
    f.close()
