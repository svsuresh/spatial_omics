from Bio import SeqIO
from gzip import open as gzopen
import sys
import datetime as dt

input_filename = sys.argv[1]
#input_filename = "test.fastq.gz"
output_filename = input_filename.rsplit(
    ".fastq.gz", 2)[0]+".processed.fastq"
print("\nInput file is " + input_filename)
print("Output file name is " + output_filename)
_start_time = dt.datetime.now().replace(microsecond=0)
print("Extraction process started at " + str(_start_time))

output_handle = open(output_filename, "w")

with gzopen(input_filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      cut_record = record[32:40] + record[70:78] + record[22:32]  # BC2 + BC1 + UMI
      SeqIO.write(cut_record, output_handle, "fastq")

output_handle.close()

#cut_record = [record[32:40] + record[70:78] + record[22:32]
#              for record in SeqIO.parse(gzopen(input_filename, "rt"), "fastq")]         # BC2 + BC1 + UMI
#
#print("Writing the file now " + str(dt.datetime.now().replace(microsecond=0)))
#
#with gzopen(output_filename, "wt") as handle:
#    SeqIO.write(cut_record, handle, "fastq")
#
_end_time = dt.datetime.now().replace(microsecond=0)
print("done " + str(_end_time))
print("Total time taken "+str(_end_time-_start_time))
