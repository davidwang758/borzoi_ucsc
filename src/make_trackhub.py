import pandas as pd
import pyBigWig
import subprocess
import uuid
import os
import sys
from gcloud_upload import *

def generate_trackhub_dir(root_dir, session_id):
    os.mkdir(root_dir + session_id)
    os.mkdir(root_dir + session_id + "/hg38")
    os.mkdir(root_dir + session_id + "/data")
    with open(root_dir + session_id + "/hub.txt", "w") as f:
        f.write("hub Borzoi Prediction Tracks " + session_id + "\n")
        f.write("shortLabel Borzoi Prediction Tracks\n")
        f.write("longLabel Borzoi prediction tracks for reference and alternative allele.\n")
        f.write("genomesFile genomes.txt\n")
        f.write("email  myEmail@address\n")
    with open(root_dir + session_id + "/genomes.txt", "w") as f:
        f.write("genome hg38\n")
        f.write("trackDb hg38/trackDb.txt")

def generate_trackDb(trackDb_file, tissues, bigwig_url_root):
    with open (trackDb_file,"w") as f:
        for t in tissues:
            f.write("track " + t + "\n")
            f.write("type bigWig\n")
            f.write("container multiWig\n")
            f.write("shortLabel " + t + "\n")
            f.write("longLabel Model predictions averaged over: " + t + "\n")
            f.write("visibility full\n")
            f.write("aggregate solidOverlay\n")
            f.write("showSubtrackColorOnUi on\n")
            f.write("maxHeightPixels 100:50:20\n")
            #f.write("viewLimits 1:100\n")
            f.write("autoScale on\n")
            f.write("priority 20\n")
            f.write("\n")
            f.write("    track  REF_" + t + "\n")
            f.write("    bigDataUrl " + bigwig_url_root + t + "_y_wt.bw\n")
            f.write("    shortLabel REF_" + t + "\n")
            f.write("    longLabel Reference allele predictions.\n")
            f.write("    graphTypeDefault bar\n")
            f.write("    type bigWig\n")
            f.write("    parent " + t + "\n")
            f.write("    color 255,0,0\n")
            f.write("\n")
            f.write("    track  ALT_" + t + "\n")
            f.write("    bigDataUrl " + bigwig_url_root + t + "_y_mut.bw\n")
            f.write("    shortLabel ALT_" + t + "\n")
            f.write("    longLabel Alternative allele predictions.\n")
            f.write("    graphTypeDefault bar\n")
            f.write("    type bigWig\n")
            f.write("    parent " + t + "\n")
            f.write("    color 0,0,255\n")
            f.write("\n")

def generate_bigwig(seq_len, chrom, center_pos, session_id, data_dir):
    start_pos = center_pos - seq_len // 2 + 16 * 32

    ref_tsv = pd.read_csv(data_dir + "y_wt.tsv",sep="\t")
    alt_tsv = pd.read_csv(data_dir + "y_mut.tsv",sep="\t")
    header_ref=f'track type=wiggle_0 name="REF" description="Reference Allele Track" visibility=full autoScale=off viewLimits=0:1000 color=0,200,100 maxHeightPixels=100:50:20 graphType=bar priority=20\nfixedStep chrom={chrom} start={start_pos} step=32 span=32'
    header_alt=f'track type=wiggle_0 name="ALT" description="Alternative Allele Track" visibility=full autoScale=off viewLimits=0:1000 color=0,100,100 maxHeightPixels=100:50:20 graphType=bar priority=20\nfixedStep chrom={chrom} start={start_pos} step=32 span=32'
    
    for t in ref_tsv.columns:
        ref_data = ref_tsv[t].values
        alt_data = alt_tsv[t].values

        with open(data_dir + session_id + f"/data/{t}_y_wt.wig","w") as f:
            f.write(header_ref + "\n")
            for i in ref_data:
                f.write(str(i) + "\n")

        with open(data_dir + session_id + f"/data/{t}_y_mut.wig","w") as f:
            f.write(header_alt + "\n")
            for i in alt_data:
                f.write(str(i) + "\n")

def convert_wig_to_bigwig(session_id, data_dir, utils_dir):
    subprocess.call(["../utils/convert_to_bigwig.sh", session_id, data_dir, utils_dir])

def upload_trackhub(bucket_name, upload_dir, local_dir):
    gcs_client = create_gcs_client()
    bucket = get_bucket(gcs_client, bucket_name)
    hubtxt_url = upload_local_directory_to_gcs(local_dir, bucket, upload_dir)
    return hubtxt_url

if __name__ == "__main__":
    session_id = str(uuid.uuid4()).replace("-", "")
    print(session_id)
    seq_len = 524288
    center_pos = 135548708
    chrom = "chr9"
    upload_dir = "ucsc/borzoi_app/" + session_id
    local_dir = "/home/davidwang/hackweek2025/hubDirectory/" + session_id
    bucket_name = "seqnn-share"
    generate_trackhub_dir("/home/davidwang/hackweek2025/hubDirectory/", session_id)
    generate_bigwig(seq_len, chrom, center_pos, session_id)
    convert_wig_to_bigwig(session_id)
    ref_tsv = pd.read_csv("/home/davidwang/hackweek2025/data/y_wt.tsv",sep="\t")
    bigwig_url_root = "https://storage.googleapis.com/seqnn-share/ucsc/borzoi_app/" + session_id + "/data/"
    generate_trackDb("/home/davidwang/hackweek2025/hubDirectory/" + session_id + "/hg38/trackDb.txt",ref_tsv.columns, bigwig_url_root)
    upload_trackhub(bucket_name, upload_dir, local_dir)
