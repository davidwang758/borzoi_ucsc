from make_trackhub import *
#from track_pred import *
import sys
import requests
import re
import uuid

def get_link(chrom,start,end, session_id, cloud_url_root):
    r = requests.post('https://genome.ucsc.edu/cgi-bin/hgHubConnect', files={"hubUrl":cloud_url_root + session_id + "/hub.txt"})
    content = r.text
    hgsid = re.findall(r'hgsid=[0-9]*_[a-zA-Z0-9]*', content)[0].split('=')[1]
    link = f'https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid={hgsid}&position={chrom}%3A{start}-{end}'
    return link

def generate_session_id():
    session_id = str(uuid.uuid4()).replace("-", "")
    return session_id

def create_ucsc_link(install_dir, chrom, center_pos, alts, offset, borzoi_session_id="", window_size):
    session_id = generate_session_id()
    if borzoi_session_id != "":
        borzoi_session_id = borzoi_session_id + "_"
    seq_len = 524288
    poses = [center_pos]
    alts = [alts]

    data_dir = install_dir + "/data/"
    util_dir = install_dir + "/utils/"
    upload_dir = "ucsc/borzoi_app/" + session_id
    local_dir = data_dir + session_id
    bucket_name = "seqnn-share"
    generate_trackhub_dir(data_dir, session_id)
    generate_bigwig(seq_len, chrom, center_pos, session_id, data_dir, borzoi_session_id)
    convert_wig_to_bigwig(session_id, data_dir, util_dir)
    ref_tsv = pd.read_csv(data_dir + borzoi_session_id + "y_wt.tsv",sep="\t")
    cloud_url_root = "https://storage.googleapis.com/seqnn-share/ucsc/borzoi_app/" 
    bigwig_url_root = cloud_url_root + session_id + "/data/"
    generate_trackDb(data_dir + session_id + "/hg38/trackDb.txt",ref_tsv.columns, bigwig_url_root)
    upload_trackhub(bucket_name, upload_dir, local_dir)
    delete_temp_dir(local_dir)

    # Remove offset if needed
    window_size = window_size // 32 * 32
    start_pos = center_pos - window_size // 2 + offset * 32
    end_pos = center_pos + window_size // 2 - offset * 32
    out_link = get_link(chrom, start_pos, end_pos, session_id, cloud_url_root)
    return out_link

if __name__ == "__main__":
    chrom = sys.argv[1]
    center_pos = int(sys.argv[2])
    alts = sys.argv[3]
    offset = int(sys.argv[4])
    install_dir = sys.argv[5] #"/home/davidwang/hackweek2025"
    borzoi_session_id = sys.argv[6]
    window_size = int(sys.argv[7])
    out_link = create_ucsc_link(install_dir, chrom, center_pos, alts, offset, borzoi_session_id, window_size)
    print(out_link)
    """
    session_id = generate_session_id()
    seq_len = 524288
    center_pos = int(sys.argv[2])
    chrom = sys.argv[1]
    poses = [int(sys.argv[2])]
    alts = [sys.argv[3]]

    install_dir = "/home/davidwang/hackweek2025"
    data_dir = install_dir + "/data/"
    util_dir = install_dir + "/utils/"
    upload_dir = "ucsc/borzoi_app/" + session_id
    local_dir = data_dir + session_id
    bucket_name = "seqnn-share"
    generate_trackhub_dir(data_dir, session_id)
    generate_bigwig(seq_len, chrom, center_pos, session_id, data_dir)
    convert_wig_to_bigwig(session_id, data_dir, util_dir)
    ref_tsv = pd.read_csv(data_dir + "y_wt.tsv",sep="\t")
    cloud_url_root = "https://storage.googleapis.com/seqnn-share/ucsc/borzoi_app/" 
    bigwig_url_root = cloud_url_root + session_id + "/data/"
    generate_trackDb(data_dir + session_id + "/hg38/trackDb.txt",ref_tsv.columns, bigwig_url_root)
    upload_trackhub(bucket_name, upload_dir, local_dir)

    # Remove offset if needed
    start_pos = center_pos - seq_len // 2 + 16 * 32
    end_pos = center_pos + seq_len // 2 - 16 * 32
    out_link = get_link(chrom, start_pos, end_pos, session_id, cloud_url_root)
    print(out_link) 
    """