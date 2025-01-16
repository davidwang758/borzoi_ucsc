from make_trackhub import *
from track_pred import *
import sys
import requests
import re
import uuid

def get_hgsid(chrom,start,end, session_id):
    r = requests.post('https://genome.ucsc.edu/cgi-bin/hgHubConnect', files={"hubUrl":"https://raw.githubusercontent.com/davidwang758/trackhub/refs/heads/main/" + session_id + "/hub.txt"})
    content = r.text
    hgsid = re.findall(r'hgsid=[0-9]*_[a-zA-Z0-9]*', content)[0].split('=')[1]
    link = f'https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid={hgsid}&position={chrom}%3A{start}-{end}'
    print(link)

def generate_session_id():
    session_id = str(uuid.uuid4()).replace("-", "")
    return session_id

if __name__ == "__main__":
    session_id = generate_session_id()
    seq_len = 524288
    center_pos = int(sys.argv[2])
    chrom = sys.argv[1]
    poses = [int(sys.argv[2])]
    alts = [sys.argv[3]]
    predict(seq_len, center_pos, chrom, poses, alts,
        fasta_file = '/home/davidwang/hackweek2025/hg38/assembly/ucsc/hg38.fa',
        params_file = '/home/davidwang/borzoi/examples/params_pred.json',
        targets_file = '/home/davidwang/borzoi/examples/targets_gtex.txt',
        model_file_root= '/home/davidwang/borzoi/examples/saved_models/f3c',
        n_reps=1, rc = True)

    generate_trackhub_dir("/home/davidwang/hackweek2025/hubDirectory/", session_id)
    generate_bigwig(seq_len, chrom, center_pos, session_id)
    convert_wig_to_bigwig(session_id)
    ref_tsv = pd.read_csv("/home/davidwang/hackweek2025/data/y_wt.tsv",sep="\t")
    bigwig_url_root = "https://raw.githubusercontent.com/davidwang758/trackhub/refs/heads/main/" + session_id + "/data/"
    generate_trackDb("/home/davidwang/hackweek2025/hubDirectory/" + session_id + "/hg38/trackDb.txt",ref_tsv.columns, bigwig_url_root)
    upload_trackhub()

    #generate_bigwig(seq_len, chrom, center_pos, session_id)
    #convert_wig_to_bigwig()
    #ref_tsv = pd.read_csv("/home/davidwang/hackweek2025/data/y_wt.tsv",sep="\t")
    #bigwig_url_root = "https://raw.githubusercontent.com/davidwang758/trackhub/refs/heads/main/data/" + session_id + "/"
    #generate_trackDb("/home/davidwang/hackweek2025/hubDirectory/hg38/trackDb.txt",ref_tsv.columns, bigwig_url_root)
    #upload_trackhub()

    start_pos = center_pos - seq_len // 2 + 16 * 32
    end_pos = center_pos + seq_len // 2 - 16 * 32
    get_hgsid(chrom, start_pos, end_pos, session_id)