import json
import h5py
import numpy as np
import pandas as pd
import pysam
from baskerville import seqnn
from baskerville import dna

def read_fasta(fasta_file):
    return pysam.Fastafile(fasta_file)

def read_params(params_file):
    with open(params_file) as params_open :
        params = json.load(params_open)
        params_model = params['model']
        params_train = params['train']
        return params_model, params_train

def read_targets(targets_file):
    targets_df = pd.read_csv(targets_file, index_col=0, sep='\t')
    target_index = targets_df.index 
    return targets_df, target_index

def read_models(model_file_root, params_model, target_index, slice_pair, n_reps, rc):
    models = []
    for rep_ix in range(n_reps) :
        model_file = model_file_root + f"{rep_ix}/train/model0_best.h5"
        seqnn_model = seqnn.SeqNN(params_model)
        seqnn_model.restore(model_file, 0)
        seqnn_model.build_slice(target_index)
        if rc :
            seqnn_model.strand_pair.append(slice_pair)
        seqnn_model.build_ensemble(rc, [0])
        
        models.append(seqnn_model)
    return models

def mut_sequence(sequence_one_hot_mut, poses, alts, start):
    for pos, alt in zip(poses, alts) :
        alt_ix = -1
        if alt == 'A' :
            alt_ix = 0
        elif alt == 'C' :
            alt_ix = 1
        elif alt == 'G' :
            alt_ix = 2
        elif alt == 'T' :
            alt_ix = 3

    sequence_one_hot_mut[pos-start-1] = 0.
    sequence_one_hot_mut[pos-start-1, alt_ix] = 1.
    return sequence_one_hot_mut

def make_seq_1hot(genome_open, chrm, start, end, seq_len):
    if start < 0:
        seq_dna = "N" * (-start) + genome_open.fetch(chrm, 0, end)
    else:
        seq_dna = genome_open.fetch(chrm, start, end)

    # Extend to full length
    if len(seq_dna) < seq_len:
        seq_dna += "N" * (seq_len - len(seq_dna))

    seq_1hot = dna.dna_1hot(seq_dna)
    return seq_1hot

def process_sequence(fasta_open, chrom, start, end, seq_len=524288):

    seq_len_actual = end - start

    # Pad sequence to input window size
    start -= (seq_len - seq_len_actual) // 2
    end += (seq_len - seq_len_actual) // 2

    # Get one-hot
    sequence_one_hot = make_seq_1hot(fasta_open, chrom, start, end, seq_len)

    return sequence_one_hot.astype("float32")

def predict_tracks(models, sequence_one_hot):

    predicted_tracks = []
    
    #Loop over model replicates
    for rep_ix in range(len(models)):

        #Predict coverage and store as float16
        yh = models[rep_ix](sequence_one_hot[None, ...])[:, None, ...].astype(
            "float16"
        )

        predicted_tracks.append(yh)

    #Concatenate across replicates
    predicted_tracks = np.concatenate(predicted_tracks, axis=1)

    return predicted_tracks

def transforms(y_wt, y_mut, tissues, track_indicies):

    track_indicies = [np.arange(0, 89).tolist()] + track_indicies
    tissues = ["GTEx_all_tissues"] + tissues
    
    track_scale = 0.01 
    track_transform = 3./4. # keep
    clip_soft = 384. #outlier removal

    preds_wt = {}
    preds_mut = {} 
    
    y_wt_curr = np.array(np.copy(y_wt), dtype=np.float32)
    y_mut_curr = np.array(np.copy(y_mut), dtype=np.float32)

    y_wt_curr /= track_scale
    y_mut_curr /= track_scale

    y_wt_curr_unclipped = (y_wt_curr - clip_soft) ** 2 + clip_soft
    y_mut_curr_unclipped = (y_mut_curr - clip_soft) ** 2 + clip_soft

    unclip_mask_wt = y_wt_curr > clip_soft
    unclip_mask_mut = y_mut_curr > clip_soft

    y_wt_curr[unclip_mask_wt] = y_wt_curr_unclipped[unclip_mask_wt]
    y_mut_curr[unclip_mask_mut] = y_mut_curr_unclipped[unclip_mask_mut]

    #Undo sqrt
    y_wt_curr = y_wt_curr ** (1. / track_transform)
    y_mut_curr = y_mut_curr ** (1. / track_transform)

    for tissue, track_index in zip(tissues, track_indicies):
        preds_wt[tissue] = np.mean(y_wt_curr[..., track_index], axis=(0, 1, 3))
        preds_mut[tissue] = np.mean(y_mut_curr[..., track_index], axis=(0, 1, 3))

    return pd.DataFrame(preds_wt), pd.DataFrame(preds_mut)

def predict(seq_len, center_pos, chrom, poses, alts,
        fasta_file, params_file, targets_file, model_file_root,
        n_reps=1, rc=True):

    start = center_pos - seq_len // 2
    end = center_pos + seq_len // 2

    params_model, params_train = read_params(params_file)
    targets_df, target_index = read_targets(targets_file)

    if rc :
        strand_pair = targets_df.strand_pair
        
        target_slice_dict = {ix : i for i, ix in enumerate(target_index.values.tolist())}
        slice_pair = np.array([
            target_slice_dict[ix] if ix in target_slice_dict else ix for ix in strand_pair.values.tolist()
        ], dtype='int32')

    models = read_models(model_file_root, params_model, target_index, slice_pair, n_reps, rc) 

    fasta_open = read_fasta(fasta_file)
    sequence_one_hot_wt = process_sequence(fasta_open, chrom, start, end)
    sequence_one_hot_mut = mut_sequence(np.copy(sequence_one_hot_wt), poses, alts, start)

    y_wt = predict_tracks(models, sequence_one_hot_wt)
    y_mut = predict_tracks(models, sequence_one_hot_mut)

    idx = []
    tissues = []
    for t in np.unique(targets_df["description"]):
        tissues.append(t)
        idx.append(np.where(targets_df["description"].values == t)[0])
    y_wt, y_mut = transforms(y_wt, y_mut, tissues, idx)

    #print(y_wt.shape)
    #print(y_wt)
    #print(y_mut.shape)
    #print(y_mut)

    y_wt.to_csv("/home/davidwang/hackweek2025/data/y_wt.tsv",index=None,sep="\t")
    y_mut.to_csv("/home/davidwang/hackweek2025/data/y_mut.tsv",index=None,sep="\t")

    #pd.DataFrame(y_wt).to_csv("/home/davidwang/hackweek2025/data/y_wt.txt",index=None,header=None)
    #pd.DataFrame(y_mut).to_csv("/home/davidwang/hackweek2025/data/y_mut.txt",index=None,header=None)


if __name__ == "__main__":
    seq_len = 524288

    center_pos = 135548708
    chrom = 'chr9'
    poses = [135548708]
    alts = ['C']
    predict(seq_len, center_pos, chrom, poses, alts,
        fasta_file = '/home/davidwang/hackweek2025/hg38/assembly/ucsc/hg38.fa',
        params_file = '/home/davidwang/borzoi/examples/params_pred.json',
        targets_file = '/home/davidwang/borzoi/examples/targets_gtex.txt',
        model_file_root= '/home/davidwang/borzoi/examples/saved_models/f3c',
        n_reps=1, rc = True)