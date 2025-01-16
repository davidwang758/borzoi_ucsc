# Generate UCSC Genome Browser Links with Borzoi Prediction Tracks

Note: The code may not fully run on another machine at the moment because the temp file hosting is currently on my git account. I will transfer this to gCloud tomorrow. No external python dependencies aside from what is installed during the Borzoi installation.

To run the code:

```
python get_ucsc_link.py <chrom> <variant_position> <alt_allele>
```

Code logic

1. Start at src/get_ucsc_link.py
2. This calls src/track_pred.py to get the model output matrix, transforms the output, aggregates the tracks with the same tissues, and outputs this as a tsv in data/y_wt.tsv and data/y_mut.tsv 
3. Then it calls src/make_trackhub.py which generates the TrackHub directory structure with relevant metadata and the wig/bigwig files. See hubDirectory/<session_id> for an example. Then it pushes this directory to another git repo (this is getting replaced with gCloud).
4. src/get_ucsc_link.py will then submit a post request with the TrackHub hub.txt location as the payload to UCSC Genome Browser. This returns this hgsid session id which is appended to the link.  

