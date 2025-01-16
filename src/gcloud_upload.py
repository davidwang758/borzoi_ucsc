import glob
import os 
from google.cloud import storage

def create_gcs_client():
    # This shouldn't need any parameters if you set up your authentication as described in the README
    return storage.Client()

def get_bucket(gcs_client, bucket_name):
    return gcs_client.bucket(bucket_name)

def upload_local_directory_to_gcs(local_path, bucket, gcs_path):
    assert os.path.isdir(local_path)
    urls = []
    for local_file in glob.glob(local_path + '/**'):
        if not os.path.isfile(local_file):
           upload_local_directory_to_gcs(local_file, bucket, gcs_path + "/" + os.path.basename(local_file))
        else:
           remote_path = os.path.join(gcs_path, local_file[1 + len(local_path):])
           blob = bucket.blob(remote_path)
           blob.upload_from_filename(local_file)
           urls.append(blob.public_url)

    #return only hub.txt path
    for url in urls:
        basename = os.path.basename(url)
        if basename == "hub.txt":
            return url
        

if __name__ == "__main__":
    gcs_client = create_gcs_client()
    bucket = get_bucket(gcs_client, "seqnn-share")
    upload_dir = "ucsc/borzoi_app/d560af1294a74cb7a8a08b509bfcbffe"
    local_dir = "/home/davidwang/hackweek2025/hubDirectory/d560af1294a74cb7a8a08b509bfcbffe"
    hubtxt_url = upload_local_directory_to_gcs(local_dir, bucket, upload_dir)
    print(hubtxt_url)


