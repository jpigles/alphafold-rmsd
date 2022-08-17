# Imports the Google Cloud Client library
from google.cloud import storage

#This downloads the item from the bucket. I am testing it on one file.
def download_blob(bucket_name, source_blob_name, destination_file_name):
    '''Downloads a blob from the bucket.'''
    #The ID of your GCS bucket
    bucket_name = 'public-datasets-deepmind-alphafold'

    # The ID of your GCS object
    source_blob_name = 'AF-Q9BXL7-F1-model_v3.cif'

    # The path to which the file should be downloaded.
    destination_file_name = "~/Documents/gsponer_lab/auotinhibition_protein_data/alphafold_files/AF-Q9BXL7-F1-model_v3.cif"

    storage_client = storage.Client()

    bucket = storage_client.bucket(bucket_name)

    # Construct a client side representation of a blob.
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)

    print(
        'Downloaded storage object {} from bucket {} to local file {}.'.format(source_blob_name, bucket_name, destination_file_name
        )
    )