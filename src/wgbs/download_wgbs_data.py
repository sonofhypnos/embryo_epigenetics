import GEOparse
import requests
import os
from ftplib import FTP
from urllib.parse import urlparse


# Define the GEO Series ID
gse_id = "GSE100272"

# Fetch GEO series
gse = GEOparse.get_GEO(geo=gse_id, destdir="geo_files")

# Function to download files from FTP URLs
def download_ftp_file(ftp_url, output_path):
    parsed_url = urlparse(ftp_url)
    ftp = FTP(parsed_url.hostname)
    ftp.login()  # assumes anonymous login
    ftp.cwd(os.path.dirname(parsed_url.path))

    with open(output_path, "wb") as file:
        ftp.retrbinary(f"RETR {os.path.basename(parsed_url.path)}", file.write)
    ftp.quit()


# Iterate over each sample in the GEO series
for gsm_name, gsm in gse.gsms.items():
    # Get the processed data file URLs
    supplementary_files = gsm.metadata.get("supplementary_file_1", [])
    print(gsm.__dict__)
    print(f"supplementary_files: {supplementary_files}")
    # Download each supplementary file
    for file_url in supplementary_files:
        # TODO: make sure we only get the cg files and not the gc files
        print(file_url)
        file_name = os.path.basename(file_url)
        file_path = os.path.join("data", file_name)
        if not os.path.exists(file_path):
            download_ftp_file(file_url, file_path)

        # # Download the file if it doesn't already exist
        #     print(f"Downloading {file_name} from {file_url}")
        #     response = requests.get(file_url)
        #     with open(file_path, "wb") as file:
        #         file.write(response.content)

print("All processed files downloaded.")
