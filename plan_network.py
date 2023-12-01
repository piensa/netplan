import requests
from datetime import datetime

datetime_str = datetime.now().strftime("%Y%m%d%H%M%S")

job_name_with_datetime = f"Job_{datetime_str}"

url = 'http://mr.qsel.columbia.edu/jobs'
headers = {
    'Accept': '*/*',
    'Accept-Language': 'en-US,en;q=0.9',
    'Connection': 'keep-alive',
    'Content-Type': 'multipart/form-data',
    'DNT': '1',
    'Origin': 'http://mr.qsel.columbia.edu',
    'Referer': 'http://mr.qsel.columbia.edu/jobs/submit',
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/116.0.0.0 Safari/537.36',
    'X-Requested-With': 'XMLHttpRequest'
}
files = {
    'job_name': (None, job_name_with_datetime),
    'model': (None, 'electrificationplanner'),
    'input_type': (None, 'file_input'),
    'zip_file': ('sample_np.zip', open('sample_np.zip', 'rb'), 'application/zip'),
    'zip_url': (None, '')
}

response = requests.post(url, headers=headers, files=files, verify=False)

print(response.text)



