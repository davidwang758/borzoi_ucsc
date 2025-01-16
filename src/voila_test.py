import requests
import re
from io import StringIO

with open("/home/davidwang/hackweek2025/data/test_annotation.bed") as file:
	genome = "hg38"
	_chr = "chr19"
	start = "58352928"
	end =  "58353327"
	file_content = file.read()
	tmp_file = StringIO(file_content)

	files = {'hgct_customText': tmp_file}
	payload = {'db': genome}

	r = requests.post('https://genome.ucsc.edu/cgi-bin/hgCustom', files=files, data=payload)

	content = r.text

	hgsid = re.findall(r'hgsid=[0-9]*_[a-zA-Z0-9]*', content)[0].split('=')[1]
	link = f'https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid={hgsid}&position={_chr}%3A{start}-{end}'
	print(link)
