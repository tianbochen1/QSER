from pathlib import Path
from hashlib import sha256
import json
import mne
import numpy as np
import requests

root = Path(__file__).resolve().parents[1]
data = root / 'data'
data.mkdir(parents=True,exist_ok=True)
records = [('open','S001R01','4743b736131a7e147c150e8b37711029b6cda5e356c4b3e8261a03cdcaaf8b0c'),('closed','S001R02','31a95e0a880e6c3d89960d9d62c144f24cc4e9f5d7e93c7f864ef61cd49e847e')]
metadata = []
for condition,record,expected in records:
    path = data/(record+'.edf')
    if not path.exists():
        response = requests.get(f'https://physionet.org/files/eegmmidb/1.0.0/S001/{record}.edf',timeout=90)
        response.raise_for_status()
        path.write_bytes(response.content)
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != expected:
        raise ValueError(f'Unexpected checksum for {record}')
    raw = mne.io.read_raw_edf(path,preload=True,verbose=False)
    channel = next(ch for ch in raw.ch_names if ch.upper().rstrip('.')=='O1')
    if raw.info['sfreq'] != 160:
        raise ValueError('The sampling frequency must be 160 Hz')
    y = raw.get_data(picks=[channel])[0]*1e6
    np.savetxt(data/(condition+'_o1.csv'),y,delimiter=',',header='microvolts',comments='')
    metadata.append(dict(condition=condition,record=record,channel=channel,sampling_frequency_hz=160,samples=len(y),sha256=digest))
(data/'eeg_metadata.json').write_text(json.dumps(metadata,indent=2))
