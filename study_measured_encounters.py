"""Run the frozen measured-data screen; no model selection on held-out outcomes."""
import json
from pathlib import Path
from measured_encounters import fetch_workbook, load_rows, screen, ARCHIVE_URL, ARCHIVE_SHA, WORKBOOK_SHA

if __name__ == '__main__':
    path = fetch_workbook()
    result = {'protocol_commit':'8adbf4f', 'source':{'url':ARCHIVE_URL,'sha256':ARCHIVE_SHA,
               'workbook_sha256':WORKBOOK_SHA},'analyses':{}}
    for name, secondary in [('primary_N',False), ('secondary_NC',True)]:
        rows, counts = load_rows(path, secondary)
        result['analyses'][name] = {'exclusions':counts, **screen(rows, secondary)}
        summary = {k:v for k,v in result['analyses'][name].items() if k != 'folds'}
        print(name, json.dumps(summary, indent=2), flush=True)
    Path('results/measured_encounters.json').write_text(json.dumps(result, indent=2)+'\n')
