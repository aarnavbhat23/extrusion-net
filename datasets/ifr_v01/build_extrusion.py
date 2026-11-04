"""Curate measured iFRAP observations; never infer missing biological truth."""
import hashlib
import json
from pathlib import Path
import re
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE_URL = 'https://data.mendeley.com/public-files/datasets/spp4x8jpdy/files/08b724e0-391c-4912-8b90-59cca05ed6c0/file_downloaded'
SOURCE_SHA256 = '5f9e49088e2218d6f35bca61f3039da72864906de711740437949e11cb47ec0b'
SOURCE_COLUMNS = ['mean_list_interp_bleached', 'mean_list_interp_unbleached',
                  'background_int', 'unfrap_cell', 'nucleus', 'frame',
                  'tagged_prot', 'depletion', 'filename']
# Frozen from identifiers and condition coverage, before downstream model fitting.
# Entire days remain together. No random nucleus/frame splits are permitted.
RESERVED_DAYS = {'20240201', '20240205', '20240429', '20240606'}
DEVELOPMENT_DAYS = {'20240116', '20240118', '20240122', '20240124',
                    '20240126', '20240207', '20240503', '20240529', '20240604'}
QUARANTINE_DAYS = {'20240220', '202402202', '20240426', '20240527'}


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def verify_source(path):
    if sha256(path) != SOURCE_SHA256:
        raise ValueError('Published source SHA-256 mismatch; refusing to build')


def acquire():
    path = ROOT / 'source/combined_dataframe.csv'
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(SOURCE_URL, headers={
            'User-Agent': 'Mozilla/5.0',
            'Referer': 'https://data.mendeley.com/datasets/spp4x8jpdy/1'})
        with urllib.request.urlopen(req, timeout=60) as incoming:
            data = incoming.read(10_000_001)
        if len(data) > 10_000_000:
            raise ValueError('Source exceeds bounded acquisition size')
        temporary = path.with_suffix('.partial')
        temporary.write_bytes(data)
        verify_source(temporary)
        temporary.replace(path)
    verify_source(path)
    return path


def filename_treatment(filename):
    """Diagnostic only: filename evidence is NOT an authoritative relabeling."""
    lower = filename.lower()
    if any(x in lower for x in ('noaux', '0h_aux')):
        return '- aux'
    if re.search(r'6h_?aux', lower):
        return '+ aux'
    if 'nodtag' in lower:
        return '- dtag'
    if re.search(r'6h_?dtag', lower):
        return '+ dtag'
    return None


def split_for_day(day):
    if day in RESERVED_DAYS:
        return 'reserved_test'
    if day in DEVELOPMENT_DAYS:
        return 'development'
    if day in QUARANTINE_DAYS:
        return 'quarantine'
    raise ValueError(f'Unreviewed acquisition-day identifier: {day}')


def validate_observations(df):
    if list(df.columns) != SOURCE_COLUMNS or df.isna().any().any():
        raise ValueError('Unexpected source schema or missing values')
    if df.duplicated(['filename', 'nucleus', 'frame']).any():
        raise ValueError('Duplicate observation identity')
    for _, group in df.groupby(['filename', 'nucleus']):
        if not np.array_equal(np.sort(group.frame), np.arange(250)):
            raise ValueError('Incomplete or duplicated frame sequence')
    shared = df.groupby(['filename', 'frame'])[['unfrap_cell', 'background_int']].nunique()
    if not shared.eq(1).all().all():
        raise ValueError('Movie controls differ between nuclei')
    if not df.groupby('filename')[['tagged_prot', 'depletion']].nunique().eq(1).all().all():
        raise ValueError('Inconsistent treatment within movie')


def author_curve(group):
    """Replay source notebook's movie-level normalization, with explicit time units."""
    by_frame = group.groupby('frame', sort=True)[SOURCE_COLUMNS[:4]].mean()
    controls = group.drop_duplicates('frame').sort_values('frame')
    denominator = (controls.unfrap_cell - controls.background_int).to_numpy()
    if np.any(denominator <= 0):
        raise ValueError('Nonpositive background-corrected reference signal')
    correction = denominator / denominator[:4].mean()
    signal = (by_frame.mean_list_interp_unbleached - by_frame.mean_list_interp_bleached).to_numpy() / correction
    if not np.isfinite(signal[5]) or signal[5] == 0:
        raise ValueError('Invalid frame-5 normalization reference')
    return pd.DataFrame({'frame': np.arange(250),
                         'time_since_first_frame_s': 30 * np.arange(250),
                         'time_since_normalization_s': 30 * (np.arange(250)-5),
                         'photobleaching_correction': correction,
                         'corrected_difference': signal,
                         'normalized_ifrap': signal / signal[5],
                         'source_postnormalization_window': np.arange(250) >= 5})


def main():
    source = acquire()
    df = pd.read_csv(source)
    validate_observations(df)
    out = ROOT/'extrusion'
    out.mkdir(exist_ok=True)
    movies, issues = [], []
    for filename, group in df.groupby('filename', sort=True):
        day = filename.split('_')[0]
        date = pd.to_datetime(day, format='%Y%m%d', errors='coerce') if len(day) == 8 else pd.NaT
        evidence = filename_treatment(filename)
        label = group.depletion.iloc[0]
        conflict = evidence is not None and evidence != label
        reasons = []
        if conflict:
            reasons.append('filename_treatment_conflict')
        if pd.isna(date):
            reasons.append('malformed_acquisition_date')
        if day == '20240220':
            reasons.append('possible_day_alias_of_malformed_202402202')
        split = split_for_day(day)
        if conflict and split != 'quarantine':
            raise ValueError('A conflicting label entered an evaluation split')
        movie_id = 'ifr_' + hashlib.sha256(filename.encode()).hexdigest()[:16]
        movie = dict(movie_id=movie_id, filename=filename, date_raw=day,
                     acquisition_date=None if pd.isna(date) else date.strftime('%Y-%m-%d'),
                     tagged_prot=group.tagged_prot.iloc[0], depletion=label,
                     condition=f'{group.tagged_prot.iloc[0]} | {label}',
                     filename_treatment_hint=evidence,
                     treatment_conflict=conflict, split=split,
                     exclusion_reasons=';'.join(reasons), nuclei=int(group.nucleus.nunique()),
                     frames=250, independent_biological_unit='unknown')
        movies.append(movie)
        for reason in reasons:
            issues.append(dict(movie_id=movie_id, filename=filename, issue=reason,
                               source_label=label, filename_hint=evidence,
                               resolution='unresolved; original values retained'))
    movies = pd.DataFrame(movies)
    result = df.copy()
    result.insert(0, 'source_row_id', np.arange(len(df)))
    result = result.merge(movies.drop(columns=['tagged_prot', 'depletion', 'nuclei', 'frames']),
                          on='filename', how='left', validate='many_to_one', sort=False)
    result['trace_id'] = result.movie_id + '_n' + result.nucleus.astype(str)
    result['time_since_first_frame_s'] = result.frame * 30
    result['reference_type'] = 'measured_processed_roi_fluorescence'
    result.to_parquet(out/'observations.parquet', index=False)
    curves = []
    for filename, group in df.groupby('filename', sort=True):
        movie = movies.set_index('filename').loc[filename]
        curve = author_curve(group)
        for key in ['movie_id', 'split', 'condition', 'date_raw']:
            curve[key] = movie[key]
        curves.append(curve)
    curves = pd.concat(curves, ignore_index=True)
    curves.to_parquet(out/'movie_curves.parquet', index=False)
    movies.to_csv(out/'movies.csv', index=False)
    pd.DataFrame(issues).to_csv(out/'source_issues.csv', index=False)
    units = result[['trace_id', 'movie_id', 'nucleus', 'date_raw', 'split', 'condition']].drop_duplicates()
    units.to_csv(out/'nuclei.csv', index=False)
    coverage = movies.pivot_table(index='condition', columns='split', values='movie_id', aggfunc='count', fill_value=0)
    coverage.to_csv(out/'condition_coverage.csv')
    split_record = dict(frozen_date='2026-10-29', stage='before downstream predictive fitting',
                        grouping='whole source acquisition day; nuclei and frames remain nested',
                        development=sorted(DEVELOPMENT_DAYS), reserved_test=sorted(RESERVED_DAYS),
                        quarantine=sorted(QUARANTINE_DAYS), calibration=None,
                        selection='metadata-only condition coverage; no outcome-based selection',
                        limitation='Public technical-day holdout, not blinded independent biological replication. WAPL - aux and PDS5 - aux have no eligible second day. No universal calibration split is available.')
    (out/'splits.json').write_text(json.dumps(split_record, indent=2)+'\n')
    audit = dict(source_sha256=sha256(source), observations=len(result), nuclei=len(units),
                 movies=len(movies), conditions=len(coverage),
                 duplicate_observation_ids=0, missing_source_values=0,
                 inconsistent_movie_controls=0, treatment_conflict_movies=int(movies.treatment_conflict.sum()),
                 malformed_date_movies=int(movies.acquisition_date.isna().sum()),
                 row_counts_by_split=result.groupby('split').size().to_dict(),
                 nuclei_by_split=units.groupby('split').size().to_dict(),
                 movies_by_split=movies.groupby('split').size().to_dict(),
                 original_values_preserved=bool(result.sort_values('source_row_id')[SOURCE_COLUMNS].reset_index(drop=True).equals(df)),
                 biological_removal_truth_available=False, exact_bleach_timestamp_available=False,
                 model_evaluations_run=0)
    (out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    manifest = dict(title='Cohesin cofactor dosage sets the rate of loop extrusion',
                    authors=['Nessim Louafi', 'Luca Giorgetti'], doi='10.17632/spp4x8jpdy.1',
                    source_url=SOURCE_URL, landing_page='https://data.mendeley.com/datasets/spp4x8jpdy/1',
                    license='CC-BY-4.0', original_file='combined_dataframe.csv',
                    sha256=sha256(source), checksum_basis='publisher SHA-256 from file information',
                    bytes=source.stat().st_size, processing_commit=(ROOT/'source/ifr_source_commit.txt').read_text().strip(),
                    changes='Added identifiers, timing, issue flags, day-grouped splits and derived movie normalization; source fluorescence and labels unchanged.')
    (out/'source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
