#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import requests


def get(url: str, dest: Path, minimum_bytes: int = 1000, tries: int = 4) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size >= minimum_bytes:
        return {'url': url, 'file': str(dest), 'bytes': dest.stat().st_size, 'cached': True}
    tmp = dest.with_name(dest.name + '.partial')
    if tmp.exists():
        tmp.unlink()
    last_error = None
    for attempt in range(tries):
        try:
            headers = {'User-Agent': 'ImmunoOncologyPilot/1.0 (research download)'}
            with requests.get(url, stream=True, timeout=(20, 90), headers=headers) as response:
                response.raise_for_status()
                size = int(response.headers.get('Content-Length', '0'))
                total = 0
                digest = hashlib.sha256()
                with tmp.open('wb') as fh:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            total += len(chunk)
                            digest.update(chunk)
                            fh.write(chunk)
                if size and total != size:
                    raise IOError(f'Incomplete download, expected {size} bytes, got {total}')
                if total < minimum_bytes:
                    raise IOError(f'File too small ({total} bytes), possible HTML or error response')
                os.replace(tmp, dest)
                return {'url': url, 'file': str(dest), 'bytes': total,
                        'sha256': digest.hexdigest(), 'cached': False}
        except Exception as exc:
            last_error = exc
            if tmp.exists():
                tmp.unlink()
            if attempt < tries - 1:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError(f'Download failed {url}: {last_error}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--url', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--minimum-bytes', type=int, default=1000)
    args = p.parse_args()
    meta = get(args.url, Path(args.output), args.minimum_bytes)
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    main()
