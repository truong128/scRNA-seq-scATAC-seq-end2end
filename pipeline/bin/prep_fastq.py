#!/usr/bin/env python3

import argparse
import tarfile
from pathlib import Path
from download import get


def unpack(archive, where):
    where.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'r:*') as tar:

        for m in tar.getmembers():
            dest = (where / m.name).resolve()
            if not dest.is_relative_to(where.resolve()) or m.islnk() or m.issym():
                raise ValueError('Unsafe tar path or link: ' + m.name)
        try:
            tar.extractall(where, filter='data')
        except TypeError:

            tar.extractall(where)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fastq-url', required=True)
    p.add_argument('--ref-url', required=True)
    args = p.parse_args()
    get(args.fastq_url, Path('fastq.tar'), 1000000)
    unpack('fastq.tar', Path('fastqs'))
    Path('fastq.tar').unlink()
    get(args.ref_url, Path('reference.tar.gz'), 1000000)
    unpack('reference.tar.gz', Path('reference'))
    Path('reference.tar.gz').unlink()
    fq = list(Path('fastqs').rglob('*R1*fastq.gz'))
    ref = list(Path('reference').rglob('reference.json'))
    if not fq or len(ref) != 1:
        raise RuntimeError(f'Invalid unpack result: {len(fq)} R1 files; {len(ref)} references')
    Path('reference_path.txt').write_text(str(ref[0].parent.resolve()))
    Path('fastq_path.txt').write_text(str(fq[0].parent.resolve()))


if __name__ == '__main__':
    main()
