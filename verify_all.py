"""Finite verification of Theorem 1.1 for every degree 4 <= n <= 2047.

    python verify_all.py [--workers N]

Steps:
  1. self-tests (tests/selftest.py);
  2. replay of certificates/joint.json  (three-parameter part, degrees 4..26);
  3. replay of certificates/blocks.json (two-parameter part, degrees 27..2047);
  4. check that every degree 4..2047 is covered exactly once.

Only the restrictions (1), (2), (3), (4), (5) of Proposition 5.1 are used.
The result is written to reports/report.json.  Exit code 0 means PASS.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from borcea import certificate            # noqa: E402
from tests import selftest                # noqa: E402


def load(path):
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--joint', default=str(ROOT / 'certificates' / 'joint.json'))
    ap.add_argument('--blocks', default=str(ROOT / 'certificates' / 'blocks.json'))
    args = ap.parse_args()
    if not 1 <= args.workers <= 64:
        raise ValueError('workers must lie between 1 and 64')

    t0 = time.time()
    print(selftest.run(), flush=True)

    joint_doc, joint_sha = load(Path(args.joint))
    joint_rows = certificate.verify_joint(joint_doc, args.workers)
    for row in joint_rows:
        print(f"degree {row['n']:>2}: PASS, {row['regions']} regions", flush=True)

    blocks_doc, blocks_sha = load(Path(args.blocks))
    block_rows = certificate.verify_blocks(blocks_doc, args.workers)
    for row in block_rows:
        print(f"degrees {row['nlo']}..{row['nhi']}: PASS, {row['regions']} regions", flush=True)

    covered = [r['n'] for r in joint_rows]
    for r in block_rows:
        covered.extend(range(r['nlo'], r['nhi'] + 1))
    if sorted(covered) != list(range(4, 2048)) or len(covered) != len(set(covered)):
        raise ValueError('degrees 4..2047 are not covered exactly once')

    report = {
        'status': 'PASS',
        'degrees': [4, 2047],
        'restrictions_used': ['(1)', '(2)', '(3)', '(4)', '(5)'],
        'joint': {'sha256': joint_sha, 'regions': sum(r['regions'] for r in joint_rows), 'per_degree': joint_rows},
        'blocks': {'sha256': blocks_sha, 'regions': sum(r['regions'] for r in block_rows), 'per_block': block_rows},
        'seconds': round(time.time() - t0, 1),
    }
    (ROOT / 'reports').mkdir(exist_ok=True)
    (ROOT / 'reports' / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: (v if k not in ('joint', 'blocks') else {'sha256': v['sha256'], 'regions': v['regions']})
                      for k, v in report.items()}, indent=2))


if __name__ == '__main__':
    main()
