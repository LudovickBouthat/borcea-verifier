"""Certificate generator (NOT part of the proof; its output is replayed by verify_all.py).

Given an input certificate, replays every leaf with the current checkers and
replaces each leaf that fails by an exact midpoint subdivision whose leaves all
pass.  Usage:

    python tools/refine.py joint  IN.json OUT.json  [--workers 2] [--maxdepth 12]
    python tools/refine.py blocks IN.json OUT.json  [--workers 2] [--maxdepth 12]

Leaves that cannot be resolved within the depth budget become
{"unresolved": true}; the verifier rejects any certificate containing one.
"""
import argparse, json, os, sys, time
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from borcea import joint, blocks
from borcea.certificate import JOINT_FORMAT, BLOCKS_FORMAT

PARTS = (12, 48, 192)
OPT = {}


# ---------- three-parameter part ----------

def jsplit(box, ax):
    mid = (box[2 * ax] + box[2 * ax + 1]) / 2
    lo, hi = list(box), list(box)
    lo[2 * ax + 1] = mid
    hi[2 * ax] = mid
    return tuple(lo), tuple(hi)


def jleaf(n, box, parts_first):
    for p in dict.fromkeys((parts_first,) + PARTS):
        name = joint.check(n, box, p)
        if name is not None:
            return {'test': name, 'parts': p}
    return None


def jrefine(n, box, depth, budget, parts_first=12):
    leaf = jleaf(n, box, parts_first)
    if leaf is not None:
        return leaf
    if depth >= OPT['maxdepth'] or budget[0] <= 0:
        return None
    budget[0] -= 1
    root = (F(n - 1), F(1), F(1))
    order = sorted(range(3), key=lambda ax: -(box[2 * ax + 1] - box[2 * ax]) / root[ax])
    for ax in (order[:2] if depth == 0 else order[:1]):
        lo, hi = jsplit(box, ax)
        x = jrefine(n, lo, depth + 1, budget)
        if x is None:
            continue
        y = jrefine(n, hi, depth + 1, budget)
        if y is None:
            continue
        return {'axis': ax, 'children': [x, y]}
    return None


def jdegree(block):
    n = block['n']
    stats = {'kept': 0, 'refined': 0, 'unresolved': 0}

    def walk(node, box):
        if 'children' in node:
            lo, hi = jsplit(box, node['axis'])
            return {'axis': node['axis'], 'children': [walk(node['children'][0], lo), walk(node['children'][1], hi)]}
        parts = node.get('parts', 12)
        name = joint.check(n, box, parts) if 'test' in node else None
        if name is not None:
            stats['kept'] += 1
            return {'test': name, 'parts': parts}
        sub = jrefine(n, box, 0, [OPT['budget']], parts)
        if sub is None:
            stats['unresolved'] += 1
            return {'unresolved': True}
        stats['refined'] += 1
        return sub

    t = time.time()
    tree = walk(block['tree'], (F(0), F(n - 1), F(0), F(1), F(0), F(1)))
    print(n, stats, '%.0fs' % (time.time() - t), flush=True)
    return {'n': n, 'tree': tree}, stats


# ---------- two-parameter part ----------

def bsplit(region, axis):
    nl, nh, al, ah, el, eh = region
    if axis == 'A':
        mid = (al + ah) / 2
        return (nl, nh, al, mid, el, eh), (nl, nh, mid, ah, el, eh)
    if axis == 'e':
        mid = (el + eh) / 2
        return (nl, nh, al, ah, el, mid), (nl, nh, al, ah, mid, eh)
    mid = (nl + nh) // 2
    return (nl, mid, al, ah, el, eh), (mid + 1, nh, al, ah, el, eh)


def bleaf(region, mesh, parts_first):
    for p in dict.fromkeys((parts_first,) + PARTS):
        name = blocks.check(*region, p, mesh)
        if name is not None:
            return {'test': name, 'parts': p}
    return None


def brefine(region, mesh, depth, budget, parts_first=12):
    leaf = bleaf(region, mesh, parts_first)
    if leaf is not None:
        return leaf
    if depth >= OPT['maxdepth'] or budget[0] <= 0:
        return None
    budget[0] -= 1
    nl, nh, al, ah, el, eh = region
    axes = ['A', 'e'] if (ah - al) / max(F(1), F(nh - 1)) >= (eh - el) else ['e', 'A']
    if nh > nl:
        axes.append('n')
    for axis in (axes if depth == 0 else axes[:1]):
        lo, hi = bsplit(region, axis)
        x = brefine(lo, mesh, depth + 1, budget)
        if x is None:
            continue
        y = brefine(hi, mesh, depth + 1, budget)
        if y is None:
            continue
        return {'axis': axis, 'children': [x, y]}
    return None


def bblock(block):
    mesh = block.get('quadrature', 'uniform')
    stats = {'kept': 0, 'refined': 0, 'unresolved': 0}

    def walk(node, region):
        if 'children' in node:
            lo, hi = bsplit(region, node['axis'])
            return {'axis': node['axis'], 'children': [walk(node['children'][0], lo), walk(node['children'][1], hi)]}
        parts = node.get('parts', 12)
        name = blocks.check(*region, parts, mesh) if 'test' in node else None
        if name is not None:
            stats['kept'] += 1
            return {'test': name, 'parts': parts}
        sub = brefine(region, mesh, 0, [OPT['budget']], parts)
        if sub is None:
            stats['unresolved'] += 1
            return {'unresolved': True}
        stats['refined'] += 1
        return sub

    t = time.time()
    nlo, nhi = block['nlo'], block['nhi']
    tree = walk(block['tree'], (nlo, nhi, F(0), F(nhi - 1), F(0), F(1)))
    print(nlo, nhi, stats, '%.0fs' % (time.time() - t), flush=True)
    out = {'nlo': nlo, 'nhi': nhi, 'tree': tree}
    if mesh != 'uniform':
        out['quadrature'] = mesh
    return out, stats


def init(opt):
    OPT.update(opt)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('part', choices=('joint', 'blocks'))
    ap.add_argument('inp')
    ap.add_argument('out')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--maxdepth', type=int, default=12)
    ap.add_argument('--budget', type=int, default=20000)
    args = ap.parse_args()
    opt = {'maxdepth': args.maxdepth, 'budget': args.budget}
    init(opt)
    doc = json.load(open(args.inp))
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init, initargs=(opt,)) as pool:
        if args.part == 'joint':
            items = sorted(doc['degrees'], key=lambda b: -len(json.dumps(b['tree'])))
            res = list(pool.map(jdegree, items))
            out = {'format': JOINT_FORMAT, 'degrees': sorted((r[0] for r in res), key=lambda b: b['n'])}
        else:
            blocks_in = [b for b in doc['blocks'] if 'tree' in b and b['nlo'] >= 27]
            res = list(pool.map(bblock, blocks_in))
            out = {'format': BLOCKS_FORMAT, 'blocks': sorted((r[0] for r in res), key=lambda b: b['nlo'])}
    total = {k: sum(r[1][k] for r in res) for k in ('kept', 'refined', 'unresolved')}
    json.dump(out, open(args.out, 'w'), separators=(',', ':'))
    print('TOTAL', total, '->', args.out)
