"""Replay of the two certificates (Section 6.1).

A certificate is a finite binary tree per degree (three-parameter part) or per
degree block (two-parameter part).  Every internal node bisects one coordinate
at its exact rational midpoint; every leaf names an exclusion and a quadrature
size.  The replay reconstructs every box from the root, recomputes an exclusion
at every leaf with the checkers of this package, and fails closed on anything
else.  Leaf names are hints only: a leaf is accepted when the checker finds
some exclusion on its box.  Induction on the tree shows that the leaves cover
the root, so a tree whose leaves are all accepted proves Lemma 6.1 for the
degrees it covers.
"""
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from . import joint, blocks

JOINT_FORMAT = 'borcea-joint-v2'
BLOCKS_FORMAT = 'borcea-blocks-v2'
MAX_DEPTH = 200


def _joint_degree(block):
    n = block['n']
    if type(n) is not int or not 4 <= n <= 26:
        raise ValueError('joint degree out of range')
    counts = {}

    def walk(node, box, depth):
        if depth > MAX_DEPTH:
            raise ValueError('depth')
        if set(node) == {'axis', 'children'}:
            ax = node['axis']
            if type(ax) is not int or ax not in (0, 1, 2) or len(node['children']) != 2:
                raise ValueError('invalid split')
            mid = (box[2 * ax] + box[2 * ax + 1]) / 2
            lo, hi = list(box), list(box)
            lo[2 * ax + 1] = mid
            hi[2 * ax] = mid
            walk(node['children'][0], tuple(lo), depth + 1)
            walk(node['children'][1], tuple(hi), depth + 1)
        elif set(node) == {'test', 'parts'}:
            parts = node['parts']
            if type(parts) is not int or not 8 <= parts <= 512:
                raise ValueError('invalid quadrature size')
            name = joint.check(n, box, parts)
            if name is None:
                raise ValueError(f'leaf not excluded: degree {n}, box {[str(x) for x in box]}')
            counts[name] = counts.get(name, 0) + 1
        else:
            raise ValueError('invalid node')

    walk(block['tree'], (F(0), F(n - 1), F(0), F(1), F(0), F(1)), 0)
    return {'n': n, 'regions': sum(counts.values()), 'exclusions': counts}


def verify_joint(doc, workers=1):
    if doc.get('format') != JOINT_FORMAT:
        raise ValueError('wrong joint certificate format')
    degrees = [b['n'] for b in doc['degrees']]
    if degrees != list(range(4, 27)):
        raise ValueError('the joint certificate must list degrees 4..26 in order')
    if workers == 1:
        rows = [_joint_degree(b) for b in doc['degrees']]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_joint_degree, doc['degrees']))
    return rows


def _blocks_block(block):
    nlo, nhi = block['nlo'], block['nhi']
    mesh = block.get('quadrature', 'uniform')
    if mesh not in ('uniform', 'geometric'):
        raise ValueError('unknown quadrature mesh')
    counts = {}

    def walk(node, nl, nh, al, ah, el, eh, depth):
        if depth > MAX_DEPTH:
            raise ValueError('depth')
        if set(node) == {'axis', 'children'}:
            ch = node['children']
            if not isinstance(ch, list) or len(ch) != 2:
                raise ValueError('invalid children')
            if node['axis'] == 'A':
                mid = (al + ah) / 2
                walk(ch[0], nl, nh, al, mid, el, eh, depth + 1)
                walk(ch[1], nl, nh, mid, ah, el, eh, depth + 1)
            elif node['axis'] == 'e':
                mid = (el + eh) / 2
                walk(ch[0], nl, nh, al, ah, el, mid, depth + 1)
                walk(ch[1], nl, nh, al, ah, mid, eh, depth + 1)
            elif node['axis'] == 'n':
                if nl == nh:
                    raise ValueError('cannot bisect a single degree')
                mid = (nl + nh) // 2
                walk(ch[0], nl, mid, al, ah, el, eh, depth + 1)
                walk(ch[1], mid + 1, nh, al, ah, el, eh, depth + 1)
            else:
                raise ValueError('invalid split axis')
        elif set(node) == {'test', 'parts'}:
            parts = node['parts']
            if type(parts) is not int or not 8 <= parts <= 4096:
                raise ValueError('invalid quadrature size')
            name = blocks.check(nl, nh, al, ah, el, eh, parts, mesh)
            if name is None:
                raise ValueError(f'leaf not excluded: degrees {nl}..{nh}, A {al}..{ah}, e {el}..{eh}')
            counts[name] = counts.get(name, 0) + 1
        else:
            raise ValueError('invalid node')

    walk(block['tree'], nlo, nhi, F(0), F(nhi - 1), F(0), F(1), 0)
    return {'nlo': nlo, 'nhi': nhi, 'quadrature': mesh, 'regions': sum(counts.values()), 'exclusions': counts}


def verify_blocks(doc, workers=1):
    if doc.get('format') != BLOCKS_FORMAT:
        raise ValueError('wrong block certificate format')
    expected = 27
    for b in doc['blocks']:
        if set(b) - {'nlo', 'nhi', 'quadrature', 'tree'}:
            raise ValueError('unexpected block keys')
        if b['nlo'] != expected or type(b['nhi']) is not int or not b['nlo'] <= b['nhi'] <= 2047:
            raise ValueError('degree blocks must be consecutive, without gap or overlap')
        expected = b['nhi'] + 1
    if expected != 2048:
        raise ValueError('the degree blocks must end at 2047')
    if workers == 1:
        rows = [_blocks_block(b) for b in doc['blocks']]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_blocks_block, doc['blocks']))
    return rows
