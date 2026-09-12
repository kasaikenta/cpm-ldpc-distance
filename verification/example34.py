"""Check the J=3,L=4 example and generate its bitmap matrix and TikZ graph.

Uses only the Python standard library. Run from the repository root:
    python3 verification/example34.py
The rank, syndrome, pairing, and graph checks do not establish the distance
lower bound; independent24.py enumerates all nonzero codewords for that check.
"""
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path
import json
import struct
import zlib

ROOT = Path(__file__).resolve().parents[1]
# Four regular hexagons are the cycles formed by colors 0 and 1.
# Coordinates are indexed by the lexicographically ordered support.
# Selected color-2 edges curve around nonincident vertex circles.
GRAPH_POSITIONS = [
    [
        -2.97861,
        2.565
    ],
    [
        -1.02139,
        1.435
    ],
    [
        2.97861,
        1.435
    ],
    [
        1.02139,
        2.565
    ],
    [
        2.0,
        -3.13
    ],
    [
        2.0,
        -0.87
    ],
    [
        -1.02139,
        -2.565
    ],
    [
        -2.97861,
        -1.435
    ],
    [
        -2.0,
        3.13
    ],
    [
        -2.0,
        0.87
    ],
    [
        2.97861,
        -2.565
    ],
    [
        1.02139,
        -1.435
    ],
    [
        -2.0,
        -3.13
    ],
    [
        -2.0,
        -0.87
    ],
    [
        -1.02139,
        2.565
    ],
    [
        -2.97861,
        1.435
    ],
    [
        2.97861,
        2.565
    ],
    [
        1.02139,
        1.435
    ],
    [
        -2.97861,
        -2.565
    ],
    [
        -1.02139,
        -1.435
    ],
    [
        2.0,
        3.13
    ],
    [
        2.0,
        0.87
    ],
    [
        2.97861,
        -1.435
    ],
    [
        1.02139,
        -2.565
    ]
]
GRAPH_CURVES = {
    "0,18": [
        [
            -4.178610000000001,
            0.855
        ],
        [
            -4.178610000000001,
            -0.855
        ]
    ],
    "2,6": [
        [
            2.210962091615905,
            -0.4640187582825712
        ],
        [
            0.029100620858714676,
            -0.9488239541920475
        ]
    ],
    "3,7": [
        [
            -0.029100620858714232,
            0.9488239541920477
        ],
        [
            -2.2109620916159045,
            0.46401875828257144
        ]
    ],
    "14,22": [
        [
            0.029100620858714232,
            0.9488239541920477
        ],
        [
            2.2109620916159045,
            0.46401875828257144
        ]
    ],
    "15,23": [
        [
            -2.210962091615905,
            -0.4640187582825712
        ],
        [
            -0.029100620858714676,
            -0.9488239541920475
        ]
    ],
    "10,16": [
        [
            4.178610000000001,
            -0.855
        ],
        [
            4.178610000000001,
            0.855
        ]
    ]
}


def build_example():
    table = json.loads((ROOT / 'joint_examples.json').read_text())
    case = next(c for c in table if (c['J'], c['L']) == (3, 4))
    E, P = case['E'], case['P']
    assert E == [[0, 0, 0, 0], [0, 1, 2, 4], [0, 7, 15, 5]] and P == 24
    J, L = len(E), len(E[0])
    H = [[int(s == (t + E[j][ell]) % P)
          for ell in range(L) for t in range(P)]
         for j in range(J) for s in range(P)]
    assert {sum(row) for row in H} == {L}
    assert {sum(row[col] for row in H) for col in range(L * P)} == {J}
    pivots = {}
    for row in H:
        word = sum(value << k for k, value in enumerate(row))
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                break
    assert len(pivots) == 70 and case['k'] == L * P - len(pivots) == 26
    support = []
    for ell in range(L):
        terms = [sum(E[j][perm[j]] for j in range(J)) % P
                 for perm in permutations([k for k in range(L) if k != ell])]
        assert len(set(terms)) == 6
        support.extend((ell, t) for t in sorted(terms))
    columns = [ell * P + t for ell, t in support]
    row_counts = [sum(row[col] for col in columns) for row in H]
    assert Counter(row_counts) == {0: 38, 2: 32, 4: 2}
    assert all(count % 2 == 0 for count in row_counts)
    checks = defaultdict(list)
    for v, (ell, t) in enumerate(support):
        for j in range(J):
            checks[j, (t + E[j][ell]) % P].append(v)
    edges = []
    for (j, s), vertices in sorted(checks.items()):
        # Consecutive vertices in lexicographic order specify the pairing Pi.
        for k in range(0, len(vertices), 2):
            u, v = vertices[k:k + 2]
            edges.append(dict(u=u, v=v, j=j, s=s))
            lu, tu = support[u]
            lv, tv = support[v]
            assert lu != lv
            assert (tv - tu - E[j][lu] + E[j][lv]) % P == 0
    assert len(edges) == 36
    neighbors = {v: set() for v in range(24)}
    colors = {v: [] for v in range(24)}
    for edge in edges:
        u, v, j = edge['u'], edge['v'], edge['j']
        neighbors[u].add(v)
        neighbors[v].add(u)
        colors[u].append(j)
        colors[v].append(j)
    assert all(sorted(value) == [0, 1, 2] for value in colors.values())
    reached, pending = {0}, [0]
    while pending:
        for v in neighbors[pending.pop()] - reached:
            reached.add(v)
            pending.append(v)
    assert len(reached) == 24
    return dict(J=J, L=L, P=P, E=E, n=96, k=26, rank=70,
                support_by_block=[[t for ell, t in support if ell == l] for l in range(L)],
                support=support, support_column_indices=columns, weight=24,
                row_weight_counts=dict(sorted(Counter(row_counts).items())),
                edges=edges, vertices=24, edge_count=36, connected=True,
                pairing_rule='Consecutive support positions in lexicographic order at each check',
                all_checks_passed=True), H


def write_png(path, cells, scale=8, block_size=None, highlighted_columns=None):
    """Lossless RGB bitmap; highlighted 1s share the graph's row colors."""
    width, height = len(cells[0]) * scale, len(cells) * scale
    colors = [(0, 114, 178), (213, 94, 0), (0, 158, 115)]
    selected = set(highlighted_columns or [])
    raw = bytearray()
    for y in range(height):
        raw.append(0)  # PNG scanline filter: none.
        for x in range(width):
            row, col = y // scale, x // scale
            value = (255, 255, 255)
            if cells[row][col]:
                if highlighted_columns is None:
                    value = (0, 0, 0)
                else:
                    value = colors[row // 24] if col in selected else (185, 185, 185)
            if block_size and ((x and x % (block_size * scale) == 0)
                               or (y and y % (block_size * scale) == 0)):
                value = (165, 165, 165)
            raw.extend(value)

    def chunk(name, data):
        return struct.pack('!I', len(data)) + name + data + struct.pack('!I', zlib.crc32(name + data) & 0xffffffff)

    png = (b'\x89PNG\r\n\x1a\n'
           + chunk(b'IHDR', struct.pack('!IIBBBBB', width, height, 8, 2, 0, 0, 0))
           + chunk(b'IDAT', zlib.compress(bytes(raw), 9)) + chunk(b'IEND', b''))
    path.write_bytes(png)


# Each support position has a unique (block color, within-block shape) pair.
MARKER_COLORS = [(0, 90, 160), (190, 70, 0), (0, 125, 90), (155, 55, 140)]
MARKER_PATHS = [
    '(0,0) circle[radius=0.5]',
    '(-0.5,-0.5) rectangle (0.5,0.5)',
    '(0,0.5) -- (-0.5,-0.5) -- (0.5,-0.5) -- cycle',
    '(0,0.5) -- (-0.5,0) -- (0,-0.5) -- (0.5,0) -- cycle',
    '(0,-0.5) -- (-0.5,0.5) -- (0.5,0.5) -- cycle',
    '(-0.18,0.5) -- (0.18,0.5) -- (0.18,0.18) -- (0.5,0.18) -- (0.5,-0.18) -- (0.18,-0.18) -- (0.18,-0.5) -- (-0.18,-0.5) -- (-0.18,-0.18) -- (-0.5,-0.18) -- (-0.5,0.18) -- (-0.18,0.18) -- cycle',
]


def marker_tex(v, position, size):
    # Both panels call exactly the same path/color renderer; only size differs.
    color, shape = divmod(v, 6)
    return (rf'\begin{{scope}}[shift={{{position}}},x={size},y={size}]' + '\n'
            + rf'\path[fill=exampleV{color}] {MARKER_PATHS[shape]};' + '\n'
            + r'\end{scope}')


def main():
    data, H = build_example()
    out = ROOT / 'figures'
    out.mkdir(exist_ok=True)
    write_png(out / 'example34_matrix.png', H, block_size=24,
              highlighted_columns=data['support_column_indices'])
    word = [int(k in data['support_column_indices']) for k in range(96)]
    write_png(out / 'example34_codeword.png', [word] * 3, block_size=24)
    (out / 'example34.json').write_text(json.dumps(data, indent=2) + '\n')
    supports = [r'\begin{aligned}']
    for ell, positions in enumerate(data['support_by_block']):
        supports.append(f'T_{ell}&=' + r'\{' + ','.join(map(str, positions)) + r'\}' + (r',\\' if ell < 3 else '.'))
    supports.append(r'\end{aligned}')
    (out / 'example34_support.tex').write_text('\n'.join(supports) + '\n')
    definitions = [rf'\definecolor{{exampleV{ell}}}{{RGB}}{{{r},{g},{b}}}'
                   for ell, (r, g, b) in enumerate(MARKER_COLORS)]
    (out / 'example34_symbols.tex').write_text('\n'.join(definitions) + '\n')
    strip = [r'% Nonzero entries are the same symbols as the graph vertices.',
             r'\draw[black!45,thin] (0,-5.45) rectangle (6.8,-5.6625);']
    for boundary in [1.7, 3.4, 5.1]:
        strip.append(rf'\draw[black!45,thin] ({boundary},-5.45) -- ({boundary},-5.6625);')
    for v, column in enumerate(data['support_column_indices']):
        x = (column + 0.5) * 6.8 / 96
        strip.append(marker_tex(v, f'({x:.8f},-5.55625)', '0.62mm'))
    (out / 'example34_codeword.tex').write_text('\n'.join(strip) + '\n')
    assert len(set((v // 6, v % 6) for v in range(24))) == 24
    assert len(GRAPH_POSITIONS) == 24
    lines = [r'\input{figures/example34_symbols.tex}', r'\definecolor{exampleJzero}{RGB}{0,114,178}',
             r'\definecolor{exampleJone}{RGB}{213,94,0}',
             r'\definecolor{exampleJtwo}{RGB}{0,158,115}',
             r'\begin{tikzpicture}[x=1cm,y=1cm,',
             r'  g0/.style={exampleJzero,line width=0.7pt},',
             r'  g1/.style={exampleJone,dashed,line width=0.85pt},',
             r'  g2/.style={exampleJtwo,densely dotted,line width=1pt,preaction={draw=white,solid,line width=2.2pt}},',
             r'  vertexlabel/.style={font=\scriptsize,inner sep=0pt,text=black}]']
    for v, (x, y) in enumerate(GRAPH_POSITIONS):
        lines.append(f'\\coordinate (v{v}) at ({x:.5f},{y:.5f});')
    for edge in data['edges']:
        u, v, j = edge['u'], edge['v'], edge['j']
        controls = GRAPH_CURVES.get(f'{u},{v}')
        if controls:
            a, b = controls
            lines.append(f'\\draw[g{j}] (v{u}) .. controls ({a[0]:.5f},{a[1]:.5f}) and ({b[0]:.5f},{b[1]:.5f}) .. (v{v});')
        else:
            lines.append(f'\\draw[g{j}] (v{u}) -- (v{v});')
    for v, (ell, t) in enumerate(data['support']):
        lines.append(f'\\fill[white] (v{v}) circle[radius=1.1mm];')
        lines.append(marker_tex(v, f'(v{v})', '2mm'))
        x, y = GRAPH_POSITIONS[v]
        cx, cy = (2.0 if x > 0 else -2.0), (2.0 if y > 0 else -2.0)
        # Put the small label 4 mm inward, in the empty hexagon interior.
        radius = ((x - cx)**2 + (y - cy)**2)**0.5
        lx, ly = x + 0.4 * (cx - x) / radius, y + 0.4 * (cy - y) / radius
        lines.append(f'\\node[vertexlabel] at ({lx:.5f},{ly:.5f}) {{$({ell},{t})$}};')
    for j, x in enumerate([-2.7, -0.7, 1.3]):
        lines.append(f'\\draw[g{j}] ({x},-4.2) -- ({x+.6},-4.2) node[right,black,font=\\small] {{$j={j}$}};')
    lines.append(r'\end{tikzpicture}')
    (out / 'example34_graph.tex').write_text('\n'.join(lines) + '\n')
    print(json.dumps({k: data[k] for k in ['J', 'L', 'P', 'rank', 'weight', 'vertices', 'edge_count', 'connected', 'row_weight_counts', 'all_checks_passed']}))


if __name__ == '__main__':
    main()
