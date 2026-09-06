"""读取CORSIKA TELESCOPE input，沿用现有NWU厘米到ENU米的转换。"""
from pathlib import Path
import argparse

import numpy as np
import pandas as pd


def read_corsika_layout(path):
    """镜号按TELESCOPE出现顺序编号；保留标签、原始坐标、半径及input行号。

    与EventIOPhotonSource.cpp及CorsikaTraceEventIOInput.cpp一致：
    北/西/上从cm转m后，ENU=(-西, 北, 上)。原点和高度来自input，不再平移。
    """
    rows = []
    for line_number, line in enumerate(Path(path).read_text(encoding='utf-8-sig').splitlines(), 1):
        fields = line.split('#', 1)[0].split()
        if not fields or fields[0].upper() != 'TELESCOPE':
            continue
        if len(fields) not in (5, 6):
            raise ValueError(f'第{line_number}行需要四个数值和可选镜名')
        try:
            north, west, up, radius = [float(item.upper().replace('D', 'E')) for item in fields[1:5]]
        except ValueError as error:
            raise ValueError(f'第{line_number}行TELESCOPE数值无效') from error
        if not np.all(np.isfinite([north, west, up, radius])) or radius <= 0:
            raise ValueError(f'第{line_number}行坐标必须有限、半径必须为正')
        index = len(rows)
        rows.append(dict(telescope_id=index+1, lactsim_index=index,
                         name=fields[5] if len(fields) == 6 else f'Tel.{index+1}',
                         corsika_north_cm=north, corsika_west_cm=west,
                         corsika_up_cm=up, radius_cm=radius,
                         east_m=-west*.01, north_m=north*.01, up_m=up*.01,
                         radius_m=radius*.01, source_input_line=line_number))
    frame = pd.DataFrame(rows)
    if len(frame) < 2:
        raise ValueError('input至少需要两台望远镜')
    if frame.name.duplicated().any():
        raise ValueError('TELESCOPE标签重复')
    return frame


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    layout = read_corsika_layout(args.input)
    layout.to_csv(args.output, index=False, float_format='%.5f')
    print(f'{len(layout)}台望远镜已转换至ENU米：{args.output}')
