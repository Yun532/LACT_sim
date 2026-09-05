"""核对两份中文说明、执行Notebook及结果记录的可追溯性；不替代科学误差检验。"""
from pathlib import Path
import hashlib
import json
import re

import nbformat

ROOT=Path(__file__).resolve().parents[1]


def main():
    for name in ['SII_PHYSICS_ZH.md','SII_IMPLEMENTATION_ZH.md']:
        path=ROOT/'docs'/name
        content=path.read_text(encoding='utf-8')
        assert content.count('$$') % 2 == 0, name
        assert content.replace('$$','').count('$') % 2 == 0, name
        for target in re.findall(r'\]\(([^)]+)\)',content):
            if target.startswith(('http:','https:','#')):
                continue
            assert (path.parent/target.split('#')[0]).exists(), (name,target)
    notebook=nbformat.read(ROOT/'notebooks/sii_complete_waveform_report.ipynb',as_version=4)
    nbformat.validate(notebook)
    code=[cell for cell in notebook.cells if cell.cell_type=='code']
    assert all(cell.execution_count is not None for cell in code)
    assert [cell.execution_count for cell in code] == list(range(1,len(code)+1))
    assert not any(out.output_type=='error' for cell in code for out in cell.outputs)
    summary=json.loads((ROOT/'validation/sii_science/summary.json').read_text(encoding='utf-8'))
    for section in ['code_sha256_lf','derived_input_sha256_lf']:
        for path,expected in summary[section].items():
            content=(ROOT/path).read_bytes().replace(b'\r\n',b'\n')
            assert hashlib.sha256(content).hexdigest()==expected, path
    assert summary['scope']['notebook_execution']=='all cells executed in this run'
    assert all(value==0 for value in summary['missing_effects_disabled'].values())
    print(f'说明链接、公式分隔符和输入哈希通过；Notebook共{len(code)}个代码单元已顺序执行且无错误。')


if __name__=='__main__':
    main()
