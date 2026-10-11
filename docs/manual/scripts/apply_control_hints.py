"""Proposes (and with --apply, inserts) tooltips for dialog widgets from the
manual's control reference (\\item[Label.] description).

python docs/manual/scripts/apply_control_hints.py            # print the proposed mapping
python docs/manual/scripts/apply_control_hints.py --apply    # wrap the widgets with dialogLayout::hinted
"""
import re
import sys

CHAPTERS = ['docs/manual/chapters/02-builder-controls.tex', 'docs/manual/chapters/03-edit-controls.tex']
# Manual subsection title (start) -> dialog source.
SECTIONS = {
    'Bulk Crystal': 'src/ui/BulkCrystalBuilderDialog.cpp',
    'Substitutional Solid Solution': 'src/ui/SubstitutionalSolidSolutionDialog.cpp',
    'CSL Grain Boundary': 'src/ui/CSLGrainBoundaryDialog.cpp',
    'Nanocrystal': 'src/ui/NanoCrystalBuilderDialog.cpp',
    'Custom Structure': 'src/ui/CustomStructureDialog.cpp',
    'Polycrystal': 'src/ui/PolyCrystalBuilderDialog.cpp',
    'Amorphous Structure': 'src/ui/AmorphousBuilderDialog.cpp',
    'Insert Dislocation': 'src/ui/DislocationBuilderDialog.cpp',
    'Add Interstitial': 'src/ui/InterstitialAtomsDialog.cpp',
    'Cell Sculptor': 'src/ui/CellSculptorDialog.cpp',
    'Merge Structures': 'src/ui/MergeStructuresDialog.cpp',
}
WIDGETS = ['Checkbox', 'InputInt', 'InputFloat', 'InputDouble', 'DragFloat', 'DragInt', 'SliderFloat', 'SliderInt',
           'InputFloat3', 'DragFloat3', 'InputInt3', 'Combo', 'RadioButton']


def clean(tex):
    t = tex.replace('\\,', ' ')
    t = re.sub(r'\\(ui|code|textbf|emph|textit)\{([^}]*)\}', r'\2', t)
    t = t.replace('\\angstrom\\ ', 'Angstrom ').replace('\\angstrom', 'Angstrom').replace('--', '-').replace('\\%', '%')
    t = re.sub(r'\$([^$]*)\$', r'\1', t)
    t = re.sub(r'\\[a-zA-Z]+\s*', '', t).replace('{', '').replace('}', '').replace('~', ' ')
    t = re.sub(r'\s+', ' ', t).strip()
    sentences = re.split(r'(?<=[.!?])\s+', t)
    out = ''
    for s in sentences:
        if len(out) + len(s) > 260 and out:
            break
        out = (out + ' ' + s).strip()
    return out


def norm(label):
    label = label.split('##')[0]
    label = re.sub(r'\(.*?\)', '', label)
    label = re.sub(r'[^a-z0-9 ]', ' ', label.lower())
    return re.sub(r'\s+', ' ', label).strip()


def items():
    result = {}
    for chapter in CHAPTERS:
        text = open(chapter, encoding='utf-8').read()
        for part in re.split(r'\\subsection\{', text)[1:]:
            title = part[:part.index('}')]
            source = next((src for key, src in SECTIONS.items() if title.startswith(key) or ('Optional ' + key) in title
                           or key in title.split(':')[0]), None)
            if not source:
                continue
            for m in re.finditer(r'\\item\[(.*?)\.\]\s*(.*?)(?=\\item\[|\\end\{description\}|\\subsection|$)', part, re.S):
                labels = [norm(x) for x in re.split(r'/|,| and ', re.sub(r'\\[a-z]+\{|\}', '', m.group(1)))]
                result.setdefault(source, []).append(([l for l in labels if l], clean(m.group(2))))
    return result


def calls(source):
    """(start, end, label) of each widget call with a literal label."""
    s = open(source, encoding='utf-8').read()
    for m in re.finditer(r'ImGui::(%s)\("([^"]*)"' % '|'.join(WIDGETS), s):
        depth, i = 0, m.start() + len('ImGui::') + len(m.group(1))
        while True:
            c = s[i]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    break
            elif c == '"':
                i = s.index('"', i + 1)
            i += 1
        yield m.start(), i + 1, m.group(2)


def main(apply):
    table = items()
    for source, entries in table.items():
        s = open(source, encoding='utf-8').read()
        edits = []
        for start, end, label in calls(source):
            if s[max(0, start - 22):start].endswith('dialogLayout::hinted('):
                continue
            key = norm(label)
            if not key:
                continue
            match = next((text for labels, text in entries if key in labels), None)
            if match is None:
                match = next((text for labels, text in entries if len(key) >= 3 and any(len(l) >= 3 and (l.startswith(key) or key.startswith(l)) for l in labels)), None)
            if match:
                edits.append((start, end, label, match))
        print('==', source, len(edits))
        for _, _, label, text in edits:
            print('  %-40s %s' % (label[:40], text[:110]))
        if apply:
            for start, end, label, text in sorted(edits, reverse=True):
                escaped = text.replace('\\', '\\\\').replace('"', '\\"')
                s = s[:start] + 'dialogLayout::hinted(' + s[start:end] + ', "' + escaped + '")' + s[end:]
            if edits and '#include "ui/DialogLayout.h"' not in s:
                s = s.replace('#include "ui/ResponsiveLayout.h"\n', '#include "ui/ResponsiveLayout.h"\n#include "ui/DialogLayout.h"\n', 1)
            open(source, 'w', encoding='utf-8', newline='\n').write(s)


main('--apply' in sys.argv)
