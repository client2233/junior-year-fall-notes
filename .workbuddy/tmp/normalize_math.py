#!/usr/bin/env python3
"""泛函分析笔记：公式块规范化
1. 所有范数竖线 \|  -> \lVert ... \rVert（消除 markdown 转义歧义）
2. 所有单行 $$...$$ 展开为标准三行块
3. 所有块级公式前后补空行
先备份到 /tmp/fa_backup，并做偶数配对校验。
"""
import glob
import os
import re
import shutil
import sys

D = "/home/LittleWatter/Documents/Obsidian Vault/大三上/泛函分析"
BK = "/tmp/fa_backup"
os.makedirs(BK, exist_ok=True)


def norm_pipes(seg):
    """把段内 \\| 按出现顺序交替换成 \\lVert / \\rVert"""
    cnt = [0]

    def rep(_m):
        cnt[0] += 1
        return r'\lVert' if cnt[0] % 2 else r'\rVert'

    return re.sub(r'\\\|', rep, seg)


def norm_inline_line(line):
    """只处理行内 $...$ 段"""
    res = []
    pos = 0
    for m in re.finditer(r'(?<!\$)\$(?!\$)(?:[^$\n]|\$(?!\$))*?(?<!\$)\$(?!\$)', line):
        res.append(line[pos:m.start()])
        res.append(norm_pipes(m.group(0)))
        pos = m.end()
    res.append(line[pos:])
    return "".join(res)


def process(text, fname, report):
    lines = text.split("\n")
    n = len(lines)

    # ---------- 第一步：\|\| 规范化 ----------
    i = 0
    while i < n:
        st = lines[i].strip()
        if st.startswith("$$"):
            if st.count("$$") >= 2:                      # 单行块
                inner = st[2:-2]
                odd = inner.count(r'\|') % 2
                if odd:
                    report.append(f"  !! {fname}:{i+1} 单行块内 \\| 为奇数，跳过")
                else:
                    indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
                    lines[i] = indent + "$$" + norm_pipes(inner) + "$$"
                i += 1
            else:                                        # 多行块
                j = i + 1
                while j < n and not lines[j].strip().endswith("$$"):
                    j += 1
                seg = "\n".join(lines[i:j + 1])
                if seg.count(r'\|') % 2:
                    report.append(f"  !! {fname}:{i+1}-{j+1} 多行块内 \\| 为奇数，跳过")
                else:
                    cnt = [0]

                    def rep(_m):
                        cnt[0] += 1
                        return r'\lVert' if cnt[0] % 2 else r'\rVert'
                    for k in range(i, j + 1):
                        lines[k] = re.sub(r'\\\|', rep, lines[k])
                i = j + 1
        else:
            lines[i] = norm_inline_line(lines[i])
            i += 1

    # ---------- 第二步：块展开 + 空行 ----------
    out = []
    i = 0
    while i < n:
        st = lines[i].strip()
        if st.startswith("$$"):
            indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
            if st.count("$$") >= 2:                      # 单行 -> 三行
                inner = st[2:-2].strip()
                if out and out[-1].strip() != "":
                    out.append("")
                out.append(indent + "$$")
                out.append(indent + inner)
                out.append(indent + "$$")
                if i + 1 < n and lines[i + 1].strip() != "":
                    out.append("")
                i += 1
            else:                                        # 多行块
                if out and out[-1].strip() != "":
                    out.append("")
                out.append(lines[i])
                j = i + 1
                while j < n:
                    out.append(lines[j])
                    if lines[j].strip().endswith("$$"):
                        break
                    j += 1
                if j + 1 < n and lines[j + 1].strip() != "":
                    out.append("")
                i = j + 1
        else:
            out.append(lines[i])
            i += 1

    # ---------- 第三步：压缩连续空行 ----------
    res = []
    for l in out:
        if l.strip() == "" and res and res[-1].strip() == "":
            continue
        res.append(l)
    return "\n".join(res)


report = []
for f in sorted(glob.glob(os.path.join(D, "*.md"))):
    name = os.path.basename(f)
    shutil.copy(f, os.path.join(BK, name))
    src = open(f, encoding="utf-8").read()
    dst = process(src, name, report)
    open(f, "w", encoding="utf-8").write(dst)
    print(f"{name}: {len(src.splitlines())} -> {len(dst.splitlines())} 行")

if report:
    print("\n".join(report))
else:
    print("\n配对校验：全部通过，无跳过")
