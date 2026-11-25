#!/usr/bin/env python3
"""Dump all text from a set of .pptx decks via the local editor_sdk MCP endpoint."""
import json
import os
import sys
import time

sys.path.insert(0, "/opt/WorkBuddy/app.asar.unpacked/resources/plugins/workbuddy-builtin/skills/tencent-local-office-edit")
import edsdk  # noqa: E402

FILES = [
    "/home/LittleWatter/Downloads/01-自修-生物学基础.pptx",
    "/home/LittleWatter/Downloads/02生物序列数据检索和获取2026秋.pptx",
    "/home/LittleWatter/Downloads/03双序列比对-1-全局比对2026秋.pptx",
    "/home/LittleWatter/Downloads/03双序列比对-2-局部比对2026秋.pptx",
    "/home/LittleWatter/Downloads/03双序列比对-3-启发式算法2026版.pptx",
]

OUT_DIR = "/home/LittleWatter/Documents/Obsidian Vault/大三上/.workbuddy/tmp/bi"
os.makedirs(OUT_DIR, exist_ok=True)
EP = edsdk._discover_endpoint()


def call(name, args, timeout=120):
    """Raw MCP call; returns parsed text payload or None on error."""
    obj = edsdk._request(EP, "tools/call", {"name": name, "arguments": args}, timeout=timeout)
    if obj.get("error"):
        return None, obj["error"].get("message", "err")
    res = obj.get("result") or {}
    if res.get("isError"):
        msg = ""
        for c in res.get("content") or []:
            msg += c.get("text", "")
        return None, msg
    for c in res.get("content") or []:
        try:
            return json.loads(c.get("text", "")), None
        except Exception:  # noqa: BLE001
            return c.get("text"), None
    return None, "empty"


def walk(shapes, acc):
    for sh in shapes or []:
        t = sh.get("text")
        if t and str(t).strip():
            acc.append((sh.get("y", 0), sh.get("x", 0), str(t).strip()))
        if sh.get("children"):
            walk(sh["children"], acc)


for idx, path in enumerate(FILES, 1):
    name = os.path.basename(path)
    call("open_file", {"file_path": path})
    n = 0
    for _ in range(120):
        info, err = call("slide_get_info", {"file_id": path})
        if isinstance(info, dict) and info.get("slide_count"):
            n = info["slide_count"]
            break
        time.sleep(1)
    lines = [f"# DECK {idx}: {name}  ({n} slides)", ""]
    print(f"deck {idx} {name}: {n} slides", flush=True)
    for p in range(n):
        page, perr = call("slide_get_page_info", {"file_id": path, "page_index": p})
        if perr or not isinstance(page, dict):
            lines.append(f"===== SLIDE {p + 1} ===== [page_info failed: {perr}]")
            continue
        acc = []
        walk(page.get("shapes"), acc)
        acc.sort(key=lambda z: (z[0], z[1]))
        body = "\n".join(a[2] for a in acc)
        # tables: pull cell text explicitly
        tbl_txt = []
        for sh in page.get("shapes") or []:
            if sh.get("type") == "table" or sh.get("geom") == "table":
                tb, terr = call("slide_get_table_info", {"file_id": path, "page_index": p, "shape_id": sh.get("id")})
                if isinstance(tb, dict):
                    for c in tb.get("cells") or []:
                        if c.get("is_merge_continue"):
                            continue
                        if (c.get("text") or "").strip():
                            tbl_txt.append(c["text"].strip())
        note, nerr = call("slide_get_notes_text", {"file_id": path, "page_index": p})
        lines.append(f"===== SLIDE {p + 1} =====")
        lines.append(body if body else "(no text shapes)")
        if tbl_txt:
            lines.append("[TABLE] " + " | ".join(tbl_txt))
        if isinstance(note, dict) and (note.get("text") or "").strip():
            lines.append("[NOTES] " + note["text"].strip())
        lines.append("")
    out = os.path.join(OUT_DIR, f"{idx:02d}.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"OK {out} chars={sum(len(x) for x in lines)}", flush=True)
