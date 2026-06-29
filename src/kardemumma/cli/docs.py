"""
generate_docs.py — Parse src/kardemumma and emit HTML + PDF function reference.

Usage:
    python3 generate_docs.py [--out-dir OUTDIR]

Outputs:
    <out-dir>/kardemumma_api.html
    <out-dir>/kardemumma_api.pdf   (requires google-chrome on PATH)
"""

from __future__ import annotations

import argparse
import ast
import inspect
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field
from html import escape
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class ParamDoc:
    name: str
    annotation: str = ""
    default: str = ""


@dataclass
class FunctionDoc:
    name: str
    is_private: bool
    params: list[ParamDoc]
    returns: str
    docstring: str
    lineno: int


@dataclass
class ClassDoc:
    name: str
    docstring: str
    methods: list[FunctionDoc]
    lineno: int


@dataclass
class ModuleDoc:
    name: str           # e.g. "prm"
    filepath: Path
    docstring: str
    functions: list[FunctionDoc]
    classes: list[ClassDoc]


# ---------------------------------------------------------------------------
# AST parsing helpers
# ---------------------------------------------------------------------------

def _annotation_str(node: Optional[ast.expr]) -> str:
    if node is None:
        return ""
    return ast.unparse(node)


def _default_str(node: Optional[ast.expr]) -> str:
    if node is None:
        return ""
    return ast.unparse(node)


def _parse_function(node: ast.FunctionDef | ast.AsyncFunctionDef) -> FunctionDoc:
    args = node.args
    all_args = args.posonlyargs + args.args + args.kwonlyargs
    if args.vararg:
        all_args.append(args.vararg)
    if args.kwarg:
        all_args.append(args.kwarg)

    # Build defaults map: align defaults to the END of positional args
    defaults_map: dict[str, str] = {}
    positional = args.posonlyargs + args.args
    n_defaults = len(args.defaults)
    for i, default in enumerate(args.defaults):
        arg = positional[len(positional) - n_defaults + i]
        defaults_map[arg.arg] = _default_str(default)
    for arg, default in zip(args.kwonlyargs, args.kw_defaults):
        if default is not None:
            defaults_map[arg.arg] = _default_str(default)

    params = []
    for arg in all_args:
        if arg.arg == "self":
            continue
        params.append(ParamDoc(
            name=arg.arg,
            annotation=_annotation_str(arg.annotation),
            default=defaults_map.get(arg.arg, ""),
        ))

    docstring = ast.get_docstring(node) or ""
    returns = _annotation_str(node.returns)
    is_private = node.name.startswith("_")

    return FunctionDoc(
        name=node.name,
        is_private=is_private,
        params=params,
        returns=returns,
        docstring=docstring,
        lineno=node.lineno,
    )


def _parse_class(node: ast.ClassDef) -> ClassDoc:
    docstring = ast.get_docstring(node) or ""
    methods: list[FunctionDoc] = []
    for item in node.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            methods.append(_parse_function(item))
    return ClassDoc(name=node.name, docstring=docstring, methods=methods, lineno=node.lineno)


def parse_module(path: Path) -> ModuleDoc:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    module_docstring = ast.get_docstring(tree) or ""

    functions: list[FunctionDoc] = []
    classes: list[ClassDoc] = []

    for node in ast.iter_child_nodes(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.append(_parse_function(node))
        elif isinstance(node, ast.ClassDef):
            classes.append(_parse_class(node))

    return ModuleDoc(
        name=path.stem,
        filepath=path,
        docstring=module_docstring,
        functions=functions,
        classes=classes,
    )


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

_CSS = """
:root {
  --accent: #3b6e8f;
  --bg: #f9fafb;
  --card: #ffffff;
  --border: #dde2e8;
  --code-bg: #f0f4f8;
  --private: #8a8a9a;
  --text: #1a1f2e;
  --small: #5a6070;
}

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  font-size: 14px;
  line-height: 1.6;
  color: var(--text);
  background: var(--bg);
  padding: 2rem 1rem;
}

.page-wrapper { max-width: 960px; margin: 0 auto; }

h1.page-title {
  font-size: 2rem;
  color: var(--accent);
  border-bottom: 3px solid var(--accent);
  padding-bottom: .5rem;
  margin-bottom: 2rem;
}

/* TOC */
.toc {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 1.2rem 1.5rem;
  margin-bottom: 2.5rem;
}
.toc h2 { font-size: 1rem; margin-bottom: .6rem; color: var(--accent); }
.toc ul { list-style: none; columns: 2; }
.toc li { margin-bottom: .25rem; }
.toc a { color: var(--accent); text-decoration: none; }
.toc a:hover { text-decoration: underline; }

/* Module section */
.module-section {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 8px;
  margin-bottom: 2.5rem;
  overflow: hidden;
}

.module-header {
  background: var(--accent);
  color: #fff;
  padding: .75rem 1.25rem;
  display: flex;
  align-items: baseline;
  gap: .75rem;
}
.module-header h2 { font-size: 1.15rem; }
.module-header .filepath { font-size: .8rem; opacity: .8; }

.module-docstring {
  padding: .9rem 1.25rem;
  color: var(--small);
  border-bottom: 1px solid var(--border);
  font-size: .85rem;
  white-space: pre-wrap;
}

/* Class block */
.class-block {
  border-top: 2px solid var(--border);
  padding: 1rem 1.25rem;
}
.class-block + .class-block { border-top-style: solid; }

.class-header {
  display: flex;
  align-items: baseline;
  gap: .5rem;
  margin-bottom: .4rem;
}
.class-header h3 { font-size: 1.05rem; color: var(--accent); }
.badge {
  font-size: .7rem;
  padding: .1rem .45rem;
  border-radius: 3px;
  font-weight: 600;
  text-transform: uppercase;
}
.badge-class  { background: #dbeafe; color: #1e4aac; }
.badge-private { background: #f3f4f6; color: var(--private); }

.class-docstring {
  font-size: .85rem;
  color: var(--small);
  margin-bottom: .8rem;
  white-space: pre-wrap;
}

/* Function / method cards */
.func-list { display: flex; flex-direction: column; gap: .6rem; }

.func-card {
  border: 1px solid var(--border);
  border-radius: 5px;
  overflow: hidden;
}
.func-card.private { opacity: .7; }

.func-signature {
  background: var(--code-bg);
  padding: .5rem .85rem;
  font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
  font-size: .82rem;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: .4rem;
  border-bottom: 1px solid var(--border);
}
.func-name   { color: #0f5ea8; font-weight: 700; }
.func-params { color: var(--text); }
.func-return { color: #6b3dab; margin-left: auto; white-space: nowrap; }

.func-body { padding: .6rem .85rem; }

.func-docstring {
  font-size: .83rem;
  color: var(--small);
  white-space: pre-wrap;
  margin-bottom: .4rem;
}

.param-table {
  width: 100%;
  border-collapse: collapse;
  font-size: .8rem;
  margin-top: .3rem;
}
.param-table th {
  text-align: left;
  font-weight: 600;
  color: var(--small);
  border-bottom: 1px solid var(--border);
  padding: .2rem .5rem;
}
.param-table td { padding: .2rem .5rem; vertical-align: top; }
.param-table tr:nth-child(even) td { background: var(--code-bg); }
.param-name { font-family: monospace; color: #0f5ea8; }
.param-type { font-family: monospace; color: #6b3dab; }
.param-default { font-family: monospace; color: #888; }

/* standalone functions */
.func-section {
  padding: .9rem 1.25rem;
  border-top: 1px solid var(--border);
}
.func-section-title {
  font-size: .75rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: .05em;
  color: var(--small);
  margin-bottom: .6rem;
}

@media print {
  body { background: white; padding: .5rem; }
  .module-section { break-inside: avoid; page-break-inside: avoid; }
  .func-card { break-inside: avoid; page-break-inside: avoid; }
  a { color: inherit; text-decoration: none; }
}
"""


def _fmt_signature(fn: FunctionDoc) -> str:
    params_html = []
    for p in fn.params:
        part = f'<span class="param-name">{escape(p.name)}</span>'
        if p.annotation:
            part += f': <span class="param-type">{escape(p.annotation)}</span>'
        if p.default:
            part += f' = <span class="param-default">{escape(p.default)}</span>'
        params_html.append(part)
    params_str = ", ".join(params_html)
    ret = f' → <span class="func-return">{escape(fn.returns)}</span>' if fn.returns else ""
    return (
        f'<span class="func-name">{escape(fn.name)}</span>'
        f'<span class="func-params">({params_str})</span>'
        f"{ret}"
    )


def _render_function_card(fn: FunctionDoc) -> str:
    private_cls = " private" if fn.is_private else ""
    docstring_html = (
        f'<div class="func-docstring">{escape(fn.docstring)}</div>'
        if fn.docstring else ""
    )

    param_rows = ""
    for p in fn.params:
        param_rows += (
            f"<tr>"
            f'<td class="param-name">{escape(p.name)}</td>'
            f'<td class="param-type">{escape(p.annotation)}</td>'
            f'<td class="param-default">{escape(p.default)}</td>'
            f"</tr>"
        )
    param_table = ""
    if param_rows:
        param_table = (
            '<table class="param-table">'
            "<thead><tr><th>Parameter</th><th>Type</th><th>Default</th></tr></thead>"
            f"<tbody>{param_rows}</tbody></table>"
        )

    return (
        f'<div class="func-card{private_cls}">'
        f'<div class="func-signature">{_fmt_signature(fn)}</div>'
        f'<div class="func-body">{docstring_html}{param_table}</div>'
        f"</div>"
    )


def _render_class(cls: ClassDoc) -> str:
    docstring_html = (
        f'<div class="class-docstring">{escape(cls.docstring)}</div>'
        if cls.docstring else ""
    )
    methods_html = '<div class="func-list">' + "".join(
        _render_function_card(m) for m in cls.methods
    ) + "</div>"

    return (
        f'<div class="class-block" id="{escape(cls.name)}">'
        f'<div class="class-header">'
        f'<h3>{escape(cls.name)}</h3>'
        f'<span class="badge badge-class">class</span>'
        f'</div>'
        f"{docstring_html}"
        f"{methods_html}"
        f"</div>"
    )


def _render_module(mod: ModuleDoc) -> str:
    try:
        rel = mod.filepath.relative_to(Path.cwd())
    except ValueError:
        rel = mod.filepath
    header = (
        f'<div class="module-header">'
        f'<h2 id="mod-{escape(mod.name)}">{escape(mod.name)}</h2>'
        f'<span class="filepath">{escape(str(rel))}</span>'
        f'</div>'
    )

    doc_html = (
        f'<div class="module-docstring">{escape(mod.docstring)}</div>'
        if mod.docstring else ""
    )

    # Public functions section
    pub_fns = [f for f in mod.functions if not f.is_private]
    priv_fns = [f for f in mod.functions if f.is_private]
    fn_sections = ""
    if pub_fns:
        cards = "".join(_render_function_card(f) for f in pub_fns)
        fn_sections += (
            f'<div class="func-section">'
            f'<div class="func-section-title">Functions</div>'
            f'<div class="func-list">{cards}</div>'
            f'</div>'
        )
    if priv_fns:
        cards = "".join(_render_function_card(f) for f in priv_fns)
        fn_sections += (
            f'<div class="func-section">'
            f'<div class="func-section-title">Internal helpers</div>'
            f'<div class="func-list">{cards}</div>'
            f'</div>'
        )

    classes_html = "".join(_render_class(c) for c in mod.classes)

    return (
        f'<div class="module-section">'
        f"{header}{doc_html}{fn_sections}{classes_html}"
        f"</div>"
    )


def _read_version() -> str:
    import re
    candidate = Path.cwd() / "pyproject.toml"
    if not candidate.exists():
        return ""
    try:
        import tomllib
        with open(candidate, "rb") as f:
            return tomllib.load(f).get("project", {}).get("version", "")
    except ImportError:
        pass
    m = re.search(r'^version\s*=\s*"([^"]+)"', candidate.read_text(), re.MULTILINE)
    return m.group(1) if m else ""


def render_html(modules: list[ModuleDoc], package_name: str) -> str:
    version = _read_version()
    version_label = f" {escape(version)}" if version else ""
    toc_items = "".join(
        f'<li><a href="#mod-{escape(m.name)}">{escape(m.name)}</a></li>'
        for m in modules
    )
    toc = (
        f'<nav class="toc"><h2>Modules</h2><ul>{toc_items}</ul></nav>'
    )

    body = "".join(_render_module(m) for m in modules)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(package_name)}{version_label} — API Reference</title>
<style>{_CSS}</style>
</head>
<body>
<div class="page-wrapper">
<h1 class="page-title">{escape(package_name)}{version_label} — API Reference</h1>
{toc}
{body}
</div>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--src-dir", default="src/kardemumma", help="Source directory to scan (default: src/kardemumma)")
    parser.add_argument("--out-dir", default="dist", help="Output directory (default: dist)")
    parser.add_argument("--package-name", default="kardemumma", help="Display name for the package")
    parser.add_argument("--no-pdf", action="store_true", help="Skip PDF generation")
    args = parser.parse_args()

    src_dir = Path(args.src_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not src_dir.exists():
        print(f"ERROR: source directory not found: {src_dir}", file=sys.stderr)
        sys.exit(1)

    py_files = sorted(f for f in src_dir.glob("*.py") if f.name != "__init__.py")
    print(f"Parsing {len(py_files)} module(s) in {src_dir} ...")

    modules: list[ModuleDoc] = []
    for path in py_files:
        print(f"  {path.name}")
        modules.append(parse_module(path))

    html_content = render_html(modules, args.package_name)

    html_path = out_dir / "kardemumma_api.html"
    html_path.write_text(html_content, encoding="utf-8")
    print(f"\nHTML written → {html_path}")

    if not args.no_pdf:
        pdf_path = out_dir / "kardemumma_api.pdf"
        chrome_bins = ["google-chrome", "chromium-browser", "chromium"]
        chrome = None
        for bin_name in chrome_bins:
            result = subprocess.run(["which", bin_name], capture_output=True, text=True)
            if result.returncode == 0:
                chrome = result.stdout.strip()
                break

        if chrome is None:
            print("WARNING: No Chrome/Chromium binary found — skipping PDF. Install google-chrome or run with --no-pdf.")
        else:
            print(f"Generating PDF via {chrome} ...")
            cmd = [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--disable-dev-shm-usage",
                f"--print-to-pdf={pdf_path.resolve()}",
                "--print-to-pdf-no-header",
                str(html_path.resolve()),
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                print(f"PDF  written → {pdf_path}")
            else:
                print(f"WARNING: Chrome exited with code {result.returncode}")
                print(result.stderr[:500])


if __name__ == "__main__":
    main()
