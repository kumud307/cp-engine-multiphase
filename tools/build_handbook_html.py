"""Build docs/TEAM_HANDBOOK.html (a self-contained web page) from docs/TEAM_HANDBOOK.md.

Embeds the key validation figures as data URIs so the page works offline and can be shared.
Run from the project root:  python tools/build_handbook_html.py   (needs: pip install markdown)
"""
import base64
import re
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent.parent
md_text = (ROOT / "docs" / "TEAM_HANDBOOK.md").read_text(encoding="utf-8")

FIGS = [
    ("fig_HE_water_ethanol", "Excess enthalpy of water + ethanol: the refit (solid) follows the sign change that VLE parameters (dashed) miss."),
    ("fig_cpE_from_HE_water_ethanol", "Cpᴱ of water + ethanol from measured Hᴱ(T), against the two NRTL parameter sets."),
    ("fig_cp_methanol_water", "Independent test: Cp of water + methanol (Tanaka et al. 2007) vs ideal mixture and Hᴱ-fitted NRTL."),
    ("fig_nanofluid_parity", "69 nanofluid DSC points: Model II (= Route A) on the diagonal, Model I far above it."),
    ("fig_freezing_points", "Freezing points of aqueous ethanol: the engine's ideal + ΔCp liquidus against Lange's data."),
    ("fig_cpT_ice_ethanol_10wt", "Apparent Cp of a 10 wt% ethanol ice slurry, split into ideal, excess and latent parts."),
]
fig_html = ['<section class="figs"><h2 id="key-figures">Key validation figures</h2><div class="fig-grid">']
for name, cap in FIGS:
    b64 = base64.b64encode((ROOT / "figures" / f"{name}.png").read_bytes()).decode()
    fig_html.append(f'<figure><img src="data:image/png;base64,{b64}" alt="{cap}"><figcaption>{cap}</figcaption></figure>')
fig_html.append("</div></section>")

body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "toc", "sane_lists"])
# insert the figures right before section 7
body = re.sub(r'(<h2 id="7-[^"]*">)', "\n".join(fig_html) + r"\n\1", body, count=1)
toc = "".join(f'<li><a href="#{m.group(1)}">{re.sub("<[^>]+>", "", m.group(2))}</a></li>'
              for m in re.finditer(r'<h2 id="([^"]+)">(.*?)</h2>', body))

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cp Engine Team Handbook</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap" rel="stylesheet">
<style>
:root {{
  /* Layout: one reading column with a sticky contents rail on wide screens; a lab-notebook feel. */
  --paper: #f7f8f6; --ink: #1b1f1e; --muted: #5b6461; --rule: #d9dedb; --accent: #1f5f8b; --tint: #eaf0f3;
  --font-display: "IBM Plex Serif", Georgia, serif;
  --font-body: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  --font-mono: "IBM Plex Mono", Consolas, monospace;
  color-scheme: light;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper: #131716; --ink: #e6eae8; --muted: #a2aba7; --rule: #2f3634; --accent: #7fb6dc; --tint: #1c2427;
    color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --paper: #131716; --ink: #e6eae8; --muted: #a2aba7; --rule: #2f3634; --accent: #7fb6dc; --tint: #1c2427;
  color-scheme: dark;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--paper); color: var(--ink); font: 16px/1.6 var(--font-body); }}
.wrap {{ display: grid; grid-template-columns: minmax(0, 1fr); gap: 32px; max-width: 1180px; margin: 0 auto; padding: 32px 20px 80px; }}
@media (min-width: 1000px) {{ .wrap {{ grid-template-columns: 240px minmax(0, 1fr); }} }}
nav {{ font-size: 14px; }}
@media (min-width: 1000px) {{ nav {{ position: sticky; top: 24px; align-self: start; max-height: calc(100vh - 48px); overflow-y: auto; }} }}
nav p {{ margin: 0 0 8px; font-weight: 600; color: var(--muted); text-transform: uppercase; letter-spacing: .06em; font-size: 12px; }}
nav ol {{ margin: 0; padding-left: 18px; display: grid; gap: 4px; }}
nav a {{ color: var(--ink); text-decoration: none; }}
nav a:hover, nav a:focus-visible {{ color: var(--accent); text-decoration: underline; }}
main {{ min-width: 0; }}
main > * {{ max-width: 72ch; }}
h1, h2, h3 {{ font-family: var(--font-display); line-height: 1.25; text-wrap: balance; }}
h1 {{ font-size: 2.1rem; margin: 0 0 8px; }}
h2 {{ font-size: 1.5rem; margin: 48px 0 12px; padding-top: 12px; border-top: 2px solid var(--accent); }}
h3 {{ font-size: 1.15rem; margin: 28px 0 8px; }}
a {{ color: var(--accent); }}
a:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
code {{ font-family: var(--font-mono); font-size: .88em; background: var(--tint); padding: 1px 5px; border-radius: 4px; }}
pre {{ background: var(--tint); border: 1px solid var(--rule); border-radius: 6px; padding: 14px 16px; overflow-x: auto; max-width: 100%; }}
pre code {{ background: none; padding: 0; }}
blockquote {{ margin: 16px 0; padding: 10px 18px; border-left: 4px solid var(--accent); background: var(--tint); }}
hr {{ border: 0; border-top: 1px solid var(--rule); margin: 28px 0; }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; font-variant-numeric: tabular-nums; }}
.tbl {{ overflow-x: auto; max-width: 100%; margin: 12px 0 20px; }}
th, td {{ border-bottom: 1px solid var(--rule); padding: 7px 10px; text-align: left; vertical-align: top; }}
th {{ font-weight: 600; border-bottom: 2px solid var(--ink); }}
main > ol > li, main > ul > li {{ margin: 4px 0; }}
.figs {{ max-width: none; }}
.fig-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }}
figure {{ margin: 0; padding: 12px; background: #ffffff; border: 1px solid var(--rule); border-radius: 6px; }}
figure img {{ max-width: 100%; height: auto; display: block; }}
figcaption {{ font-size: 13px; color: #3c4442; margin-top: 8px; }}
.meta {{ color: var(--muted); margin: 0 0 24px; }}
</style>
</head>
<body>
<div class="wrap">
<nav aria-label="Contents"><p>Contents</p><ol>{toc}</ol></nav>
<main>
<p class="meta">Generated from docs/TEAM_HANDBOOK.md in the project repository · CL 202 Topic 5 · demo on 20 Nov 2026</p>
{body}
</main>
</div>
</body>
</html>
"""
html = re.sub(r"(<table>.*?</table>)", r'<div class="tbl">\1</div>', html, flags=re.S)
out = ROOT / "docs" / "TEAM_HANDBOOK.html"
out.write_text(html, encoding="utf-8")
print("written", out, f"{len(html) / 1e6:.2f} MB")
