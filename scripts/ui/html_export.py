import os
import html
import webbrowser
import logging
from datetime import datetime
from pathlib import Path

from scripts.core.session_tmp import session_path

try:
    from scripts.core.i18n import tr
except ImportError:
    try:
        from core.i18n import tr
    except ImportError:
        try:
            from i18n import tr
        except ImportError:
            def tr(text: str) -> str:
                return text

from scripts.core.comments import LOW_COVERAGE, result_comments

try:
    import PySimpleGUI as sg
except ImportError:
    sg = None


def provenance_text(provenance):
    """Ligne de provenance affichée en pied de page de l'export."""
    p = provenance or {}
    parts = [f"TGV {p.get('tgv_version', '?')} ({p.get('tgv_commit', '?')})"]
    if p.get("trgt_version"):
        parts.append(f"TRGT {p['trgt_version']}")
    if p.get("catalog"):
        parts.append(f"catalogue {p['catalog']}")
    if p.get("thresholds_sha256"):
        parts.append(f"clinical_thresholds.yaml SHA-256 {p['thresholds_sha256'][:16]}")
    return " · ".join(parts)


def generate_html_table(headers, rows, sample_name, run_id=None, low_depth_threshold=None, provenance=None):
    logging.info(f"Generating export HTML table for patient '{sample_name}' with {len(rows)} selected rows.")

    headers = list(headers)
    comments_label = tr("Commentaires")
    if comments_label not in headers:
        headers.append(comments_label)

    thead_html = (
        "<tr><th class='col-checkbox'></th>"
        + "".join(f"<th>{tr(h)}</th>" for h in headers)
        + "</tr>"
    )

    tbody_rows = []
    for r_dict in rows:
        row_html = "<tr>"
        row_html += "<td class='col-checkbox'><input type='checkbox' class='row-checkbox' checked></td>"
        
        r_obj = r_dict.get("Result_obj")
        comments = result_comments(r_obj, low_depth_threshold) if r_obj else []

        is_classif_modified = False
        if r_obj:
            if getattr(r_obj, "classification1_bio", None) or getattr(r_obj, "classification2_bio", None):
                is_classif_modified = True
        if r_dict.get("Classification") != r_dict.get("Classification_auto"):
            is_classif_modified = True

        is_genotype_modified = False
        if r_obj:
            if getattr(r_obj, "genotype1_bio", None) is not None or getattr(r_obj, "genotype2_bio", None) is not None:
                is_genotype_modified = True
        gt_current = r_dict.get("Génotype") or r_dict.get("Genotype", "")
        if gt_current != r_dict.get("Genotype_auto"):
            is_genotype_modified = True

        for h in headers:
            if h in ("Commentaires", comments_label):
                # Seul le libellé fixe est traduit : les notes sont des données copiées vers le SIL
                comment_val = "; ".join(tr(c) if c == LOW_COVERAGE else c for c in comments)
                row_html += f"<td><span contenteditable='true' class='comment-input' data-placeholder=\"{tr('Ajouter un commentaire...')}\">{comment_val}</span></td>"
            else:
                val = r_dict.get(h, '')
                classes = []
                if h == "Classification" and is_classif_modified:
                    classes.append("modified-cell")
                elif h in ("Génotype", "Genotype") and is_genotype_modified:
                    classes.append("modified-cell")
                
                class_attr = f" class='{' '.join(classes)}'" if classes else ""
                row_html += f"<td{class_attr}>{tr(str(val))}</td>"

        row_html += "</tr>"
        tbody_rows.append(row_html)
    tbody_html = "".join(tbody_rows)

    page = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{tr("Résultats TGV")} - {sample_name}</title>
<style>
  :root {{
    --bg:        #f7f7f8;
    --surface:   #ffffff;
    --surface2:  #fdf4f7;
    --border:    #e8dde2;
    --accent:    #c8336a;
    --accent-lt: #fdf0f4;
    --accent2:   #9e1f4f;
    --muted:     #8a7480;
    --text:      #1a1015;
    --text-dim:  #5a4050;
    --success:   #2a7a4a;
    --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    --font-mono: SFMono-Regular, Consolas, "Liberation Mono", Menlo, Monaco, Courier, monospace;
  }}

  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    background: var(--bg);
    color: var(--text);
    font-family: var(--font-sans);
    font-size: 14px;
    line-height: 1.6;
    min-height: 100vh;
  }}

  .header {{
    background: var(--surface);
    border-bottom: 3px solid var(--accent);
    padding: 18px 40px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 24px;
  }}
  .header-brand {{
    display: flex;
    align-items: center;
    gap: 16px;
  }}
  .header-logo-mark {{
    width: 36px;
    height: 36px;
    border-radius: 50%;
    background: linear-gradient(135deg, #e8457a 0%, #9e1f4f 100%);
    flex-shrink: 0;
  }}
  .header-titles h1 {{
    font-size: 16px;
    font-weight: 600;
    color: var(--text);
    letter-spacing: -0.2px;
    line-height: 1.2;
  }}
  .header-titles .subtitle {{
    font-family: var(--font-mono);
    font-size: 11px;
    color: var(--muted);
    margin-top: 2px;
  }}

  .main {{ padding: 32px 40px; max-width: 1600px; margin: 0 auto; }}
  .section {{ margin-bottom: 44px; }}

  .section-header {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 16px;
    padding-bottom: 10px;
    border-bottom: 1px solid var(--border);
  }}
  .section-title {{
    font-size: 13px;
    font-weight: 600;
    color: var(--text-dim);
    text-transform: uppercase;
    letter-spacing: 0.6px;
  }}

  .table-wrapper {{
    overflow-x: auto;
    border: 1px solid var(--border);
    border-radius: 6px;
    background: var(--surface);
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12.5px;
    white-space: nowrap;
  }}
  thead th {{
    background: var(--surface2);
    color: var(--text-dim);
    font-weight: 600;
    font-size: 10.5px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    padding: 9px 14px;
    text-align: left;
    border-bottom: 2px solid var(--border);
  }}
  tbody tr {{ border-bottom: 1px solid var(--border); }}
  tbody tr:last-child {{ border-bottom: none; }}
  tbody tr:hover {{ background: var(--surface2); }}
  tbody td {{
    padding: 8px 14px;
    color: var(--text);
    font-family: var(--font-mono);
    font-size: 12px;
  }}
  
  tbody td:nth-child(2) {{
    color: var(--accent2);
    font-weight: 500;
    font-family: var(--font-sans);
    font-size: 12.5px;
  }}

  td.modified-cell {{
    position: relative;
    background-color: var(--accent-lt) !important;
    border-bottom: 1.5px dashed var(--accent) !important;
  }}
  td.modified-cell::after {{
    content: " ✎";
    color: var(--accent);
    font-weight: bold;
    font-size: 11px;
    font-family: var(--font-sans);
    pointer-events: none;
  }}

  .comment-input {{
    display: block;
    width: 100%;
    min-width: 220px;
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 4px 8px;
    font-family: var(--font-sans);
    font-size: 12px;
    background: var(--surface);
    outline: none;
  }}
  .comment-input:focus {{
    border-color: var(--accent);
    background: #fff;
  }}
  .comment-input:empty::before {{
    content: attr(data-placeholder);
    color: var(--muted);
    font-style: italic;
  }}

  .col-checkbox {{ width: 45px; text-align: center; padding: 8px; }}
  .row-checkbox {{ transform: scale(1.15); cursor: pointer; accent-color: var(--accent); }}

  .btn-container {{
    margin-top: 20px;
    display: flex;
    gap: 10px;
    flex-wrap: wrap;
  }}
  
  .btn {{
    display: inline-flex;
    align-items: center;
    gap: 6px;
    padding: 8px 16px;
    border-radius: 4px;
    font-size: 11px;
    font-family: var(--font-mono);
    font-weight: 500;
    cursor: pointer;
    text-transform: uppercase;
  }}
  
  .btn-primary {{ background: var(--accent-lt); border: 1px solid #d9a0b4; color: var(--accent2); }}
  .btn-primary:hover {{ background: var(--accent); border-color: var(--accent); color: #fff; }}
  .btn-primary.copied {{ background: var(--success); border-color: var(--success); color: #fff; }}
  .btn-secondary {{ background: var(--surface); border: 1px solid var(--border); color: var(--text-dim); }}
  .btn-secondary:hover {{ background: var(--surface2); color: var(--accent2); }}

  .footer {{
    margin-top: 56px;
    padding: 18px 40px;
    border-top: 1px solid var(--border);
    display: flex;
    align-items: center;
    justify-content: space-between;
  }}
  .provenance {{
    font-family: var(--font-mono);
    font-size: 11px;
    color: var(--muted);
  }}
  .footer-dot {{
    width: 8px; height: 8px;
    border-radius: 50%;
    background: linear-gradient(135deg, #e8457a, #9e1f4f);
    display: inline-block;
    margin-right: 8px;
  }}
</style>
</head>
<body>

<div class="header">
  <div class="header-brand">
    <div class="header-logo-mark"></div>
    <div class="header-titles">
      <h1>TRGT Global Viewer &mdash; {tr("Export Data")}</h1>
      <div class="subtitle">{sample_name}</div>
    </div>
  </div>
</div>

<main class="main">
  <section class="section">
    <div class="section-header">
      <div class="section-title">{tr("Résultats TGV")}</div>
    </div>
    
    <div class="table-wrapper">
      <table>
        <thead>
          {thead_html}
        </thead>
        <tbody>
          {tbody_html}
        </tbody>
      </table>
    </div>

    <div class="btn-container">
      <button type="button" class="btn btn-secondary" onclick="setAllCheckboxes(true)">{tr("Tout cocher")}</button>
      <button type="button" class="btn btn-secondary" onclick="setAllCheckboxes(false)">{tr("Tout décocher")}</button>
      <button type="button" id="btn-copy" class="btn btn-primary" onclick="copySelectedRows()">
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" width="12" height="12">
          <rect x="9" y="9" width="13" height="13" rx="2"/>
          <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/>
        </svg>
        <span>{tr("Copier la sélection")}</span>
      </button>
    </div>
  </section>
</main>

<footer class="footer">
  <div><span class="footer-dot"></span><span>TRGT Global Viewer &bull; {sample_name}</span></div>
  <div class="provenance">{html.escape(provenance_text(provenance))}</div>
</footer>

<script>
  function setAllCheckboxes(checked) {{
    document.querySelectorAll('.row-checkbox').forEach(cb => cb.checked = checked);
  }}

  async function copySelectedRows() {{
    const runId = "{run_id}";
    const sampleName = "{sample_name}";
    const rows = document.querySelectorAll('table tbody tr');
    let tsvContent = '';
    let hasSelected = false;

    rows.forEach(row => {{
      const checkbox = row.querySelector('.row-checkbox');
      if (checkbox && checkbox.checked) {{
        hasSelected = true;
        const cells = Array.from(row.querySelectorAll('td')).slice(1).map(td => {{
          const tempTd = td.cloneNode(true);
          tempTd.classList.remove('modified-cell');
          const commentInput = tempTd.querySelector('.comment-input');
          let textVal = commentInput ? commentInput.innerText : tempTd.innerText;
          textVal = textVal.replace(/[\\u26A0\\uFE0F]/g, '');
          return textVal.trim();
        }});
        const rowData = [runId, sampleName, ...cells];
        tsvContent += rowData.join('\\t') + '\\n';
      }}
    }});

    if (!hasSelected) {{
      alert("{tr("Aucune ligne n'est sélectionnée pour la copie.")}");
      return;
    }}

    try {{
      await navigator.clipboard.writeText(tsvContent);
      const btn = document.getElementById('btn-copy');
      btn.classList.add('copied');
      const textSpan = btn.querySelector('span');
      const originalText = textSpan.innerText;
      textSpan.innerText = "{tr("Copié !")}";
      setTimeout(() => {{
        btn.classList.remove('copied');
        textSpan.innerText = originalText;
      }}, 2000);
    }} catch (err) {{
      alert("{tr("Erreur lors de la copie : ")}" + err);
    }}
  }}
</script>
</body>
</html>
"""
    return page


def save_and_open_html(html_content, name="tgv_export"):
    try:
        # Répertoire de session (supprimé à la fermeture), nom unique par export
        tmp_path = session_path("exports", f"{os.path.basename(name)}_{datetime.now():%Y%m%d_%H%M%S_%f}.html")
        with open(tmp_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        uri = Path(tmp_path).as_uri()
        webbrowser.open(uri)
    except Exception as e:
        if sg:
            sg.popup_error(f"{tr('Erreur')} :\n{e}")
        else:
            logging.error(f"Error opening HTML: {e}")