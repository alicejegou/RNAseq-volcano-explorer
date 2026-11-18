# Volcano explorer: overlay several DESeq2 results and link genes
# Upload a base DESeq2 table and overlay tables, then add gene lists.
# Runs in the browser with Shinylive: no server, data never leaves the
# user's computer.

library(shiny)
library(ggplot2)
library(dplyr)
library(scales)
library(plotly)

# HELPERS ______________________________________________________________

# Default colors of the experiments: base first, then one per overlay
default_cols <- c("#41ba43", "#c20845", "#1f77b4", "#ff7f0e",
                  "#9467bd", "#8c564b", "#e377c2", "#17becf")
default_col <- function(i) {
  default_cols[(i - 1) %% length(default_cols) + 1]
}

# Default colors of the gene lists (rings, links, labels)
layer_cols <- c("#111111", "#e6194b", "#4363d8", "#f58231",
                "#911eb4", "#42d4f4", "#f032e6", "#9a6324")

# What a gene list can contain
src_choices <- c(
  "Gene list" = "list",
  "DEG common to all experiments" = "common",
  "Top N DEG of base, by padj" = "top_p",
  "Top N DEG of base, by |log2FC|" = "top_fc"
)

`%||%` <- function(a, b) if (is.null(a)) b else a

# Safe input id built from an experiment name
exp_id <- function(nm) gsub("[^A-Za-z0-9_]", "_", nm)

# Numeric input value, or a default when empty (NA / NULL)
num_or <- function(v, default) {
  if (is.null(v) || is.na(v)) default else v
}

# Read a csv/tsv/;-separated table. The separator is guessed from the
# header line; row names (write.table style) become a column "X".
read_table_auto <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  seps  <- c(",", "\t", ";")
  n_sep <- vapply(seps, function(s) {
    lengths(regmatches(first, gregexpr(s, first, fixed = TRUE)))
  }, integer(1))
  df <- read.table(
    path, sep = seps[which.max(n_sep)], header = TRUE, quote = "\"",
    check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""
  )
  if (.row_names_info(df) > 0) {
    df <- cbind(X = rownames(df), df)
    rownames(df) <- NULL
  }
  if (names(df)[1] == "") names(df)[1] <- "X"
  df
}

# First column name matching a pattern, else a default
guess_col <- function(cols, pattern, default) {
  hit <- grep(pattern, cols, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) hit[1] else default
}

# Guess the gene ID, log2FC and padj columns of one table
guess_cols <- function(cols) {
  list(
    id   = guess_col(cols, "^(x|gene|gene_?id|id|name|symbol|locus)$",
                     cols[1]),
    lfc  = guess_col(cols, "log2_?fold|log2_?fc|^logfc$|^lfc$",
                     cols[min(2, length(cols))]),
    padj = guess_col(cols, "^padj$|^adj|fdr|q_?val", cols[length(cols)])
  )
}

# Keep only gene / log2FC / padj under standard names
std_table <- function(df, c_id, c_lfc, c_padj) {
  df %>%
    transmute(
      gene = as.character(.data[[c_id]]),
      lfc  = suppressWarnings(as.numeric(.data[[c_lfc]])),
      padj = suppressWarnings(as.numeric(.data[[c_padj]]))
    ) %>%
    distinct(gene, .keep_all = TRUE) %>%
    mutate(nlp = -log10(padj))
}

# Gene list from text (space, comma, ; or newline separated)
read_gene_list <- function(txt = "") {
  g <- unlist(strsplit(txt, "[[:space:],;]+"))
  unique(g[g != ""])
}

# TRUE only for significant genes (NA -> FALSE)
is_deg <- function(padj, lfc, p, l) (padj < p & abs(lfc) > l) %in% TRUE

# Format one table column as text for display
fmt_cell <- function(v) {
  if (is.numeric(v)) v <- signif(v, 3)
  out <- as.character(v)
  out[is.na(v)] <- ""
  out
}

# HTML table whose rows carry the gene and its plot coordinates (JSON),
# so that hovering a row draws a ring around the gene on the plot.
# Headers are clickable (sort) and carry a filter box (client side JS).
gene_table <- function(r, pts, main = FALSE) {
  types <- vapply(r, function(v) {
    if (is.logical(v)) "lgl" else if (is.numeric(v)) "num" else "chr"
  }, character(1))
  cells <- lapply(r, fmt_cell)
  body <- lapply(seq_len(nrow(r)), function(i) {
    p <- pts[[r$gene[i]]]
    tags$tr(
      `data-gene` = r$gene[i],
      `data-pts`  = if (is.null(p)) "[]" else as.character(p),
      unname(lapply(cells, function(col) tags$td(col[i])))
    )
  })
  heads <- lapply(seq_along(r), function(j) {
    tags$th(
      `data-type` = types[j],
      tags$div(
        class = "sort-lab",
        title = "Click: sort. Shift+click: add a sort level",
        names(r)[j], tags$span(class = "sort-ind")
      ),
      if (types[j] == "lgl") {
        tags$select(
          class = "flt-in",
          tags$option(value = "", "all"),
          tags$option(value = "TRUE", "TRUE"),
          tags$option(value = "FALSE", "FALSE"),
          tags$option(value = "NA", "NA")
        )
      } else {
        tags$input(
          type = "text", class = "flt-in",
          placeholder = if (types[j] == "num") ">1 <5" else "contains..."
        )
      }
    )
  })
  tags$div(
    class = if (main) "tbl-wrap main" else "tbl-wrap",
    tags$div(
      class = "tbl-count",
      tags$span(class = "tbl-n", paste0(nrow(r), " rows")), " ",
      tags$a(href = "#", class = "tbl-reset", "reset sort / filters")
    ),
    tags$div(
      class = "tbl-box",
      style = "max-height: 320px; overflow: auto; position: relative;",
      tags$table(
        class = "table table-condensed table-striped gene-table",
        tags$thead(tags$tr(heads)),
        tags$tbody(body)
      )
    ),
    if (main) tags$button(
      type = "button", class = "btn btn-default btn-sm",
      style = "margin-top: 6px;", onclick = "tblDownload(this);",
      "Download filtered table (.csv)"
    ),
    # A rebuilt table has no filters: tell the plot
    if (main) tags$script("if (window.tblSync) window.tblSync();")
  )
}

# One gene list block (name, color, display switches, source, genes).
# Inputs are named name_<id>, lcol_<id>, show_<id>, links_<id>,
# names_<id>, src_<id>, genes_<id>, lfile_<id>, n_<id>.
layer_ui <- function(id, name, color, text = "") {
  on_rm <- sprintf(
    "Shiny.setInputValue('rm_layer', '%s', {priority: 'event'});", id)
  on_col <- sprintf("Shiny.setInputValue('lcol_%s', this.value)", id)
  tags$div(
    id = paste0("layer_", id),
    style = paste0("border: 1px solid #ccc; border-radius: 4px;",
                   " padding: 6px; margin-bottom: 6px;"),
    fluidRow(
      column(7, textInput(paste0("name_", id), NULL, name, width = "100%")),
      column(3, tags$input(type = "color", value = color,
                           onchange = on_col)),
      column(2, tags$button(type = "button",
                            class = "btn btn-default btn-sm",
                            onclick = on_rm, "x"))
    ),
    fluidRow(
      column(4, checkboxInput(paste0("show_", id),  "Highlight", TRUE)),
      column(4, checkboxInput(paste0("links_", id), "Links", TRUE)),
      column(4, checkboxInput(paste0("names_", id), "Names", FALSE))
    ),
    selectInput(paste0("src_", id), NULL, src_choices, width = "100%"),
    conditionalPanel(
      sprintf("input.src_%s == 'list'", id),
      textAreaInput(paste0("genes_", id), NULL, text, rows = 4,
                    width = "100%",
                    placeholder = "Mp1g02140\nMp4g13300"),
      fileInput(paste0("lfile_", id), "...or upload a list",
                width = "100%")
    ),
    conditionalPanel(
      sprintf("input.src_%s == 'top_p' || input.src_%s == 'top_fc'",
              id, id),
      numericInput(paste0("n_", id), "N", 20, min = 1)
    )
  )
}

# FRONT-END CODE _______________________________________________________

app_css <- "
.gene-table th { position: sticky; top: 0; background: #fff; z-index: 2;
                 min-width: 90px; }
.gene-table tr[data-gene]:hover td,
.gene-table tr.hl-row td { background-color: #fff3b0 !important; }
.sort-lab { cursor: pointer; user-select: none; white-space: nowrap; }
.flt-in { width: 100%; font-weight: normal; font-size: 11px; }
.tbl-count { font-size: 12px; color: #555; margin: 4px 0; }
"

# Browser-side code: table hover rings, selection reset, downloads,
# table sort / filter
app_js <- r"---(
function volcanoUnHL() {
  var gd = document.getElementById('volcano');
  if (!gd || !gd.data) return;
  var idx = [];
  gd.data.forEach(function(t, i) {
    if (t.name === '__hl') idx.push(i);
  });
  if (idx.length > 0) Plotly.deleteTraces(gd, idx);
}

function volcanoHL(pts) {
  var gd = document.getElementById('volcano');
  if (!gd || !gd.data || pts.length === 0) return;
  volcanoUnHL();
  // WebGL trace added last: drawn on top of every other trace
  Plotly.addTraces(gd, {
    x: pts.map(function(p) { return p[0]; }),
    y: pts.map(function(p) { return p[1]; }),
    type: 'scattergl', mode: 'markers', name: '__hl',
    hoverinfo: 'skip', showlegend: false,
    marker: {
      symbol: 'circle-open', size: 18, color: 'black',
      line: {color: 'black', width: 3}
    }
  });
}

function volcanoClearSel() {
  var gd = document.getElementById('volcano');
  if (gd && gd.data) {
    Plotly.restyle(gd, {selectedpoints: [null]});
  }
}

// Downloads are done in the browser: no downloadHandler, no service
// worker needed (works on static hosting such as GitHub Pages)
function saveText(filename, text) {
  var blob = new Blob([text], {type: 'text/csv;charset=utf-8'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(function() { URL.revokeObjectURL(a.href); }, 60000);
}

(function waitShiny() {
  if (window.Shiny && Shiny.addCustomMessageHandler) {
    Shiny.addCustomMessageHandler('save_text', function(m) {
      saveText(m.filename, m.text);
    });
  } else {
    setTimeout(waitShiny, 100);
  }
})();

function volcanoPNG() {
  var gd = document.getElementById('volcano');
  if (gd && gd.data) {
    Plotly.downloadImage(gd, {
      format: 'png', filename: 'volcano',
      width: 1400, height: 1000, scale: 2
    });
  }
}

function sameRow(tr, e) {
  return e.relatedTarget && tr.contains(e.relatedTarget);
}

document.addEventListener('mouseover', function(e) {
  var tr = e.target.closest ? e.target.closest('tr[data-gene]') : null;
  if (tr && !sameRow(tr, e)) {
    volcanoHL(JSON.parse(tr.getAttribute('data-pts')));
  }
});

document.addEventListener('mouseout', function(e) {
  var tr = e.target.closest ? e.target.closest('tr[data-gene]') : null;
  if (tr && !sameRow(tr, e)) volcanoUnHL();
});

// TABLE SORT / FILTER (client side, per table)
function tblNumMatch(s, expr) {
  var x = s === '' ? NaN : parseFloat(s);
  if (isNaN(x)) return false;
  return expr.split(/\s+/).every(function(tok) {
    var m = tok.match(/^(abs)?(>=|<=|>|<|=)?(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)$/i);
    if (!m) return true;
    var v = parseFloat(m[3]), y = m[1] ? Math.abs(x) : x;
    switch (m[2] || '=') {
      case '>':  return y > v;
      case '>=': return y >= v;
      case '<':  return y < v;
      case '<=': return y <= v;
      default:   return y === v;
    }
  });
}

function tblMatch(s, f) {
  if (f.t === 'lgl') return (s === '' ? 'NA' : s) === f.v;
  if (f.t === 'num') return tblNumMatch(s, f.v);
  return s.toLowerCase().indexOf(f.v.toLowerCase()) >= 0;
}

function tblApply(tbl) {
  var tb = tbl.tBodies[0];
  var rows = Array.prototype.slice.call(tb.rows);
  if (rows.length && rows[0]._i === undefined) {
    rows.forEach(function(tr, i) { tr._i = i; });
  }
  var ths = Array.prototype.slice.call(tbl.querySelectorAll('thead th'));
  var filters = [];
  ths.forEach(function(th, c) {
    var inp = th.querySelector('.flt-in');
    var v = inp ? inp.value.trim() : '';
    if (v !== '') filters.push({c: c, t: th.getAttribute('data-type'), v: v});
  });
  var shown = 0;
  rows.forEach(function(tr) {
    var ok = filters.every(function(f) {
      return tblMatch(tr.children[f.c].textContent.trim(), f);
    });
    tr.style.display = ok ? '' : 'none';
    if (ok) shown++;
  });
  var sorts = tbl._sort || [];
  rows.sort(function(a, b) {
    for (var k = 0; k < sorts.length; k++) {
      var s = sorts[k];
      var va = a.children[s.c].textContent.trim();
      var vb = b.children[s.c].textContent.trim();
      if (va === '' || vb === '') {
        if (va === vb) continue;
        return va === '' ? 1 : -1;   // empty cells always last
      }
      var d = s.t === 'num' ? parseFloat(va) - parseFloat(vb)
                            : va.localeCompare(vb);
      if (d !== 0) return s.dir * d;
    }
    return a._i - b._i;
  });
  rows.forEach(function(tr) { tb.appendChild(tr); });
  ths.forEach(function(th, c) {
    var ind = th.querySelector('.sort-ind');
    var idx = -1;
    sorts.forEach(function(s, k) { if (s.c === c) idx = k; });
    ind.textContent = idx < 0 ? '' :
      ' ' + (sorts[idx].dir === 1 ? '\u25B2' : '\u25BC') +
      (sorts.length > 1 ? (idx + 1) : '');
  });
  var n = tbl.closest('.tbl-wrap').querySelector('.tbl-n');
  if (n) n.textContent = shown + ' / ' + rows.length + ' rows';
  if (tbl.closest('.tbl-wrap').classList.contains('main')) tblSync();
}

document.addEventListener('click', function(e) {
  var rst = e.target.closest ? e.target.closest('.tbl-reset') : null;
  if (rst) {
    e.preventDefault();
    var t = rst.closest('.tbl-wrap').querySelector('table');
    t._sort = [];
    t.querySelectorAll('.flt-in').forEach(function(i) { i.value = ''; });
    tblApply(t);
    return;
  }
  var lab = e.target.closest ? e.target.closest('.sort-lab') : null;
  if (!lab) return;
  var th = lab.parentNode, tbl = th.closest('table');
  var c = Array.prototype.indexOf.call(th.parentNode.children, th);
  var s = tbl._sort || [];
  var cur = s.filter(function(x) { return x.c === c; })[0];
  var next = !cur ? 1 : (cur.dir === 1 ? -1 : 0);   // asc -> desc -> off
  if (e.shiftKey) {
    if (cur && next === 0) s = s.filter(function(x) { return x.c !== c; });
    else if (cur) cur.dir = next;
    else s.push({c: c, dir: 1, t: th.getAttribute('data-type')});
  } else {
    s = next === 0 ? [] :
      [{c: c, dir: next, t: th.getAttribute('data-type')}];
  }
  tbl._sort = s;
  tblApply(tbl);
});

document.addEventListener('input', function(e) {
  if (e.target.classList && e.target.classList.contains('flt-in')) {
    tblApply(e.target.closest('table'));
  }
});

// MAIN TABLE FILTERS -> PLOT
// When a filter is active, the plot only shows the listed genes that
// pass it (the genes are sent to R, which redraws the plot).
window._fltSent = 'off';

function tblSync() {
  var cb = document.getElementById('link_flt');
  var tbl = document.querySelector('.tbl-wrap.main table');
  var on = false, genes = [];
  if (tbl && cb && cb.checked) {
    on = Array.prototype.some.call(
      tbl.querySelectorAll('.flt-in'),
      function(i) { return i.value.trim() !== ''; }
    );
    if (on) {
      Array.prototype.forEach.call(tbl.tBodies[0].rows, function(tr) {
        if (tr.style.display !== 'none') {
          genes.push(tr.getAttribute('data-gene'));
        }
      });
    }
  }
  var key = on ? genes.slice().sort().join('\u0001') : 'off';
  clearTimeout(window._fltT);
  if (key === window._fltSent) return;
  window._fltT = setTimeout(function() {
    window._fltSent = key;
    if (window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue('tbl_vis', {on: on, genes: genes});
    }
  }, on ? 400 : 0);
}

document.addEventListener('change', function(e) {
  if (e.target && e.target.id === 'link_flt') tblSync();
});

// Send the visible rows (in displayed order) to R for the csv download
function tblDownload(btn) {
  var tbl = btn.closest('.tbl-wrap').querySelector('table');
  var genes = Array.prototype.filter.call(tbl.tBodies[0].rows, function(tr) {
    return tr.style.display !== 'none';
  }).map(function(tr) { return tr.getAttribute('data-gene'); });
  Shiny.setInputValue('dl_flt', {genes: genes, t: Math.random()},
                      {priority: 'event'});
}
)---"

# Plot point hover -> highlighted row(s) in the table(s)
plot_js <- r"---(
function(el, x) {
  function clearRows() {
    document.querySelectorAll('tr.hl-row').forEach(function(r) {
      r.classList.remove('hl-row');
    });
  }
  if (el.removeAllListeners) {
    el.removeAllListeners('plotly_hover');
    el.removeAllListeners('plotly_unhover');
  }
  el.on('plotly_hover', function(d) {
    var t = d.points[0].text;
    if (!t) return;
    var gene = String(t).split('<br>')[0];
    clearRows();
    var rows = document.querySelectorAll('tr[data-gene]');
    for (var i = 0; i < rows.length; i++) {
      if (rows[i].getAttribute('data-gene') === gene) {
        rows[i].classList.add('hl-row');
        var box = rows[i].closest('.tbl-box');
        if (box) box.scrollTop = rows[i].offsetTop - box.clientHeight / 2;
      }
    }
  });
  el.on('plotly_unhover', clearRows);
  if (window.tblSync) window.tblSync();
}
)---"

# UI ___________________________________________________________________
ui <- fluidPage(
  tags$head(tags$style(HTML(app_css)), tags$script(HTML(app_js))),
  titlePanel("Volcano explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Data"),
      uiOutput("ui_base"),
      uiOutput("ui_over"),
      uiOutput("ui_restrict"),
      fluidRow(
        column(6, actionButton("clear_over", "Clear overlay",
                               width = "100%")),
        column(6, actionButton("clear_all", "Clear all", width = "100%"))
      ),
      uiOutput("exp_settings"),
      h4("2. Gene lists to highlight and link"),
      tags$div(
        id = "layers_box",
        layer_ui("L1", "List 1", layer_cols[1]),
        tags$div(id = "layers_end")
      ),
      actionButton("add_layer", "+ Add a gene list"),
      h4("3. Display"),
      textInput("leg_title", "Legend title", "Experiment"),
      fluidRow(
        column(6, numericInput("padj", "padj cutoff", 0.05, step = 0.01)),
        column(6, numericInput("lfc", "|log2FC| cutoff", 1, step = 0.5))
      ),
      fluidRow(
        column(6, numericInput("xmin", "X min (empty = auto)", NA)),
        column(6, numericInput("xmax", "X max (empty = auto)", NA))
      ),
      numericInput("ymax", "Y max (empty = auto)", NA),
      sliderInput("lab_size", "Gene name size", min = 2, max = 6,
                  value = 3, step = 0.5),
      checkboxInput("link_flt", "Plot only the genes passing the main table filters",
                    TRUE),
      fluidRow(
        column(6, sliderInput("ring", "List ring thickness", min = 0,
                              max = 2, value = 0.5, step = 0.1)),
        column(6, sliderInput("pt_size", "Listed point size", min = 1,
                              max = 5, value = 2, step = 0.5))
      ),
      h4("4. Export"),
      tags$button(
        type = "button", class = "btn btn-default",
        onclick = "volcanoPNG();", "Plot (.png)"
      ),
      actionButton("dl_tbl", "Highlighted genes (.csv)"),
      actionButton("dl_sel", "Selected on plot (.csv)"),
      actionButton("dl_all", "Full merged table (.csv)")
    ),
    mainPanel(
      plotlyOutput("volcano", height = 600),
      uiOutput("summary"),
      uiOutput("tbl_ui"),
      h4("Selected on the plot (lasso / box select)"),
      textOutput("sel_summary"),
      actionButton("sel_to_list", "Use selection as gene list"),
      tags$button(
        type = "button", class = "btn btn-default",
        onclick = paste0(
          "volcanoClearSel();",
          "Shiny.setInputValue('sel_clear', Math.random(),",
          " {priority: 'event'});"
        ),
        "Clear selection"
      ),
      uiOutput("sel_ui")
    )
  )
)

# SERVER _______________________________________________________________
server <- function(input, output, session) {

  # FILE INPUTS ________________________________________________________
  # A fileInput cannot be emptied from the server: the clear buttons bump
  # a counter that is part of the input ids, which rebuilds fresh inputs.
  # over_files accumulates the overlay files; layers holds the id and the
  # default color of every gene list block.
  rv <- reactiveValues(
    all = 0, over = 0, over_files = NULL,
    counter = 1, layers = c(L1 = layer_cols[1])
  )

  output$ui_base <- renderUI({
    fileInput(paste0("file_base_", rv$all), "Base experiment (DESeq2 table)")
  })
  output$ui_over <- renderUI({
    fileInput(paste0("files_overlay_", rv$over),
              "Add overlay experiment(s)", multiple = TRUE)
  })
  output$ui_restrict <- renderUI({
    fileInput(paste0("restrict_", rv$all),
              "Restrict to genes of a list (optional)")
  })

  f_base     <- reactive(input[[paste0("file_base_", rv$all)]])
  f_over     <- reactive(input[[paste0("files_overlay_", rv$over)]])
  f_restrict <- reactive(input[[paste0("restrict_", rv$all)]])

  # Each overlay upload is added to the previous ones (same name = replace)
  observeEvent(f_over(), {
    new <- f_over()[, c("name", "datapath")]
    old <- rv$over_files
    if (!is.null(old)) old <- old[!old$name %in% new$name, , drop = FALSE]
    rv$over_files <- rbind(old, new)
  })

  # Remove one overlay (index in the overlay store)
  observeEvent(input$rm_over, {
    i  <- as.integer(input$rm_over)
    ov <- rv$over_files
    if (!is.null(ov) && i >= 1 && i <= nrow(ov)) {
      rv$over_files <- ov[-i, , drop = FALSE]
    }
  })

  observeEvent(input$clear_over, {
    rv$over <- rv$over + 1
    rv$over_files <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$clear_all, {
    rv$all  <- rv$all + 1
    rv$over <- rv$over + 1
    rv$over_files <- NULL
    for (id in names(rv$layers)) {
      updateTextAreaInput(session, paste0("genes_", id), value = "")
    }
  }, ignoreInit = TRUE)

  # GENE LIST BLOCKS ___________________________________________________
  # An uploaded list file is appended to the text box of its block
  register_layer <- function(id) {
    observeEvent(input[[paste0("lfile_", id)]], {
      f   <- input[[paste0("lfile_", id)]]
      txt <- paste(readLines(f$datapath), collapse = "\n")
      cur <- input[[paste0("genes_", id)]] %||% ""
      all <- c(cur, txt)
      updateTextAreaInput(session, paste0("genes_", id),
                          value = paste(all[nzchar(all)], collapse = "\n"))
    })
  }

  add_layer <- function(name = NULL, text = "") {
    rv$counter <- rv$counter + 1
    id  <- paste0("L", rv$counter)
    col <- layer_cols[length(rv$layers) %% length(layer_cols) + 1]
    rv$layers <- c(rv$layers, setNames(col, id))
    insertUI(
      selector = "#layers_end", where = "beforeBegin",
      ui = layer_ui(id, name %||% paste("List", rv$counter), col, text)
    )
    register_layer(id)
  }

  register_layer("L1")
  observeEvent(input$add_layer, add_layer(), ignoreInit = TRUE)

  observeEvent(input$rm_layer, {
    id <- input$rm_layer
    removeUI(selector = paste0("#layer_", id), immediate = TRUE)
    rv$layers <- rv$layers[names(rv$layers) != id]
  })

  # DATA _______________________________________________________________
  # Base file first, then overlay files; names come from the file names
  raw_list <- reactive({
    req(f_base())
    f  <- f_base()[, c("name", "datapath")]
    ov <- rv$over_files
    if (!is.null(ov) && nrow(ov) > 0) f <- rbind(f, ov)
    nms <- make.unique(tools::file_path_sans_ext(f$name))
    setNames(lapply(f$datapath, read_table_auto), nms)
  })

  base_name <- reactive(names(raw_list())[1])
  sel_exps  <- reactive(names(raw_list()))

  # Legend names chosen in the UI (file name when left empty)
  exp_labels <- reactive({
    lab <- vapply(sel_exps(), function(nm) {
      v <- input[[paste0("lab_", exp_id(nm))]]
      if (is.null(v) || !nzchar(trimws(v))) nm else trimws(v)
    }, character(1))
    make.unique(unname(lab))
  })

  # Per-experiment settings: legend name, color picker, column mapping
  output$exp_settings <- renderUI({
    raw <- raw_list()
    # Keep the current value when the panel is redrawn
    keep <- function(id, default) {
      v <- isolate(input[[id]])
      if (is.null(v)) default else v
    }
    rows <- lapply(seq_along(raw), function(i) {
      nm   <- names(raw)[i]
      id   <- exp_id(nm)
      cols <- names(raw[[i]])
      g    <- guess_cols(cols)
      tags$div(
        style = "border-top: 1px solid #ddd; padding: 6px 0;",
        if (i > 1) tags$button(
          type = "button", class = "btn btn-default btn-xs",
          style = "float: right;",
          onclick = sprintf(
            "Shiny.setInputValue('rm_over', %d, {priority: 'event'});",
            i - 1),
          "Remove"
        ),
        tags$strong(paste0(if (i == 1) "Base: " else "Overlay: ", nm)),
        tags$input(
          type = "color", value = keep(paste0("col_", id), default_col(i)),
          style = "margin-left: 8px;",
          onchange = sprintf("Shiny.setInputValue('col_%s', this.value)", id)
        ),
        textInput(paste0("lab_", id), "Legend name",
                  value = keep(paste0("lab_", id), nm), width = "100%"),
        tags$details(
          tags$summary("Columns"),
          selectInput(paste0("id_", id), "Gene ID", cols,
                      selected = keep(paste0("id_", id), g$id),
                      width = "100%"),
          selectInput(paste0("lfc_", id), "log2FC", cols,
                      selected = keep(paste0("lfc_", id), g$lfc),
                      width = "100%"),
          selectInput(paste0("padj_", id), "padj", cols,
                      selected = keep(paste0("padj_", id), g$padj),
                      width = "100%")
        )
      )
    })
    tagList(h4("Colors, names and columns"), rows)
  })

  # Standardised tables; columns are guessed per file, editable in the UI
  tabs <- reactive({
    raw <- raw_list()
    lapply(setNames(names(raw), names(raw)), function(nm) {
      df   <- raw[[nm]]
      cols <- names(df)
      g    <- guess_cols(cols)
      id   <- exp_id(nm)
      pick <- function(key) {
        v <- input[[paste0(key, "_", id)]]
        if (!is.null(v) && v %in% cols) v else g[[key]]
      }
      out <- std_table(df, pick("id"), pick("lfc"), pick("padj"))
      validate(need(
        any(is.finite(out$lfc)),
        paste0("No numeric log2FC in '", nm, "': check its columns")
      ))
      out
    })
  })

  # Optional restriction to a list of genes (e.g. genes with a DAP peak)
  restricted <- reactive({
    tb <- tabs()
    f  <- f_restrict()
    rg <- if (is.null(f)) character(0) else
      read_gene_list(paste(readLines(f$datapath), collapse = "\n"))
    if (length(rg) > 0) {
      tb <- lapply(tb, function(d) filter(d, gene %in% rg))
    }
    tb
  })

  # GENE LISTS _________________________________________________________
  # Settings of every gene list block (name, color, source, typed genes)
  layer_defs <- reactive({
    ids <- names(rv$layers)
    lapply(seq_along(ids), function(i) {
      id <- ids[i]
      nm <- trimws(input[[paste0("name_", id)]] %||% "")
      list(
        id    = id,
        name  = if (nzchar(nm)) nm else paste("List", i),
        color = input[[paste0("lcol_", id)]] %||% rv$layers[[id]],
        src   = input[[paste0("src_", id)]] %||% "list",
        n     = max(1, num_or(input[[paste0("n_", id)]], 20)),
        genes = read_gene_list(input[[paste0("genes_", id)]] %||% "")
      )
    })
  })

  # Gene set of every list (needs the data for the DEG-based sources)
  layers_eval <- reactive({
    tb    <- restricted()
    exps  <- sel_exps()
    b     <- tb[[base_name()]]
    b_deg <- b[is_deg(b$padj, b$lfc, input$padj, input$lfc), ]
    common <- Reduce(intersect, lapply(exps, function(e) {
      d <- tb[[e]]
      d$gene[is_deg(d$padj, d$lfc, input$padj, input$lfc)]
    }))
    lapply(layer_defs(), function(L) {
      L$set <- switch(
        L$src,
        list   = L$genes,
        common = common,
        top_p  = b_deg %>% arrange(padj) %>%
          slice_head(n = L$n) %>% pull(gene),
        top_fc = b_deg %>% arrange(desc(abs(lfc))) %>%
          slice_head(n = L$n) %>% pull(gene)
      )
      L
    })
  })

  all_genes <- reactive({
    unique(as.character(unlist(lapply(layers_eval(), function(l) l$set))))
  })

  # For every listed gene, the names of the lists it belongs to
  membership <- reactive({
    g <- all_genes()
    L <- layers_eval()
    lists <- vapply(g, function(x) {
      hit <- vapply(L, function(l) x %in% l$set, logical(1))
      paste(vapply(L[hit], function(l) l$name, character(1)),
            collapse = ", ")
    }, character(1))
    tibble(gene = g, lists = unname(lists))
  })

  # One row per gene, one block of columns per experiment
  wide <- reactive({
    tb   <- restricted()
    exps <- sel_exps()
    leg  <- exp_labels()
    parts <- lapply(seq_along(exps), function(i) {
      d <- tb[[exps[i]]]
      out <- data.frame(
        gene = d$gene, lfc = d$lfc, padj = d$padj,
        deg = is_deg(d$padj, d$lfc, input$padj, input$lfc)
      )
      names(out)[-1] <- paste0(c("lfc", "padj", "deg"), "_", leg[i])
      out
    })
    Reduce(function(x, y) full_join(x, y, by = "gene"), parts)
  })

  # Crossing table of all listed genes
  report <- reactive({
    out <- membership() %>% left_join(wide(), by = "gene")
    lfc_cols <- grep("^lfc_", names(out), value = TRUE)
    out$n_found <- rowSums(!is.na(out[lfc_cols]))
    out$status <- case_when(
      out$n_found == 0 ~ "not found",
      out$n_found < length(lfc_cols) ~ "found in some experiments only",
      TRUE ~ "found in all"
    )
    out
  })

  # Plot coordinates (JSON) of some genes, used by the table rows to
  # draw a ring on the plot when hovered
  gene_pts <- function(g) {
    tb <- restricted()
    d  <- bind_rows(lapply(sel_exps(), function(e) {
      filter(tb[[e]], gene %in% g, is.finite(lfc), is.finite(nlp))
    }))
    lapply(split(d, d$gene), function(x) {
      jsonlite::toJSON(unname(Map(c, x$lfc, x$nlp)), digits = NA)
    })
  }

  # PLOT _______________________________________________________________
  build_plot <- function() {
    validate(need(f_base(), "Upload the base experiment first."))
    tb     <- restricted()
    exps   <- sel_exps()
    leg    <- exp_labels()
    bn     <- leg[1]
    layers <- layers_eval()

    # Active filters of the main table: keep only the genes passing them
    fv <- input$tbl_vis
    if (isTRUE(fv$on)) {
      keep_g <- as.character(unlist(fv$genes))
      layers <- lapply(layers, function(l) {
        l$set <- intersect(l$set, keep_g)
        l
      })
    }

    # Genes without finite values (padj NA or 0) cannot be drawn
    long <- bind_rows(lapply(seq_along(exps), function(i) {
      mutate(tb[[exps[i]]], exp = factor(leg[i], levels = leg))
    })) %>%
      left_join(membership(), by = "gene") %>%
      filter(is.finite(lfc), is.finite(nlp)) %>%
      mutate(
        deg = is_deg(padj, lfc, input$padj, input$lfc),
        tip = paste0(
          gene, "<br>", exp, "<br>log2FC: ", round(lfc, 2),
          "<br>padj: ", signif(padj, 3),
          ifelse(is.na(lists), "", paste0("<br>Lists: ", lists))
        )
      )

    # Colors of the experiments (default palette until one is picked)
    pal <- setNames(
      vapply(seq_along(exps), function(i) {
        v <- input[[paste0("col_", exp_id(exps[i]))]]
        if (is.null(v)) default_col(i) else v
      }, character(1)),
      leg
    )

    p <- ggplot(long, aes(lfc, nlp, color = exp, alpha = deg)) +
      geom_point(mapping = aes(text = tip, key = gene), size = 0.6) +
      scale_color_manual(values = pal, name = input$leg_title) +
      scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.8),
                         name = "DEG") +
      geom_hline(yintercept = -log10(input$padj), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1) * input$lfc, linetype = "dashed") +
      labs(x = "log2 fold change", y = "-log10(padj)") +
      theme_minimal()

    # Per-list display switches (inputs are NULL until the block is drawn)
    flag <- function(prefix, id, default) {
      v <- input[[paste0(prefix, "_", id)]]
      if (is.null(v)) default else isTRUE(v)
    }
    show_l  <- vapply(layers, function(l) flag("show",  l$id, TRUE),
                      logical(1))
    links_l <- vapply(layers, function(l) flag("links", l$id, TRUE),
                      logical(1))
    names_l <- vapply(layers, function(l) flag("names", l$id, FALSE),
                      logical(1))
    vis <- layers[show_l]   # lists whose highlight is visible

    # Halo of each visible gene list (concentric if a gene is in several)
    pt_size <- input$pt_size %||% 2
    ring    <- input$ring %||% 0.5
    halos   <- list()
    for (i in seq_along(vis)) {
      hl <- filter(long, gene %in% vis[[i]]$set)
      if (nrow(hl) == 0 || ring <= 0) next
      rk <- unname(vapply(hl$gene, function(x) {
        sum(vapply(vis[seq_len(i)], function(l) x %in% l$set, logical(1)))
      }, numeric(1)))
      for (r in unique(rk)) {
        halos[[length(halos) + 1]] <- list(
          data  = hl[rk == r, ],
          size  = pt_size + 2 * ring * r,
          color = vis[[i]]$color
        )
      }
    }
    # Outermost halo first, so that the inner ones stay visible
    sizes <- vapply(halos, function(h) h$size, numeric(1))
    for (h in halos[order(-sizes)]) {
      p <- p + geom_point(
        data = h$data,
        mapping = aes(x = lfc, y = nlp, text = tip, key = gene),
        inherit.aes = FALSE, size = h$size, color = h$color
      )
    }

    # Listed points of visible lists, in the color of their experiment
    vis_genes <- unique(as.character(unlist(lapply(vis, function(l) l$set))))
    hl_all <- filter(long, gene %in% vis_genes)
    if (nrow(hl_all) > 0) {
      p <- p + geom_point(
        data = hl_all,
        mapping = aes(x = lfc, y = nlp, color = exp, text = tip, key = gene),
        inherit.aes = FALSE, size = pt_size, alpha = 0.9
      )
    }

    # Links and names, each switchable per list
    labeled <- character(0)
    for (i in seq_along(layers)) {
      if (!links_l[i] && !names_l[i]) next
      L  <- layers[[i]]
      hl <- filter(long, gene %in% L$set)
      if (nrow(hl) == 0) next

      if (links_l[i]) {
        # Segments from the base experiment to each other experiment
        seg <- hl %>%
          filter(exp != bn) %>%
          select(gene, xend = lfc, yend = nlp) %>%
          inner_join(
            hl %>% filter(exp == bn) %>% select(gene, x = lfc, y = nlp),
            by = "gene"
          )
        if (nrow(seg) > 0) {
          p <- p + geom_segment(
            data = seg, aes(x = x, y = y, xend = xend, yend = yend),
            inherit.aes = FALSE, color = L$color,
            linewidth = 0.4, alpha = 0.7
          )
        }
      }

      if (names_l[i]) {
        # One label per gene, in the color of its first list
        lab <- hl %>%
          filter(!gene %in% labeled) %>%
          arrange(exp != bn) %>%
          distinct(gene, .keep_all = TRUE)
        labeled <- c(labeled, lab$gene)
        p <- p + geom_text(
          data = lab, mapping = aes(x = lfc, y = nlp, label = gene),
          inherit.aes = FALSE, color = L$color,
          size = input$lab_size %||% 3, vjust = -0.8
        )
      }
    }

    # Zoom without removing points; empty fields keep automatic limits
    custom <- c(input$xmin, input$xmax, input$ymax)
    if (nrow(long) > 0 && any(!is.na(custom))) {
      x_rng <- range(long$lfc)
      p <- p + coord_cartesian(
        xlim = c(num_or(input$xmin, x_rng[1]),
                 num_or(input$xmax, x_rng[2])),
        ylim = c(0, num_or(input$ymax, max(long$nlp)))
      )
    }
    p
  }

  # OUTPUTS ____________________________________________________________
  # Remove toWebGL() if labels render badly (slower with many genes)
  output$volcano <- renderPlotly({
    ggplotly(build_plot(), tooltip = "text", source = "volcano") %>%
      toWebGL() %>%
      event_register("plotly_selected") %>%
      htmlwidgets::onRender(plot_js)
  })

  output$tbl_ui <- renderUI({
    r <- report()
    req(nrow(r) > 0)
    gene_table(r, gene_pts(r$gene), main = TRUE)
  })

  # One summary line per gene list
  output$summary <- renderUI({
    w  <- wide()
    lc <- grep("^lfc_", names(w), value = TRUE)
    in_all <- w$gene[rowSums(!is.na(w[lc])) == length(lc)]
    tagList(lapply(layers_eval(), function(l) {
      n_all <- sum(l$set %in% in_all)
      n_any <- sum(l$set %in% w$gene)
      tags$div(
        tags$span(style = paste0("color: ", l$color, ";"), "\u25CF "),
        tags$strong(l$name), ": ", length(l$set), " genes (",
        n_all, " in all experiments, ", n_any - n_all,
        " in some only, ", length(l$set) - n_any, " not found)"
      )
    }))
  })

  # LASSO SELECTION ____________________________________________________
  # Selection snapshot hidden by the "Clear selection" button
  cleared <- reactiveVal(NULL)
  observeEvent(input$sel_clear, {
    cleared(event_data("plotly_selected", source = "volcano"))
  })

  # Genes under the lasso / box select (the plot sends them via "key")
  sel_genes <- reactive({
    ev <- event_data("plotly_selected", source = "volcano")
    if (is.null(ev) || !"key" %in% names(ev)) return(character(0))
    if (identical(ev, cleared())) return(character(0))
    unique(na.omit(as.character(ev$key)))
  })

  # Merged table restricted to the selected genes
  selected_df <- reactive({
    filter(wide(), gene %in% sel_genes())
  })

  output$sel_summary <- renderText({
    n <- length(sel_genes())
    if (n == 0) {
      "No selection yet: use the lasso or box select tool on the plot."
    } else if (n > 200) {
      paste0(n, " genes selected (first 200 shown, download for all)")
    } else {
      paste0(n, " genes selected")
    }
  })

  output$sel_ui <- renderUI({
    r <- selected_df()
    req(nrow(r) > 0)
    r <- head(r, 200)
    gene_table(r, gene_pts(r$gene))
  })

  # Turn the selection into a new gene list
  observeEvent(input$sel_to_list, {
    g <- sel_genes()
    req(length(g) > 0)
    add_layer("Selection", paste(g, collapse = "\n"))
  })

  # DOWNLOADS __________________________________________________________
  # Done in the browser (no downloadHandler): it needs a service worker,
  # which is not always available on static hosting.
  send_csv <- function(df, filename) {
    txt <- paste(capture.output(write.csv(df, row.names = FALSE)),
                 collapse = "\n")
    session$sendCustomMessage("save_text",
                              list(filename = filename, text = txt))
  }
  observeEvent(input$dl_tbl, send_csv(report(), "highlighted_genes.csv"))
  observeEvent(input$dl_sel, send_csv(selected_df(), "selected_on_plot.csv"))
  observeEvent(input$dl_all, send_csv(wide(), "merged_table.csv"))
  # Main table as displayed: current filters and sort order
  observeEvent(input$dl_flt, {
    g <- as.character(unlist(input$dl_flt$genes))
    r <- report()
    send_csv(r[match(g, r$gene), ], "filtered_table.csv")
  })
}

shinyApp(ui, server)