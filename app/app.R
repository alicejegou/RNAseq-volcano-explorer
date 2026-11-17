# Volcano explorer: overlay several DESeq2 results and link genes
# Upload a base DESeq2 table and overlay tables, then add genes.
# Runs in the browser with Shinylive: no server, data never leaves the
# user's computer.

library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)
library(plotly)

# HELPERS ______________________________________________________________

# Default colors: base first, then one per overlay (recycled)
default_cols <- c("#41ba43", "#c20845", "#1f77b4", "#ff7f0e",
                  "#9467bd", "#8c564b", "#e377c2", "#17becf")
default_col <- function(i) {
  default_cols[(i - 1) %% length(default_cols) + 1]
}

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

# Gene list from a text box and/or a file (space, comma, ; or newline)
read_gene_list <- function(txt = "", file = NULL) {
  if (!is.null(file)) {
    txt <- paste(txt, paste(readLines(file$datapath), collapse = "\n"))
  }
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

# FRONT-END CODE _______________________________________________________

app_css <- "
.gene-table th { position: sticky; top: 0; background: #fff; }
.gene-table tr[data-gene]:hover td,
.gene-table tr.hl-row td { background-color: #fff3b0 !important; }
"

# Table row hover -> ring around the gene's points on the plot
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
  Plotly.addTraces(gd, {
    x: pts.map(function(p) { return p[0]; }),
    y: pts.map(function(p) { return p[1]; }),
    type: 'scatter', mode: 'markers', name: '__hl',
    hoverinfo: 'skip', showlegend: false,
    marker: {
      size: 16, color: 'rgba(0,0,0,0)',
      line: {color: 'black', width: 3}
    }
  });
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
)---"

# Plot point hover -> highlighted row in the table (runs on the widget)
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
        var box = document.getElementById('tbl_box');
        if (box) box.scrollTop = rows[i].offsetTop - box.clientHeight / 2;
        break;
      }
    }
  });
  el.on('plotly_unhover', clearRows);
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
      h4("2. Genes to highlight and link"),
      radioButtons(
        "mode", NULL,
        choices = c(
          "My gene list" = "list",
          "DEG common to all experiments" = "common",
          "Top N DEG of base, by padj" = "top_p",
          "Top N DEG of base, by |log2FC|" = "top_fc"
        )
      ),
      conditionalPanel(
        "input.mode == 'list'",
        textAreaInput("genes", NULL, rows = 6,
                      placeholder = "Mp1g02140\nMp4g13300"),
        uiOutput("ui_gfile")
      ),
      conditionalPanel(
        "input.mode == 'top_p' || input.mode == 'top_fc'",
        numericInput("ntop", "N", 20, min = 1)
      ),
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
      h4("4. Export"),
      downloadButton("dl_png", "Plot (.png)"),
      downloadButton("dl_tbl", "Highlighted genes (.csv)"),
      downloadButton("dl_all", "Full merged table (.csv)")
    ),
    mainPanel(
      plotlyOutput("volcano", height = 600),
      textOutput("summary"),
      uiOutput("tbl_ui")
    )
  )
)

# SERVER _______________________________________________________________
server <- function(input, output, session) {
  
  # FILE INPUTS ________________________________________________________
  # A fileInput cannot be emptied from the server: the clear buttons bump
  # a counter that is part of the input ids, which rebuilds fresh inputs.
  rv <- reactiveValues(all = 0, over = 0)
  
  observeEvent(input$clear_over, {
    rv$over <- rv$over + 1
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_all, {
    rv$all  <- rv$all + 1
    rv$over <- rv$over + 1
    updateTextAreaInput(session, "genes", value = "")
  }, ignoreInit = TRUE)
  
  output$ui_base <- renderUI({
    fileInput(paste0("file_base_", rv$all), "Base experiment (DESeq2 table)")
  })
  output$ui_over <- renderUI({
    fileInput(paste0("files_overlay_", rv$over), "Overlay experiment(s)",
              multiple = TRUE)
  })
  output$ui_restrict <- renderUI({
    fileInput(paste0("restrict_", rv$all),
              "Restrict to genes of a list (optional)")
  })
  output$ui_gfile <- renderUI({
    fileInput(paste0("gfile_", rv$all), "...or upload a list")
  })
  
  f_base     <- reactive(input[[paste0("file_base_", rv$all)]])
  f_over     <- reactive(input[[paste0("files_overlay_", rv$over)]])
  f_restrict <- reactive(input[[paste0("restrict_", rv$all)]])
  f_gene     <- reactive(input[[paste0("gfile_", rv$all)]])
  
  # DATA _______________________________________________________________
  # Base file first, then overlay files; names come from the file names
  raw_list <- reactive({
    req(f_base())
    f <- f_base()[, c("name", "datapath")]
    if (!is.null(f_over())) {
      f <- rbind(f, f_over()[, c("name", "datapath")])
    }
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
    rg <- read_gene_list("", f_restrict())
    if (length(rg) > 0) {
      tb <- lapply(tb, function(d) filter(d, gene %in% rg))
    }
    tb
  })
  
  user_genes <- reactive(read_gene_list(input$genes, f_gene()))
  
  # Genes to highlight, depending on the chosen mode
  target_genes <- reactive({
    tb    <- restricted()
    exps  <- sel_exps()
    b     <- tb[[base_name()]]
    b_deg <- b[is_deg(b$padj, b$lfc, input$padj, input$lfc), ]
    switch(
      input$mode,
      list = user_genes(),
      common = Reduce(intersect, lapply(exps, function(e) {
        d <- tb[[e]]
        d$gene[is_deg(d$padj, d$lfc, input$padj, input$lfc)]
      })),
      top_p = b_deg %>% arrange(padj) %>%
        slice_head(n = input$ntop) %>% pull(gene),
      top_fc = b_deg %>% arrange(desc(abs(lfc))) %>%
        slice_head(n = input$ntop) %>% pull(gene)
    )
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
  
  # Crossing table of the highlighted genes
  report <- reactive({
    out <- tibble(gene = target_genes()) %>%
      left_join(wide(), by = "gene")
    lfc_cols <- grep("^lfc_", names(out), value = TRUE)
    out$n_found <- rowSums(!is.na(out[lfc_cols]))
    out$status <- case_when(
      out$n_found == 0 ~ "not found",
      out$n_found < length(lfc_cols) ~ "found in some experiments only",
      TRUE ~ "found in all"
    )
    out
  })
  
  # Plot coordinates of each highlighted gene (JSON), used by the table
  # rows to draw a ring on the plot when hovered
  tbl_pts <- reactive({
    tb <- restricted()
    g  <- target_genes()
    d  <- bind_rows(lapply(sel_exps(), function(e) {
      filter(tb[[e]], gene %in% g, is.finite(lfc), is.finite(nlp))
    }))
    lapply(split(d, d$gene), function(x) {
      jsonlite::toJSON(unname(Map(c, x$lfc, x$nlp)), digits = NA)
    })
  })
  
  # PLOT _______________________________________________________________
  # Build the volcano plot
  #   interactive = TRUE  : plotly version (hover text, plain labels)
  #   interactive = FALSE : static version (ggrepel labels), used for PNG
  build_plot <- function(interactive = FALSE) {
    validate(need(f_base(), "Upload the base experiment first."))
    tb   <- restricted()
    exps <- sel_exps()
    leg  <- exp_labels()
    bn   <- leg[1]
    
    # Genes without finite values (padj NA or 0) cannot be drawn
    long <- bind_rows(lapply(seq_along(exps), function(i) {
      mutate(tb[[exps[i]]], exp = factor(leg[i], levels = leg))
    })) %>%
      filter(is.finite(lfc), is.finite(nlp)) %>%
      mutate(
        deg = is_deg(padj, lfc, input$padj, input$lfc),
        tip = paste0(gene, "<br>", exp, "<br>log2FC: ", round(lfc, 2),
                     "<br>padj: ", signif(padj, 3))
      )
    
    # Colors chosen in the UI (default palette until one is picked)
    pal <- setNames(
      vapply(seq_along(exps), function(i) {
        v <- input[[paste0("col_", exp_id(exps[i]))]]
        if (is.null(v)) default_col(i) else v
      }, character(1)),
      leg
    )
    
    hl <- filter(long, gene %in% target_genes())
    
    # Segments from the base experiment to each other experiment
    seg <- hl %>%
      filter(exp != bn) %>%
      select(gene, xend = lfc, yend = nlp) %>%
      inner_join(
        hl %>% filter(exp == bn) %>% select(gene, x = lfc, y = nlp),
        by = "gene"
      )
    
    # One label per gene, taken from the base experiment when possible
    lab <- hl %>%
      arrange(exp != bn) %>%
      distinct(gene, .keep_all = TRUE)
    
    # The "text" aesthetic (hover) is only valid for plotly
    map_all <- if (interactive) aes(text = tip) else aes()
    map_hl <- if (interactive) {
      aes(x = lfc, y = nlp, color = exp, text = tip)
    } else {
      aes(x = lfc, y = nlp, color = exp)
    }
    
    # ggrepel is not supported by plotly: plain labels when interactive
    lab_map <- aes(x = lfc, y = nlp, label = gene)
    lab_layer <- if (interactive) {
      geom_text(data = lab, mapping = lab_map, inherit.aes = FALSE,
                size = 3, vjust = -0.8)
    } else {
      geom_text_repel(data = lab, mapping = lab_map, inherit.aes = FALSE,
                      size = 2.5, max.overlaps = 30)
    }
    
    p <- ggplot(long, aes(lfc, nlp, color = exp, alpha = deg)) +
      geom_point(mapping = map_all, size = 0.6) +
      scale_color_manual(values = pal, name = input$leg_title) +
      scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.8),
                         name = "DEG") +
      geom_hline(yintercept = -log10(input$padj), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1) * input$lfc, linetype = "dashed") +
      labs(x = "log2 fold change", y = "-log10(padj)") +
      theme_minimal()
    
    if (nrow(hl) > 0) {
      p <- p +
        geom_segment(
          data = seg, aes(x = x, y = y, xend = xend, yend = yend),
          inherit.aes = FALSE, color = "#474747",
          linewidth = 0.3, alpha = 0.6
        ) +
        geom_point(
          data = hl, mapping = map_hl,
          inherit.aes = FALSE, size = 2, alpha = 0.9
        ) +
        lab_layer
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
    ggplotly(build_plot(interactive = TRUE), tooltip = "text") %>%
      toWebGL() %>%
      htmlwidgets::onRender(plot_js)
  })
  
  # Hand-made table: each row carries its gene and plot coordinates
  output$tbl_ui <- renderUI({
    r <- report()
    req(nrow(r) > 0)
    pts   <- tbl_pts()
    cells <- lapply(r, fmt_cell)
    body <- lapply(seq_len(nrow(r)), function(i) {
      p <- pts[[r$gene[i]]]
      tags$tr(
        `data-gene` = r$gene[i],
        `data-pts`  = if (is.null(p)) "[]" else as.character(p),
        unname(lapply(cells, function(col) tags$td(col[i])))
      )
    })
    tags$div(
      id = "tbl_box",
      style = "max-height: 320px; overflow: auto; position: relative;",
      tags$table(
        class = "table table-condensed table-striped gene-table",
        tags$thead(tags$tr(lapply(names(r), tags$th))),
        tags$tbody(body)
      )
    )
  })
  
  output$summary <- renderText({
    r <- report()
    paste0(
      nrow(r), " genes: ", sum(r$status == "found in all"),
      " found in all experiments, ",
      sum(r$status == "found in some experiments only"),
      " in some only, ", sum(r$status == "not found"), " not found"
    )
  })
  
  output$dl_png <- downloadHandler(
    filename = function() "volcano.png",
    content = function(file) {
      ggsave(file, build_plot(), width = 7, height = 5, dpi = 300)
    }
  )
  output$dl_tbl <- downloadHandler(
    filename = function() "highlighted_genes.csv",
    content = function(file) write.csv(report(), file, row.names = FALSE)
  )
  output$dl_all <- downloadHandler(
    filename = function() "merged_table.csv",
    content = function(file) write.csv(wide(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)