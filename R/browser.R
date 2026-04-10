
#' Plot coverages along genome
#'
#' Visualizes genomic coverage profiles with gene annotations. Supports both
#' the legacy interface (fstart/fend/fchr) and the preferred GRanges interface.
#'
#' @param ranges A GRanges object specifying one or more genomic regions to plot.
#'   If provided, this is the preferred way to specify regions. If multiple ranges
#'   are provided, one plot is created per range. If \code{ranges} has names or
#'   a "name" metadata column, these are used as plot titles.
#' @param fstart Deprecated. Start of the genomic window (bp). Use \code{ranges} instead.
#' @param fend Deprecated. End of the genomic window (bp). Use \code{ranges} instead.
#' @param fchr Deprecated. Chromosome of the genomic window. Use \code{ranges} instead.
#' @param profs Named list of coverages to be plotted (RleList objects, including SimpleRleList and CompressedRleList).
#' @param cols Colors of the profiles.
#' @param ann A named list of annotation dataframes. Each dataframe should have
#'   columns: chr, start, end, col (optional), label (optional).
#' @param ylabel The label of the y axis.
#' @param ylims List of ylimits for plotting the coverages: list(c(min, max), ...).
#' @param txdb A TxDb transcription database used for plotting gene models.
#' @param ftitle Title of the plot. If NA and ranges has names, uses range name.
#' @param collapse Collapse gene models to longest transcript (default = TRUE).
#' @param with.genes.highlited Vector of gene IDs that should be highlighted.
#' @param plot.labels Plot gene labels (default = TRUE).
#' @param grid Plot background grid (default = FALSE).
#' @param with.average Plot average (default = FALSE).
#'
#' @return For a single range: invisibly returns a gTree grob of the plot.
#'   For multiple ranges: returns a named list of gTree grobs (one per range).
#'   In both cases, plots are also drawn to the active graphics device.
#'
#' @details
#' The function draws coverage profiles as filled polygons overlaid with gene
#' model annotations from the TxDb database. When multiple ranges are provided,
#' each range is plotted on a new page.
#'
#' The \code{ranges} argument is the preferred way to specify genomic regions.
#' The legacy arguments (fstart, fend, fchr) are deprecated but still supported
#' for backwards compatibility. If both \code{ranges} and legacy arguments are
#' provided, a warning is issued and \code{ranges} takes precedence.
#'
#' @examples
#' \dontrun{
#' # Using GRanges (preferred)
#' library(GenomicRanges)
#' regions <- GRanges("chrX", IRanges(c(1660000, 1800000), c(1720000, 1850000)))
#' names(regions) <- c("Region A", "Region B")
#' plotProfiles(ranges = regions, profs = list(MSL2 = cov), txdb = txdb)
#'
#' # Legacy interface (deprecated)
#' plotProfiles(fstart = 1660000, fend = 1720000, fchr = "chrX",
#'              profs = list(MSL2 = cov), txdb = txdb)
#' }
#'
#' @export

plotProfiles <- function(ranges = NULL, fstart = NULL, fend = NULL, fchr = NULL,
                         profs, cols = c(), ann = NULL, ylabel = "coverage",
                         ylims = list(), txdb, ftitle = NA, collapse = TRUE,
                         with.genes.highlited = c(), plot.labels = TRUE,
                         grid = FALSE, with.average = FALSE) {

  require(grid)
  require(IRanges)
  require(GenomicRanges)
  require(HilbertVis)
  require(RColorBrewer)
  require(AnnotationDbi)

  # Handle input: GRanges vs legacy arguments
  if (!is.null(ranges)) {
    # Validate ranges
    if (!inherits(ranges, "GRanges")) {
      # Try to coerce
      tryCatch({
        ranges <- as(ranges, "GRanges")
      }, error = function(e) {
        stop("'ranges' must be a GRanges object or coerceable to GRanges: ", e$message)
      })
    }

    # Warn if legacy args also provided
    if (!is.null(fstart) || !is.null(fend) || !is.null(fchr)) {
      warning("'ranges' provided; ignoring deprecated arguments 'fstart', 'fend', 'fchr'")
    }

    n_ranges <- length(ranges)
  } else {
    # Legacy interface
    if (is.null(fstart) || is.null(fend) || is.null(fchr)) {
      stop("Either 'ranges' (preferred) or all of 'fstart', 'fend', 'fchr' must be provided")
    }
    ranges <- GRanges(fchr, IRanges(fstart, fend))
    n_ranges <- 1
  }

  # Get range names for titles
  range_names <- names(ranges)
  if (is.null(range_names) && "name" %in% colnames(mcols(ranges))) {
    range_names <- mcols(ranges)$name
  }

  # Single range: use ftitle if provided, else range name
  # Multiple ranges: use range names if available
  titles <- if (n_ranges == 1 && !is.na(ftitle)) {
    ftitle
  } else if (!is.null(range_names)) {
    range_names
  } else {
    rep(NA, n_ranges)
  }

  # Plot each range and collect grobs
  grobs <- vector("list", n_ranges)
  grob_names <- if (!is.null(range_names)) {
    range_names
  } else {
    paste0("range_", seq_len(n_ranges))
  }
  names(grobs) <- grob_names

  for (i in seq_len(n_ranges)) {
    r <- ranges[i]
    grobs[[i]] <- .plotProfiles_oneRange(
      fchr = as.character(seqnames(r)),
      fstart = start(r),
      fend = end(r),
      profs = profs,
      cols = cols,
      ann = ann,
      ylabel = ylabel,
      ylims = ylims,
      txdb = txdb,
      ftitle = titles[i],
      collapse = collapse,
      with.genes.highlited = with.genes.highlited,
      plot.labels = plot.labels,
      grid = grid,
      with.average = with.average
    )
  }

  # Return

  if (n_ranges == 1) {
    invisible(grobs[[1]])
  } else {
    invisible(grobs)
  }
}


#' Internal: plot a single genomic range
#' @keywords internal
.plotProfiles_oneRange <- function(fchr, fstart, fend, profs, cols, ann, ylabel,
                                    ylims, txdb, ftitle, collapse,
                                    with.genes.highlited, plot.labels, grid,
                                    with.average) {

  ############################################################################
  # Validate args
  ############################################################################
  if (!(fchr %in% unlist(lapply(profs, names)))) {
    stop("chromosome '", fchr, "' not found in coverage profiles")
  }

  txdb_seqnames <- tryCatch(
    levels(seqnames(genes(txdb))),
    error = function(e) character(0)
  )
  if (length(txdb_seqnames) > 0 && !(fchr %in% txdb_seqnames)) {
    warning("chromosome '", fchr, "' not found in txdb annotation")
  }

  options(scipen = 100)

  ############################################################################
  # HELPER FUNCTION DEFINITIONS
  ############################################################################
  profile.xscale <- function(xscale) {
    main.ticks <- grid.pretty(xscale)
    dist.main <- main.ticks[2] - main.ticks[1]
    minor.ticks <- seq(
      from = main.ticks[1] - dist.main,
      to = main.ticks[length(main.ticks)] + dist.main,
      by = dist.main / 5
    )
    minor.ticks <- minor.ticks[!(minor.ticks %in% main.ticks)]
    minor.ticks <- minor.ticks[minor.ticks > xscale[1] & minor.ticks < xscale[2]]
    for (i in minor.ticks) {
      grid.lines(c(i, i), c(1, 0.9), default.units = "native", gp = gpar(col = 1))
    }
    for (i in main.ticks) {
      grid.lines(c(i, i), c(1, 0.8), default.units = "native", gp = gpar(col = 1))
      grid.text(
        paste(formatC(i / 1000, format = "f", digits = 0,
                      big.mark = ".", decimal.mark = ","), "K"),
        i, 0.7, default.units = "native",
        just = c("center", "top"), gp = gpar(cex = 0.4)
      )
    }
    grid.lines(xscale, 1, default.units = "native", gp = gpar(col = 1))
  }

  bg.grid <- function(xscale) {
    main.ticks <- grid.pretty(xscale)
    dist.main <- main.ticks[2] - main.ticks[1]
    minor.ticks <- seq(
      from = main.ticks[1] - dist.main,
      to = main.ticks[length(main.ticks)] + dist.main,
      by = dist.main / 5
    )
    minor.ticks <- minor.ticks[!(minor.ticks %in% main.ticks)]
    minor.ticks <- minor.ticks[minor.ticks > xscale[1] & minor.ticks < xscale[2]]
    for (i in minor.ticks) {
      grid.lines(c(i, i), c(0, 1), default.units = "native",
                 gp = gpar(col = "lightgrey", lwd = 0.5))
    }
    for (i in main.ticks) {
      grid.lines(c(i, i), c(0, 1), default.units = "native",
                 gp = gpar(col = "lightgrey", lwd = 0.5))
    }
  }

  ts.grid.rect <- function(start, end, bottom, top, col = 1, lcol = 1, lwd = 1) {
    grid.rect(
      x = start, y = bottom,
      width = end - start, height = top - bottom,
      just = c("left", "bottom"), default.units = "native",
      gp = gpar(fill = col, col = lcol, lwd = lwd)
    )
  }

  plot.genes <- function(vp, strand, frame.genes, bumps, gene.height, label.shift,
                         collapse, with.genes.highlited, plot.labels,
                         gene_tx_map, frame_transcripts, exon_tx_map, tx_to_gene) {
    exon <- 2.0
    gene.frame <- frame.genes[strand(frame.genes) == strand]
    pushViewport(vp)
    if (length(gene.frame) > 0) {
      for (i in 1:length(gene.frame)) {
        current_gene_id <- names(gene.frame)[i]
        txts.names <- gene_tx_map$TXNAME[gene_tx_map$GENEID == current_gene_id]
        if (length(txts.names) == 0) next
        txts <- frame_transcripts[mcols(frame_transcripts)$tx_name %in% txts.names]
        if (length(txts) > 0) {
          if (collapse) {
            txts <- txts[which.max(width(txts))]
          }
          for (k in 1:length(txts)) {
            tx_name_k <- mcols(txts)$tx_name[k]
            gene.id <- tx_to_gene[tx_name_k]
            col <- ifelse(gene.id %in% with.genes.highlited, "#FFAAAA", "#F2F2F2")
            lcol <- ifelse(gene.id %in% with.genes.highlited, "#FF0000", "#000000")

            if (strand == "-") {
              t1 <- (bumps * 10) - ((k - 1) * gene.height + gene.height / 2 - label.shift / 2)
            } else {
              t1 <- ((k - 1) * gene.height + gene.height / 2) - (label.shift / 2)
            }
            ts.grid.rect(start(txts[k]), end(txts[k]), t1, t1, lcol = lcol, lwd = 2)

            exon_data <- exon_tx_map[exon_tx_map$TXNAME == tx_name_k, ]
            if (nrow(exon_data) > 0) {
              for (l in 1:nrow(exon_data)) {
                ts.grid.rect(exon_data$EXONSTART[l], exon_data$EXONEND[l],
                             t1 - exon, t1 + exon, col = col, lcol = lcol, lwd = 1)
              }
            }
            if (plot.labels) {
              label <- gene.id
              if (strand == "-") {
                grid.text(label, end(txts)[k], t1 - label.shift - 0.5,
                          default.units = "native", just = c("right", "bottom"),
                          gp = gpar(cex = 0.5))
              } else {
                grid.text(label, end(txts)[k], t1 + label.shift + 0.4,
                          default.units = "native", just = c("right", "top"),
                          gp = gpar(cex = 0.5))
              }
            }
          }
        }
      }
    }
    popViewport()
  }

  plot.profiles <- function(profs, fchr, fstart, fend, ylims, cols, xscale, with.average) {
    vsize <- as.integer(dev.size()[1] * 150)
    margins <- unit(0.19, "lines")
    font.size.label <- 0.7
    panel.background <- "#cccccc25"

    if (length(cols) == 0) {
      cols <- RColorBrewer::brewer.pal(ifelse(length(profs) > 2, length(profs), 3), "Dark2")
    }

    for (i in 1:length(profs)) {
      if (is(profs[[i]], "RleList")) {
        xl <- c(fstart, seq(fstart, fend, length.out = vsize), fend)
        yl <- c(0, HilbertVis::shrinkVector(
          as.vector(profs[[i]][[fchr]])[fstart:fend], newLength = vsize), 0)
        yl[is.na(yl)] <- 0

        if (length(ylims) == 0) {
          vmax <- max(pretty(0:ceiling(max(yl, na.rm = TRUE))))
        } else {
          vmax <- ylims[[i]][2]
        }

        # Clipped viewport
        prof1 <- viewport(
          x = 0, y = unit((length(profs) - i) / length(profs), "npc") + 2 * margins,
          w = 1, h = unit(1 / length(profs), "npc") - 2 * margins,
          just = c("left", "bottom"), xscale = xscale, yscale = c(0, vmax), clip = "on"
        )
        pushViewport(prof1)
        grid.rect(gp = gpar(fill = panel.background, lty = 0))
        grid.lines(xscale, 0, default.units = "native", gp = gpar(col = "grey66", lwd = 0.5))
        grid.polygon(x = xl, y = yl, default.units = "native",
                     gp = gpar(col = "#55555570", lwd = 0.5, fill = NA))
        grid.polygon(x = xl, y = yl, default.units = "native",
                     gp = gpar(col = NA, fill = cols[[i]]))
        grid.text(names(profs)[[i]], unit(0.4, "lines"), unit(1, "npc") - unit(0.3, "lines"),
                  just = c("left", "top"), gp = gpar(cex = font.size.label, font = 1))
        popViewport()

        # Unclipped viewport for border
        prof1 <- viewport(
          x = 0, y = unit((length(profs) - i) / length(profs), "npc") + 2 * margins,
          w = 1, h = unit(1 / length(profs), "npc") - 2 * margins,
          just = c("left", "bottom"), xscale = xscale, yscale = c(0, vmax), clip = "off"
        )
        pushViewport(prof1)
        grid.rect(gp = gpar(lwd = 1, fill = NA))
        popViewport()

        # Y-axis
        yax <- viewport(
          x = unit(-0.2, "lines"),
          y = unit((length(profs) - i) / length(profs), "npc") + 2 * margins,
          w = unit(0.2, "lines"), h = unit(1 / length(profs), "npc") - 2 * margins,
          just = c("left", "bottom"), xscale = xscale, yscale = c(0, vmax), clip = "off"
        )
        pushViewport(yax)
        grid.yaxis(gp = gpar(cex = 0.45))
        popViewport()
      }
    }
  }

  plot.annotation <- function(ann, fchr, fstart, fend, xscale) {
    margins <- unit(0.05, "lines")
    for (i in 1:length(ann)) {
      ann1 <- viewport(
        x = 0, y = unit((length(ann) - i) / length(ann), "npc") + 2 * margins,
        w = 1, h = unit(1 / length(ann), "npc") - 2 * margins,
        just = c("left", "bottom"), xscale = xscale, clip = "on", yscale = c(0, 1)
      )
      pushViewport(ann1)
      ca <- ann[[i]]
      ca <- ca[ca$chr == fchr & ca$end > fstart & ca$start < fend, ]
      if (nrow(ca) > 0) {
        for (k in 1:nrow(ca)) {
          col_val <- if (is.null(ca$col[k])) 1 else as.character(ca$col[k])
          ts.grid.rect(ca$start[k], ca$end[k], 0, 1, col = col_val, lcol = col_val, lwd = 1)
          if ("label" %in% names(ca)) {
            grid.text(as.character(ca$label[k]),
                      unit(ca$start[k], "native") - unit(0.003, "npc"),
                      unit(0.5, "npc"), just = c("right", "center"),
                      gp = gpar(cex = 0.5, font = 1))
          }
        }
      }
      popViewport()

      ann1 <- viewport(
        x = 0, y = unit((length(ann) - i) / length(ann), "npc") + 2 * margins,
        w = 1, h = unit(1 / length(ann), "npc") - 2 * margins,
        just = c("left", "bottom"), xscale = xscale, clip = "off", yscale = c(0, 1)
      )
      pushViewport(ann1)
      grid.text(names(ann)[i], unit(-0.01, "npc"), unit(0.5, "npc"),
                just = c("right", "center"), gp = gpar(cex = 0.5, font = 1))
      popViewport()
    }
  }

  ############################################################################
  # VAL DEFINITIONS
  ############################################################################
  gene.height <- 10
  gene.lines <- 1.3
  label.shift <- 6
  scale.y.offset <- 1.2
  xscale <- c(fstart, fend)

  ############################################################################
  # START PLOT
  ############################################################################
  grid.newpage()
  main <- viewport(
    x = unit(3, "lines"), y = unit(1, "lines"),
    width = unit(1, "npc") - unit(4, "lines"),
    height = unit(1, "npc") - unit(3, "lines"),
    just = c("left", "bottom"), xscale = xscale
  )
  pushViewport(main)

  # Background grid
  if (grid) {
    bg.grid(xscale)
  }

  # Get genes in frame
  frameRange <- GRanges(fchr, IRanges(fstart, fend))
  frame.genes <- subsetByOverlaps(genes(txdb), frameRange)
  frame.genes.rev <- frame.genes[strand(frame.genes) == "-"]
  frame.genes.fwd <- frame.genes[strand(frame.genes) == "+"]

  # Pre-fetch all gene-to-transcript mappings in one batch query
  gene_ids_in_frame <- names(frame.genes)
  if (length(gene_ids_in_frame) > 0) {
    gene_tx_map <- suppressMessages(
      AnnotationDbi::select(txdb, keys = gene_ids_in_frame,
                            columns = "TXNAME", keytype = "GENEID")
    )
  } else {
    gene_tx_map <- data.frame(GENEID = character(0), TXNAME = character(0))
  }

  # Pre-fetch all transcripts for the frame
  all_tx_names <- unique(gene_tx_map$TXNAME)
  if (length(all_tx_names) > 0) {
    frame_transcripts <- transcripts(txdb, filter = list(tx_name = all_tx_names))
  } else {
    frame_transcripts <- GRanges()
  }

  # Pre-fetch all exons for transcripts in frame
  if (length(all_tx_names) > 0) {
    exon_tx_map <- suppressMessages(
      AnnotationDbi::select(txdb, keys = all_tx_names,
                            columns = c("EXONID", "EXONSTART", "EXONEND"),
                            keytype = "TXNAME")
    )
  } else {
    exon_tx_map <- data.frame(TXNAME = character(0), EXONID = integer(0),
                               EXONSTART = integer(0), EXONEND = integer(0))
  }

  # Pre-compute tx_name to gene_id reverse mapping
  tx_to_gene <- setNames(gene_tx_map$GENEID, gene_tx_map$TXNAME)

  # Max number of transcript isoforms (using pre-fetched data)
  get.no.tracks <- function(genes, gene_tx_map) {
    if (length(genes) == 0) return(1)
    gene_ids <- names(genes)
    counts <- table(gene_tx_map$GENEID[gene_tx_map$GENEID %in% gene_ids])
    if (length(counts) == 0) return(1)
    max(counts)
  }

  rbumps <- ifelse(collapse, 1, get.no.tracks(frame.genes.rev, gene_tx_map))
  fbumps <- ifelse(collapse, 1, get.no.tracks(frame.genes.fwd, gene_tx_map))

  # X-scale
  scale <- viewport(
    x = 0, y = unit(gene.lines * rbumps + 0.5, "lines") - unit(scale.y.offset, "lines"),
    w = 1, h = unit(1.5, "lines"), clip = "off",
    just = c("left", "bottom"), name = "scale", xscale = xscale
  )
  pushViewport(scale)
  profile.xscale(xscale)
  popViewport()

  # Reverse strand genes
  rev <- viewport(
    x = 0, y = 0, w = 1, h = unit(gene.lines * rbumps, "lines"),
    clip = "on", just = c("left", "bottom"), name = "rev",
    xscale = xscale, yscale = c(-1, gene.height * rbumps + 1)
  )
  plot.genes(rev, "-", frame.genes, rbumps, gene.height, label.shift, collapse,
             with.genes.highlited, plot.labels, gene_tx_map, frame_transcripts,
             exon_tx_map, tx_to_gene)

  # Forward strand genes
  sval <- if (is.null(ann)) 0.5 else 0.7
  fwd <- viewport(
    x = 0, y = unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines"),
    w = 1, h = unit(gene.lines * fbumps, "lines"),
    clip = "on", just = c("left", "bottom"), name = "fwd",
    xscale = xscale, yscale = c(-1, gene.height * fbumps + 1)
  )
  plot.genes(fwd, "+", frame.genes, fbumps, gene.height, label.shift, collapse,
             with.genes.highlited, plot.labels, gene_tx_map, frame_transcripts,
             exon_tx_map, tx_to_gene)

  # Annotations
  if (!is.null(ann)) {
    annot <- viewport(
      x = 0, y = unit(gene.lines * rbumps + 0.5, "lines") + unit(0.5, "lines") +
        unit(gene.lines * fbumps, "lines"),
      w = 1, h = unit(0.5, "lines"),
      just = c("left", "bottom"), name = "annot", xscale = xscale
    )
    pushViewport(annot)
    plot.annotation(ann, fchr, fstart, fend, xscale)
    popViewport()
  }

  # Profiles
  sval <- if (is.null(ann)) 0.2 else 0.7
  prof <- viewport(
    x = 0,
    y = unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines") +
      unit(gene.lines * fbumps, "lines"),
    w = 1,
    h = unit(1, "npc") - (unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines") +
                            unit(gene.lines * fbumps, "lines")),
    just = c("left", "bottom"), name = "prof", xscale = xscale
  )
  pushViewport(prof)
  plot.profiles(profs, fchr, fstart, fend, ylims, cols, xscale, with.average)
  popViewport()

  # Y-axis label
  prof.ylab <- viewport(
    x = unit(-1.8, "lines"),
    y = unit(gene.lines * rbumps + 0.5, "lines") + unit(0.7, "lines") +
      unit(gene.lines * fbumps, "lines") + unit(0.2, "lines"),
    w = unit(0.4, "lines"),
    h = unit(1, "npc") - (unit(gene.lines * rbumps + 0.5, "lines") + unit(0.7, "lines") +
                            unit(gene.lines * fbumps, "lines") + unit(0.2, "lines")),
    just = c("left", "bottom"), clip = "off"
  )
  pushViewport(prof.ylab)
  grid.text(ylabel, 0.3, 0.5, gp = gpar(font = 1, cex = 0.7), rot = 90)
  popViewport()

  # Title
  tit <- viewport(
    x = 0, y = unit(1, "npc"), w = 1, h = unit(1, "lines"),
    just = c("left", "bottom"), clip = "off"
  )
  pushViewport(tit)
  title_text <- if (is.na(ftitle)) {
    paste0(fchr, ":", formatC(fstart, digits = 0, format = "f",
                               big.mark = ".", decimal.mark = ","),
           "..", formatC(fend, digits = 0, format = "f",
                         big.mark = ".", decimal.mark = ","))
  } else {
    ftitle
  }
  grid.text(title_text, 0.5, 0.5, gp = gpar(font = 1, cex = 0.7),
            default.units = "native")
  popViewport()

  popViewport()  # main

  # Capture and return grob
  grob <- grid.grab(wrap = TRUE)
  invisible(grob)
}
