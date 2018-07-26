

#' Plot coverages along genome
#'
#' @param fstart start of the genomic window (bp)
#' @param fend   end of the genomic window (bp)
#' @param fchr   chromosome of the genomic window
#' @param profs  list of coverages to be plotted
#' @param cols colors of the profiles
#' @param ann a dataframe of annotations
#' @param ylable the label of the y axis
#' @param ylims list of ylimits for plotting the coverages. list(c(min, max), c(min, max))
#' @param txdb a transcription database used for plotting gene models
#' @param ftitle title of the plot
#' @param collapse collapse gene models (default = TRUE)
#' @param with.genes.highlited vector of gene ids that should be highlighted
#' @param plot.labels plot labels lf gene (default = TRUE)
#' @param grid plot grid (default = TRUE)
#' @param with.average plot average (default = FALSE)
#'
#' @export

plotProfiles <-
  function(fstart, fend, fchr, profs, cols = c(), ann=NULL, ylabel = "coverage", ylims = list(), txdb, ftitle = NA,
           collapse = T, with.genes.highlited = c(), plot.labels = T, grid = F, with.average = F) {

    require(grid)
    require(IRanges)
    require(HilbertVis)
    require(RColorBrewer)
    require(AnnotationDbi)

    ########################################################################################################################################################################################################
    # test args
    ########################################################################################################################################################################################################

    if ( !(fchr %in% unlist(lapply(profs, names))) )
      stop("chromosome provided in parameter 'fchr' does not match chromosome names in coverages")
    if ( !(fchr %in% levels(seqnames(genes(txdb)) )))
       warning("chromosome provided in parameter 'fchr' does not match chromosome names in annotation object")

    options(scipen = 100)

    ########################################################################################################################################################################################################
    # FUNCTION DEFINITIONS
    ########################################################################################################################################################################################################
    profile.xscale <- function(xscale) {
      #grid.rect(gp = gpar(col = "grey"))
      main.ticks <- grid.pretty(xscale)
      dist.main <- main.ticks[2] - main.ticks[1]
      minor.ticks <- seq(
        from = main.ticks[1] - dist.main,
        to = main.ticks[length(main.ticks)] + dist.main,
        by = dist.main / 5
      )
      minor.ticks <- minor.ticks[!(minor.ticks %in% main.ticks)]
      minor.ticks <-
        minor.ticks[minor.ticks > xscale[1] & minor.ticks < xscale[2]]
      for (i in minor.ticks) {
        grid.lines(c(i, i),
                   c(1, 0.9),
                   default.units = "native",
                   gp = gpar(col = 1))
      }
      for (i in main.ticks) {
        grid.lines(c(i, i),
                   c(1, 0.8),
                   default.units = "native",
                   gp = gpar(col = 1))
        grid.text(
          paste(
            formatC(
              i / 1000,
              format = "f",
              digits = 0,
              big.mark = ".",
              decimal.mark = ","
            ),
            "K"
          ),
          i,
          0.7,
          default.units = "native",
          just = c("center", "top"),
          gp = gpar(cex = 0.4)
        )
      }
      grid.lines(xscale, 1, default.units = "native", gp = gpar(col = 1))
      #grid.xaxis(name="scale", gp=gpar(cex=0.5))

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
      minor.ticks <-
        minor.ticks[minor.ticks > xscale[1] & minor.ticks < xscale[2]]
      for (i in minor.ticks) {
        grid.lines(
          c(i, i),
          c(0, 1),
          default.units = "native",
          gp = gpar(col = "lightgrey", lwd = 0.5)
        )
      }
      for (i in main.ticks) {
        grid.lines(
          c(i, i),
          c(0, 1),
          default.units = "native",
          gp = gpar(col = "lightgrey", lwd = 0.5)
        )
      }
    }
    gff.get.parent <- function(x) {
      rx <- regexpr("Parent=([^;]*)", text = x, perl = T)
      substring(
        x,
        attr(rx, "capture.start"),
        attr(rx, "capture.start") + attr(rx, "capture.length") - 1
      )
      #	res <- substring(strsplit(x,';')[[1]][3],8,1000)
      #	if (is.na(res)) {res<-""}
      #	res
    }
    ts.grid.rect <- function(start, end,  bottom, top, col = 1, lcol = 1, lwd = 1) {
      grid.rect(
        x = start,
        y = bottom,
        width = end - start,
        height = top - bottom,
        just = c("left", "bottom"),
        default.units = "native",
        gp = gpar(
          fill = col,
          col = lcol,
          lwd = lwd
        )
      )
    }
    filter.longest.txt <- function(txts) {
      txts[which.max(width(txts))]
    }
    plot.points <- function(layout.frame, profile.frame, col = 1) {
      grid.points(
        layout.frame$pos,
        profile.frame,
        pch = 16,
        gp = gpar(col = col, cex = 0.1)
      )
    }
    plot.bars <- function(design.frame, profile.frame, fcol = 1, bcol = 0) {
      for (i in 1:nrow(design.frame)) {
        grid.rect(
          x = design.frame$pos[i],
          y = 0,
          width = design.frame$length[i],
          height = profile.frame[i],
          just = c("left", "bottom"),
          default.units = "native",
          gp = gpar(fill = bcol, col = fcol)
        )
      }
    }
    plot.bars2 <- function(design.frame, profile.frame, y.min, fcol = 1, bcol = 0) {
      for (i in 1:nrow(design.frame)) {
        grid.rect(
          x = design.frame$pos[i],
          y = y.min,
          width = design.frame$length[i],
          height = profile.frame[i] - y.min,
          just = c("left", "bottom"),
          default.units = "native",
          gp = gpar(fill = bcol, col = fcol)
        )
      }
    }
    averageWindow <- function(pos, design, profile, window.size = 500) {
      window <-
        design$pos + design$length >= pos - window.size &
        design$pos <= pos + window.size
      median(profile[window])
    }
    plot.genes <- function(vp, strand, frame.genes, bumps, gene.height, label.shift,
                           collapse,  with.genes.highlited = c(), plot.labels) {
      #browser()
      exon <- 2.0
      gene.frame <- frame.genes[strand(frame.genes) == strand]
      pushViewport(vp)
      #grid.rect(gp = gpar(col = "grey"))
      if (length(gene.frame) > 0) {
        for (i in 1:length(gene.frame)) {
          txts.names <-
            suppressMessages(AnnotationDbi::select(
              txdb,
              keys = names(gene.frame)[i],
              columns = "TXNAME",
              keytype = "GENEID"
            )$TXNAME)
          txts <- transcripts(txdb, filter = list(tx_name = txts.names))
          if (length(txts) > 0) {
            if (collapse) {
              filter.longest.txt <- function(txts) {
                txts[which.max(width(txts))]
              }

              txts <- filter.longest.txt(txts)
            }
            for (k in 1:length(txts)) {
              gene.id <-
                suppressMessages(
                  AnnotationDbi::select(
                    txdb,
                    keys = elementMetadata(txts[, 2])[k, 1],
                    columns = "GENEID",
                    keytype = "TXNAME"
                  )$GENEID
                )
              col <-
                ifelse(gene.id %in% with.genes.highlited,
                       "#FFAAAA",
                       "#F2F2F2")
              lcol <-
                ifelse(gene.id %in% with.genes.highlited,
                       "#FF0000",
                       "#000000")

              if (strand == "-") {
                t1  <-
                  (bumps * 10) - ((k - 1) * gene.height + gene.height / 2 - label.shift /
                                    2)
              }  else {
                t1  <- ((k - 1) * gene.height + gene.height / 2) - (label.shift / 2)
              } #fw
              ts.grid.rect(start(txts[k]),
                           end(txts[k]),
                           t1,
                           t1,
                           lcol = lcol,
                           lwd = 2)

              exons <-
                exons(txdb, filter = list(tx_name = elementMetadata(txts[, 2])[1, 1]))

              if (length(exons) > 0) {
                for (l in 1:length(exons)) {
                  ts.grid.rect(
                    start(exons[l]),
                    end(exons[l]),
                    t1 - exon,
                    t1 + exon,
                    col = col,
                    lcol = lcol,
                    lwd = 1
                  )
                }
              }
              if (plot.labels) {
                label <- gene.id
                if (strand == "-") {
                  grid.text(
                    label,
                    end(txts)[k],
                    t1 - label.shift - 0.5,
                    default.units = "native",
                    just = c("right", "bottom"),
                    gp = gpar(cex = 0.5)
                  )
                }
                else {
                  grid.text(
                    label,
                    end(txts)[k],
                    t1 + label.shift + 0.4,
                    default.units = "native",
                    just = c("right", "top"),
                    gp = gpar(cex = 0.5)
                  )
                }
              }
            }
          }
        }
      }
      popViewport()
    }
    plot.profiles <- function(profs, fchr, fstart, fend, ylims, cols,
                              xscale, with.average = F) {
      vsize <- as.integer(dev.size()[1] * 150)

      #margins <- unit(0.04, "npc")
      margins <- unit(0.19, "lines")
      font.size.label <- 0.7
      panel.background <- "#cccccc25"

      if (length(cols) == 0) {
        cols <-
          RColorBrewer::brewer.pal((ifelse(length(profs) > 2, length(profs), 3)) , "Dark2")
      }

      for (i in 1:length(profs)) {
        if (class(profs[[i]]) == "SimpleRleList") {
          #			framel <- window(profs[[i]][[fchr]], start= fstart,end= fend)
          #			xl <- c(0,start(framel),max(end(framel))) + fstart
          #			yl <- c(0,runValue(framel),0)

          xl <- c(fstart, seq(fstart, fend, length.out = vsize), fend)
          yl <-
            c(0, HilbertVis::shrinkVector(as.vector(profs[[i]][[fchr]])[fstart:fend], newLength =
                                vsize), 0)
          yl[is.na(yl)] <- 0
          #browser()
          if (length(ylims) == 0) {
            vmax <-
              max(pretty(0:ceiling(max(
                yl, na.rm = T
              ))))
          } else {
            vmax <- ylims[[i]][2]
          }

          #clipped
          prof1 <-
            viewport(
              x = 0,
              y = unit((length(profs) - i) * 1 / length(profs), "npc") + 2 * margins,
              w = 1,
              h = unit(1 / length(profs), "npc") - 2 * margins,
              just = c("left", "bottom"),
              xscale = xscale,
              yscale = c(0, vmax),
              clip = "on"
            )
          pushViewport(prof1)
          grid.rect(gp = gpar(fill = panel.background, lty = 0))
          #grid.yaxis(gp=gpar(cex=0.5))
          #grid.lines(xscale, 0, default.units="native",gp=gpar(col="grey66", lwd=0.5))
          #grid.polygon(x=xl, y=yl,default.units="native",gp=gpar(col="#33333370", fill=paste(cols[[i]],"95", sep="")))
          grid.lines(xscale,
                     0,
                     default.units = "native",
                     gp = gpar(col = "grey66", lwd = 0.5))
          grid.polygon(
            x = xl,
            y = yl,
            default.units = "native",
            gp = gpar(
              col = "#55555570",
              lwd = 0.5,
              fill = NA
            )
          )
          grid.polygon(
            x = xl,
            y = yl,
            default.units = "native",
            gp = gpar(col = NA, fill = cols[[i]])
          )
          grid.text(
            names(profs)[[i]],
            unit(0.4, "lines"),
            unit(1, "npc") - unit(0.3, "lines"),
            just = c("left", "top"),
            gp = gpar(cex = font.size.label, font = 1)
          )
          popViewport()

          ##unclipped
          prof1 <-
            viewport(
              x = 0,
              y = unit((length(profs) - i) * 1 / length(profs), "npc") + 2 * margins,
              w = 1,
              h = unit(1 / length(profs), "npc") - 2 * margins,
              just = c("left", "bottom"),
              xscale = xscale,
              yscale = c(0, vmax),
              clip = "off"
            )
          pushViewport(prof1)
          grid.rect(gp = gpar(lwd = 1, fill = NA))
          popViewport()

          yax <-
            viewport(
              x = unit(-0.2, "lines"),
              y = unit((length(profs) - i) * 1 / length(profs), "npc") + 2 * margins,
              w = unit(0.2, "lines"),
              h = unit(1 / length(profs), "npc") - 2 * margins,
              just = c("left", "bottom"),
              xscale = xscale,
              yscale = c(0, vmax),
              clip = "off"
            )
          pushViewport(yax)
          #			at <- pretty(seq(0, floor(vmax)), n=5)[c(1,3,5)]
          grid.yaxis(gp = gpar(cex = 0.45))#, at=at)
          #grid.text("label", unit(-2.5, "lines") , 0.5, just="center", rot=90, gp=gpar(cex=font.size.label, font=2))
          popViewport()


        }
      }

    }
    plot.annotation <- function(ann, fchr, fstart, fend, xscale) {

      margins <- unit(0.05, "lines")
      #grid.rect(gp = gpar(col = "#333333"))

      for (i in 1:length(ann)) {
        ann1 <- viewport(x = 0, y = unit((length(ann)-i) * 1/length(ann), "npc") + 2*margins,
                         w = 1, h = unit(1/length(ann), "npc")- 2*margins,
                         just = c("left", "bottom"), xscale=xscale, clip="on",
                         yscale=c(0, 1))
        pushViewport(ann1)

        ca <- ann[[i]]
        ca <- ca[ca$chr==fchr & ca$end>fstart & ca$start<fend,]
        #browser()
        if (nrow(ca) > 0) {
          for (k in 1:nrow(ca)) {
            ts.grid.rect(ca$start[k], ca$end[k], 0, 1, col=ifelse(is.null(ca$col[k]),1,as.character(ca$col[k])),
                         lcol=ifelse(is.null(ca$col[k]),1,as.character(ca$col[k])), lwd=1)
            grid.text(as.character(ca$label[k]), unit(ca$start[k],"native")-unit(0.003,"npc"), unit(0.5, "npc"), just=c("right", "center"), gp=gpar(cex=0.5, font=1))

          }
        }

        popViewport()

        ann1 <- viewport(x = 0, y = unit((length(ann)-i) * 1/length(ann), "npc") + 2*margins,
                         w = 1, h = unit(1/length(ann), "npc")- 2*margins,
                         just = c("left", "bottom"), xscale=xscale, clip="off",
                         yscale=c(0, 1))
        pushViewport(ann1)
        grid.text(names(ann)[i], unit(-0.01, "npc"), unit(0.5, "npc"), just=c("right", "center"), gp=gpar(cex=0.5, font=1))

        popViewport()


      }


    }

    ########################################################################################################################################################################################################
    # VAL DEFINITIONS
    ########################################################################################################################################################################################################
    ## DEFINES
    gene.height <- 10
    gene.lines <- 1.3
    exon <- 2.4
    label.shift <- 6
    scale.y.offset <- 1.2
    vsize <- as.integer(dev.size()[1] * 150)
    xscale <- c(fstart, fend)

    seqlevels(txdb) <- fchr
    ### START PLOT

    grid.newpage()
    main <-
      viewport(
        x = unit(3, "lines"),
        y = unit(1, "lines"),
        width = unit(1, "npc") - unit(4, "lines"),
        height = unit(1, "npc") - unit(3, "lines"),
        just = c("left", "bottom"),
        xscale = xscale
      )
    pushViewport(main)

    ########
    # bkgnd grid
    ########

    if (grid) {
      bg.grid(xscale)
    }

    frameRange <- GRanges(fchr, IRanges(fstart, fend))

    all.genes <- genes(txdb)
    frame.genes <- subsetByOverlaps(genes(txdb), frameRange)
    frame.txts <- subsetByOverlaps(transcripts(txdb), frameRange)
    frame.genes.rev <- frame.genes[strand(frame.genes) == "-"]
    frame.genes.fwd <- frame.genes[strand(frame.genes) == "+"]

    ## how many bumps? one gene

    #gff.frame <- get.allfeatures.in.frame(gff, fchr, fstart, fend)
    #gene.frame.rev <- get.namedOrientation(get.namedFeatures(gff.frame, "gene"), "-")
    #gene.frame.fwd <- get.namedOrientation(get.namedFeatures(gff.frame, "gene"), "+")

    ## max number of transcript isoforms
    get.no.tracks <- function(txdb, genes) {
      max(sapply(names(genes), function(x) {
        nrow(suppressMessages(
          AnnotationDbi::select(
            txdb,
            keys = x,
            columns = "TXNAME",
            keytype = "GENEID"
          )
        ))
      }))
    }

    rbumps <- ifelse(collapse, 1, get.no.tracks(txdb, frame.genes.rev))
    fbumps <- ifelse(collapse, 1, get.no.tracks(txdb, frame.genes.fwd))

    ########
    #xscale
    ########
    scale <-
      viewport(
        x = 0,
        y = unit(gene.lines * rbumps + 0.5, "lines") - unit(scale.y.offset, "lines"),
        w = 1,
        h = unit(1.5, "lines"),
        clip = "off",
        just = c("left", "bottom"),
        name = "scale",
        xscale = xscale
      )
    pushViewport(scale)
    profile.xscale(xscale)
    popViewport()

    ########
    # genes track
    ########

    rev <-
      viewport(
        x = 0,
        y = 0,
        w = 1,
        h = unit(gene.lines * rbumps, "lines"),
        clip = "on",
        just = c("left", "bottom"),
        name = "rev",
        xscale = xscale,
        yscale = c(-1, gene.height * rbumps + 1)
      )


    plot.genes(
      vp = rev,
      strand = "-",
      frame.genes,
      bumps = rbumps,
      gene.height,
      label.shift,
      collapse,
      with.genes.highlited = c(),
      plot.labels
    )


    ## switch annot
    if (is.null(ann)) sval=0.5 else sval=0.7
    fwd <-
      viewport(
        x = 0,
        y = unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines"),
        w = 1,
        h = unit(gene.lines * fbumps, "lines"),
        clip = "on",
        just = c("left", "bottom"),
        name = "fwd",
        xscale = xscale,
        yscale = c(-1, gene.height * fbumps + 1)
      )
    plot.genes(
      vp = fwd,
      strand = "+",
      frame.genes,
      bumps = fbumps,
      gene.height,
      label.shift,
      collapse,
      with.genes.highlited = c(),
      plot.labels
    )

        ## switch annot
    if (!is.null(ann)) {
      ########
      #annot
      ########
      annot <- viewport(x = 0, y = unit(gene.lines*rbumps+0.5, "lines") + unit(0.5, "lines") +
                          unit(gene.lines*fbumps, "lines"), w = 1,
                        h = unit(0.5, "lines"),
                        just = c("left", "bottom"), name = "annot", xscale=xscale)

      pushViewport(annot)
      plot.annotation(ann, fchr, fstart, fend, xscale)
      popViewport()

    }

    if (is.null(ann)) sval=0.2 else sval=0.7
     ########
    #prof
    ########
    prof <-
      viewport(
        x = 0,
        y = unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines") +
          unit(gene.lines * fbumps, "lines"),
        w = 1,
        h = unit(1, "npc") - (
          unit(gene.lines * rbumps + 0.5, "lines") + unit(sval, "lines") +
            unit(gene.lines * fbumps, "lines")
        ),
        just = c("left", "bottom"),
        name = "fwd",
        xscale = xscale
      )


    pushViewport(prof)
    #grid.rect(gp = gpar(col = "grey"))
    plot.profiles(profs, fchr, fstart, fend, ylims, cols, xscale, with.average)
    popViewport() # profiles

    prof.ylab <-
      viewport(
        x = unit(-1.8, "lines"),
        y = unit(gene.lines * rbumps + 0.5, "lines") + unit(0.7, "lines") +
          unit(gene.lines * fbumps, "lines") + unit(0.2, "lines"),
        w = unit(0.4, "lines"),
        h = unit(1, "npc") - (
          unit(gene.lines * rbumps + 0.5, "lines") + unit(0.7, "lines") +
            unit(gene.lines * fbumps, "lines") + unit(0.2, "lines")
        ),
        just = c("left", "bottom"),
        clip = "off"
      )


    pushViewport(prof.ylab)
    grid.text(ylabel,
              0.3,
              0.5,
              gp = gpar(font = 1, cex = 0.7),
              rot = 90)
    popViewport() # profiles


    ########
    #tit
    ########
    tit <-
      viewport(
        x = 0,
        y = unit(1, "npc"),
        w = 1,
        h = unit(1, "lines"),
        just = c("left", "bottom"),
        clip = "off"
      )
    pushViewport(tit)
    #grid.rect(gp = gpar(col = "grey"))
    grid.text(
      ifelse(
        is.na(ftitle),
        paste(
          fchr,
          ":",
          formatC(
            fstart,
            digits = 0,
            format = "f",
            big.mark = ".",
            decimal.mark = ","
          ),
          "..",
          formatC(
            fend,
            digits = 0,
            format = "f",
            big.mark = ".",
            decimal.mark = ","
          ),
          sep = ""
        ),
        ftitle
      ),
      0.5,
      0.5,
      gp = gpar(font = 1, cex = 0.7),
      default.units = "native"
    )
    popViewport()

    popViewport() # main
    #main
  }


