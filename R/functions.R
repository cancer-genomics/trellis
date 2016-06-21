myGrob <- function(x, r=unit(0.25, "inches")){
  tg <- textGrob(x)
  cg <- circleGrob(r=r)
  boxedText <- gTree(children=gList(tg, cg))
}

plotGrob <- function(g, x, y, fill="gray80", color="steelblue"){
  grid.draw(editGrob(g, vp=viewport(x, y), gp=gpar(fill=fill)))
  grid.draw(editGrob(g, vp=viewport(x, y),
                     gp=gpar(fill="transparent", color=color)))
}

ar <- arrow(ends="last", length=unit(0.075, "inches"), type="closed")

addArrows <- function(x, x2, y, y2, h=unit(0.3, "inches")){
  grid.move.to(unit(x[1], "npc") + h, unit(y[1], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[1], "npc"), arrow=ar,
               gp=gpar(fill="black"))

  ## PR -> D
  grid.move.to(unit(x[1], "npc") + h, unit(y[2], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[1], "npc") - 1/3*h, arrow=ar,
               gp=gpar(fill="black"))

  ## RD -> A
  grid.move.to(unit(x[1], "npc") + h, unit(y[1], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[2], "npc") + 1/5*h, arrow=ar,
               gp=gpar(fill="black"))

  ## PR -> A
  grid.move.to(unit(x[1], "npc") + h, unit(y[2], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[2], "npc") - 1/5*h, arrow=ar,
               gp=gpar(fill="black"))

  ## PR -> I
  grid.move.to(unit(x[1], "npc") + h, unit(y[2], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[3], "npc") + 1/5*h, arrow=ar,
               gp=gpar(fill="black"))

  ## PR -> T
  grid.move.to(unit(x[1], "npc") + h, unit(y[2], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[4], "npc") + 1/5*h, arrow=ar,
               gp=gpar(fill="black"))

  ## SR -> T
  grid.move.to(unit(x[1], "npc") + h, unit(y[3], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[4], "npc") - 1/5*h, arrow=ar,
               gp=gpar(fill="black"))

  ## SR -> I
  grid.move.to(unit(x[1], "npc") + h, unit(y[3], "npc"))
  grid.line.to(unit(x2[1], "npc") - h, unit(y2[3], "npc") - 1/5*h, arrow=ar,
               gp=gpar(fill="black"))
}

integrationFig <- function(x, y, x2, y2, grobs){
  plotGrob(grobs[["RD"]], x=x[1], y=y[1])
  plotGrob(grobs[["PR"]], x=x[1], y=y[2])
  plotGrob(grobs[["SR"]], x=x[1], y=y[3])
  plotGrob(grobs[["D"]], x=x2[1], y=y2[1], fill="lightblue")
  plotGrob(grobs[["A"]], x=x2[1], y=y2[2], fill="salmon")
  plotGrob(grobs[["I"]], x=x2[1], y=y2[3], fill="orange")
  plotGrob(grobs[["T"]], x=x2[1], y=y2[4], fill="beige")
  addArrows(x=x, x2=x2, y=y, y2=y2)
}

exampleRpData <- function(bviews, path){
  bv <- bviews[, "CGOV33T.bam"]
  id <- colnames(bv)
  sv <- readRDS(file.path(path, "data/segment/1deletions/CGOV33T.bam.rds"))
  region <- GRanges("chr1", IRanges(82350000, 82410000))
  i <- subjectHits(findOverlaps(region, variant(sv)))
  sv2 <- sv[i]
  prp <- proper(sv2)
  irp <- improper(sv2)
  frp <- c(granges(first(prp)), granges(first(irp)))
  lrp <- c(granges(last(prp)), granges(last(irp)))
  dat <- data.frame(start=start(frp), end=end(frp),
                    start2=start(lrp), end2=end(lrp))
  dat <- dat[order(dat$start), ]
  dat$y <- seq_len(nrow(dat))
  dat
}

methodsGrobs <- function(){
  grobRD <- myGrob("RD")
  grobPR <- myGrob("PR")
  grobSR <- myGrob("SR")
  grobD <- myGrob("D")
  grobA <- myGrob("A")
  grobI <- myGrob("I")
  grobT <- myGrob("T")
  grobs <- list(RD=grobRD, PR=grobPR,
              SR=grobSR, D=grobD,
              A=grobA, I=grobI,
              T=grobT)
}

rdData <- function(path, bview){
  region <- GRanges("chr2", IRanges(215e6, 235e6))
  pview <- PreprocessViews2(bview)
  paths(pview) <- file.path(path, "data/preprocess/3background_adj/CGOV33T.bam.rds")
  setScale(pview) <- 1000
  index <- subjectHits(findOverlaps(region, rowRanges(pview)))
  pview <- pview[index, ]
  segs <- readRDS(file.path(path, "data/segment/0cbs/CGOV33T.bam.rds"))
  segs <- subsetByOverlaps(segs, region)
  df.segs <- data.frame(x=start(segs)/1e6, xend=end(segs)/1e6,
                        y=segs$seg.mean)
  df <- data.frame(r=assays(pview)[,1], start=start(rowRanges(pview)))
  df <- df[seq(1, nrow(df), 20), ]
  med <- median(df$r, na.rm=TRUE)
  df$r <- df$r - med
  df.segs$y <- df.segs$y - med
  df.segs$xend[5] <-  end(region)/1e6
  df.segs$x[1] <-  start(region)/1e6
  list(data=df, segs=df.segs, region=region)
}
