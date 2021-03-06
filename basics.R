# Data
## Data package
dataDirectory <- file.path("EpigeneticsCSAMA2015","inst","bedfiles")

# Reading the filtered ChIP-seq reads
library(GenomicRanges)
library(rtracklayer)
library(IRanges)


input <- import.bed(file.path(dataDirectory, 'ES_input_filtered_ucsc_chr6.bed'))
rep1 <- import.bed(file.path(dataDirectory, 'H3K27ac_rep1_filtered_ucsc_chr6.bed'))
rep2 <- import.bed(file.path(dataDirectory, 'H3K27ac_rep2_filtered_ucsc_chr6.bed'))


# Preparation of the ChIP-seq and control samples: read extension
library(chipseq)

prepareChIPseq <- function(reads){
  frag.len <- median(estimate.mean.fraglen(reads))
  cat(paste0('Median fragment size for this library is ',round(frag.len)))
  reads.extended <- resize(reads, width = frag.len)
  return(trim(reads.extended))
}

input <- prepareChIPseq(input)
rep1 <- prepareChIPseq(rep1)
rep2 <- prepareChIPseq(rep2)

# Binning the ChIP-seq and control
## Generation of bins
load("~/R/chip-seq/EpigeneticsCSAMA2015/data/si.RData")
binsize = 200
bins = GenomicRanges::tileGenome(si['chr6'], tilewidth=binsize,
                  cut.last.tile.in.chrom=TRUE)
## Binning
BinChIPseq = function( reads, bins ){
  
  mcols(bins)$score = countOverlaps( bins, reads ) 
  return( bins ) 
}

input.200bins <- BinChIPseq( input, bins )
rep1.200bins <- BinChIPseq( rep1, bins )
rep2.200bins <- BinChIPseq( rep2, bins )

## Exporting binned data
export(input.200bins, 
       con='input_chr6.bedGraph',
       format = "bedGraph")
export(rep1.200bins, 
       con='H3K27ac_rep1_chr6.bedGraph',
       format = "bedGraph")
export(rep2.200bins, 
       con='H3K27ac_rep2_chr6.bedGraph',
       format = "bedGraph")

# Visualisation of ChIP-seq data with Gviz
library(Gviz)
load("~/R/chip-seq/EpigeneticsCSAMA2015/data/bm.RData")
AT = GenomeAxisTrack( )
plotTracks(c( bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto", 
           cex.title=1, fontsize=10 )

input.track = DataTrack(input.200bins, 
                        strand="*", genome="mm9", col.histogram='gray',
                        fill.histogram='black', name="Input", col.axis="black",
                        cex.axis=0.4, ylim=c(0,150))

rep1.track = DataTrack(rep1.200bins, 
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 1", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))

rep2.track = DataTrack(rep2.200bins, 
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 2", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))

plotTracks(c(input.track, rep1.track, rep2.track, bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto", 
           type="histogram", cex.title=0.7, fontsize=10 )

# ChIP-seq peaks
## Peaks - basic analysis in R
peaks.rep1 = import.bed(file.path(dataDirectory,'Rep1_peaks_ucsc_chr6.bed'))
peaks.rep2 = import.bed(file.path(dataDirectory,'Rep2_peaks_ucsc_chr6.bed'))
peaks1.track = AnnotationTrack(peaks.rep1, 
                               genome="mm9", name='Peaks Rep. 1',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)
peaks2.track = AnnotationTrack(peaks.rep2, 
                               genome="mm9", name='Peaks Rep. 2',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)

plotTracks(c(input.track, rep1.track, peaks1.track,
             rep2.track, peaks2.track, bm, AT),
           from=122630000, to=122700000,
           transcriptAnnotation="symbol", window="auto", 
           type="histogram", cex.title=0.7, fontsize=10 )

## Venn diagrams
ovlp = findOverlaps( peaks.rep1, peaks.rep2 )
ov = min( length(unique( queryHits(ovlp) )), length(unique( subjectHits(ovlp) ) ) )
library(VennDiagram)

draw.pairwise.venn( 
  area1=length(peaks.rep1),
  area2=length(peaks.rep2), 
  cross.area=ov, 
  category=c("rep1", "rep2"), 
  fill=c("steelblue", "blue3"), 
  cat.cex=0.7)

enriched.regions = Reduce(subsetByOverlaps, list(peaks.rep1, peaks.rep2))
enr.reg.track = AnnotationTrack(enriched.regions,
                                genome="mm9", name='Enriched regions',
                                chromosome='chr6',
                                shape='box',fill='green3',size=2)

plotTracks(c(input.track, rep1.track, peaks1.track,
             rep2.track, peaks2.track, enr.reg.track, 
             bm, AT),
           from=122630000, to=122700000,
           transcriptAnnotation="symbol", window="auto", 
           type="histogram", cex.title=0.5, fontsize=10 )

## Isolation of promoters overlapping H3K27ac peaks
### Identification of promoters
load("~/R/chip-seq/EpigeneticsCSAMA2015/data/egs.RData")
egs$TSS = ifelse( egs$strand == "1", egs$start_position, egs$end_position )
promoter_regions = 
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = egs$TSS - 200,
                            end = egs$TSS + 200 ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$external_gene_id)
### Overlapping promoters with H3K27ac enriched regions
ovlp2 = findOverlaps( enriched.regions, promoter_regions )

cat(sprintf( "%d of %d promoters are overlapped by an enriched region.",
             length( unique(subjectHits(ovlp2)) ), length( promoter_regions ) ) )

ovlp2b = findOverlaps( promoter_regions, enriched.regions )

cat(sprintf( "%d of %d enriched regions overlap a promoter.",
             length( unique( subjectHits(ovlp2b) ) ), length( enriched.regions ) ) )

promotor_total_length = sum(width(reduce(promoter_regions)))
promotor_fraction_of_chromosome_6 = promotor_total_length / seqlengths(si)["chr6"]
binom.test( length( unique( subjectHits( ovlp2b ) ) ), length( enriched.regions ), promotor_fraction_of_chromosome_6 )
pos.TSS = egs[ unique( queryHits( findOverlaps( promoter_regions, enriched.regions ) ) ),]

## Analysis of the distribution of H3K27ac around a subset of gene promoters
tiles = sapply( 1:nrow(pos.TSS), function(i)
  if( pos.TSS$strand[i] == "1" )
    pos.TSS$TSS[i] + seq( -1000, 900, length.out=20 )
  else
    pos.TSS$TSS[i] + seq( 900, -1000, length.out=20 ) )

tiles = GRanges(tilename = paste( rep( pos.TSS$ensembl_gene_id, each=20), 1:20, sep="_" ),
                seqnames = Rle( rep(paste0('chr', pos.TSS$chromosome_name), each=20) ), 
                ranges = IRanges(start = as.vector(tiles),
                                 width = 100),
                strand = Rle(rep("*", length(as.vector(tiles)))),
                seqinfo=si)

H3K27ac.p = countOverlaps( tiles, rep1) +
  countOverlaps( tiles, rep2 )

H3K27ac.p.matrix = matrix( H3K27ac.p, nrow=nrow(pos.TSS), 
                           ncol=20, byrow=TRUE )

colors = colorRampPalette(c('white','red','gray','black'))(100) 

layout(mat=matrix(c(1,2,0,3), 2, 2), 
       widths=c(2,2,2), 
       heights=c(0.5,5,0.5,5), TRUE)


par(mar=c(2,2,0.75,0.5))
image(seq(0, max(H3K27ac.p.matrix), length.out=100), 1,
      matrix(seq(0, max(H3K27ac.p.matrix), length.out=100),100,1),
      col = colors,
      xlab='Distance from TSS', ylab='',
      main='Number of reads', yaxt='n',
      lwd=3, axes=TRUE)
box(col='black', lwd=2)
image(x=seq(-1000, 1000, length.out=20),
      y=1:nrow(H3K27ac.p.matrix),
      z=t(H3K27ac.p.matrix[order(rowSums(H3K27ac.p.matrix)),]), 
      col=colors,
      xlab='Distance from TSS (bp)',
      ylab='Promoters', lwd=2)
box(col='black', lwd=2)
abline(v=0, lwd=1, col='gray')
plot(x=seq(-1000, 1000, length.out=20),
     y=colMeans(H3K27ac.p.matrix),
     ty='b', pch=19,
     col='red4',lwd=2,
     ylab='Mean tag count',
     xlab='Distance from TSS (bp)')
abline(h=seq(1,100,by=5),
       v=seq(-1000, 1000, length.out=20),
       lwd=0.25, col='gray')
box(col='black', lwd=2)