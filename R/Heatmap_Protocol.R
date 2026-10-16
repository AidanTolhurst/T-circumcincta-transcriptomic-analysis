#Library calls
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

#Parameter definition
INPUT_FILE <- "C:/Users/tolhu/OneDrive/Documents/Cantacessi_Research/Data/merged_data.xlsx"  # output of the Python merge script
SHEET_NAME <- "Sheet1"
ID_COL     <- "gene_id"


STAGE_COMPARISONS <- c("L3_v_xL3", "xL3_v_L4", "L4_v_aMale", "L4_v_aFemale")
STAGE_LABELS       <- c("L3 vs xL3", "xL3 vs L4", "L4 vs aMale", "L4 vs aFemale")
LOG2FC_COLS        <- paste0("log2FoldChange_", STAGE_COMPARISONS)

AMP_PREDICTOR_COLS  <- c("ampir prob_AMP", "multipep antimicrobial", "multipep antibacterial")
PREDICTOR_LABELS    <- c("ampir\n(prob_AMP)", "multipep\n(antimicrobial)", "multipep\n(antibacterial)")
PREDICTOR_THRESHOLD <- 0  #Threshold for marking, 0 to allow for all retrieved values to be mapped

ROW_FILTER_THRESHOLD <- 0.9  # Threshold for gene inclusion

PLOT_TITLE  <- "Pairwise comparison of transcriptomics across life stages of T. circumcincta"

# Reproducible output filenames: controlled by environment variable REPRO_TAG or auto-generated
# tag uses UTC timestamp and, when available, the current git commit short SHA.
REPRO_TAG <- Sys.getenv("REPRO_TAG", unset = "")
.git_sha <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) character(0))
.git_sha <- if (length(.git_sha) == 0) NULL else .git_sha
if (REPRO_TAG == "") {
  .timestamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  .tag <- if (!is.null(.git_sha)) paste0(.timestamp, "_", .git_sha) else .timestamp
} else {
  .tag <- REPRO_TAG
}
OUTPUT_FILE <- paste0("Heatmap_rowfilter", ROW_FILTER_THRESHOLD, "_", .tag, ".pdf")

FONT_FAMILY <- "Helvetica"     

STANDARD_COLUMN_WIDTHS_IN <- c("1 column" = 3.5, "1.5 column" = 5, "2 column" = 7.2)
MAX_FIGURE_HEIGHT_IN <- 6   # MCP's hard cap on figure height - shared with Heatmap_Trajectory_Protocol.r

# This heatmap is one panel of a larger composite figure, not a standalone MCP submission
# figure - it gets resized (along with the other panels) when the composite is assembled,
# so hitting the column-width / 6 in limits on THIS script's own PDF doesn't matter, and
# every filtered gene is drawn at full cell size below. Set TRUE only if this script ever
# needs to produce a standalone, submission-ready figure on its own.
ENFORCE_MCP_PAGE_LIMITS <- FALSE

#Loading data
df <- read_excel(INPUT_FILE, sheet = SHEET_NAME)

cat("Loaded", nrow(df), "rows.\n")
cat("Columns:", paste(colnames(df), collapse = ", "), "\n\n")
cat("Stage labels: S1 =", STAGE_COMPARISONS[1], ", S2 =", STAGE_COMPARISONS[2],
    ", S3 =", STAGE_COMPARISONS[3], ", S4 =", STAGE_COMPARISONS[4], "\n\n")

df[AMP_PREDICTOR_COLS] <- lapply(df[AMP_PREDICTOR_COLS], function(x) ifelse(is.na(x), 0, x))


n_before_filter <- nrow(df)
keep <- apply(df[AMP_PREDICTOR_COLS] > ROW_FILTER_THRESHOLD, 1, any)
df <- df[keep, ]
cat("Filtered genes: kept", nrow(df), "of", n_before_filter,
    "(> ", ROW_FILTER_THRESHOLD, " in at least one of: ",
    paste(AMP_PREDICTOR_COLS, collapse = ", "), ")\n", sep = "")

# Order genes by their strongest AMP evidence across predictors, highest first
df$max_amp_score <- apply(df[AMP_PREDICTOR_COLS], 1, max, na.rm = TRUE)
df <- df %>% arrange(desc(max_amp_score))


ROW_HEIGHT_IN <- 0.15                        
header_margin_in <- 1.3

# MCP fixes a standalone figure at <=6 in tall and >=7pt type, so row count and cell size
# trade off directly there. That constraint only applies when ENFORCE_MCP_PAGE_LIMITS is
# on; otherwise every filtered gene is drawn at the full target row height, since this
# panel's overall size will be set later at composite-assembly time.
if (ENFORCE_MCP_PAGE_LIMITS) {
  max_rows_at_target <- floor((MAX_FIGURE_HEIGHT_IN - header_margin_in) / ROW_HEIGHT_IN)
  if (nrow(df) > max_rows_at_target) {
    cat(sprintf(
      "Note: %d genes pass the filter, but only %d fit at the target row height (%.2f in) within MCP's %.0f in max height. ",
      nrow(df), max_rows_at_target, ROW_HEIGHT_IN, MAX_FIGURE_HEIGHT_IN
    ))
    cat(sprintf("Keeping the top %d genes by AMP evidence - raise ROW_FILTER_THRESHOLD to control this list directly instead.\n",
                max_rows_at_target))
    df <- df[seq_len(max_rows_at_target), ]
  }
}
effective_row_height_in <- ROW_HEIGHT_IN
FIGURE_HEIGHT_IN <- max(nrow(df) * effective_row_height_in + header_margin_in, 3)

gene_ids <- df[[ID_COL]]

#Matrices
fc_vals <- as.matrix(df[LOG2FC_COLS])
rownames(fc_vals) <- gene_ids
colnames(fc_vals) <- STAGE_LABELS

#Colour
fc_range <- max(abs(fc_vals), na.rm = TRUE)
fc_col <- colorRamp2(c(-fc_range, 0, fc_range), c("#2166AC", "white", "#B2182B"))

# Build one 4-column (S1-S4) log2FC block per AMP predictor. A gene's log2FC is only
build_predictor_heatmap <- function(predictor_col, predictor_label, block_index) {
  flagged <- df[[predictor_col]] > PREDICTOR_THRESHOLD
  mat <- fc_vals
  mat[!flagged, ] <- NA
  
  legend_params <- list(
    title     = "log2FC",
    title_gp  = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    labels_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY)
  )
  
  Heatmap(
    mat,
    name = paste0("log2FC_", block_index),   # internal id only - legend title set above
    col = fc_col,
    na_col = "grey95",
    rect_gp = gpar(col = "grey50", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    height = unit(nrow(mat) * effective_row_height_in, "in"),
    show_row_names = (block_index == 1),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
    column_title = predictor_label,
    column_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
    column_names_rot = 45,   # was 90 - stage-comparison labels now sit at a 45-degree angle
    column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
    heatmap_legend_param = legend_params,
    show_heatmap_legend = (block_index == 1),  # one shared legend for all blocks
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.na(mat[i, j])) {
        grid.text("-", x, y, gp = gpar(fontsize = 8, col = "grey60", fontfamily = FONT_FAMILY))
      }
    },
    width = unit(length(STAGE_LABELS) * 0.9, "cm")
  )
}

ht_list <- NULL
for (k in seq_along(AMP_PREDICTOR_COLS)) {
  ht <- build_predictor_heatmap(AMP_PREDICTOR_COLS[k], PREDICTOR_LABELS[k], k)
  ht_list <- if (is.null(ht_list)) ht else ht_list + ht
}

# Width
heatmap_width_cm <- length(AMP_PREDICTOR_COLS) * length(STAGE_LABELS) * 0.9 + 2
legend_width_cm  <- 4  # legend sits vertically on the right
natural_width_in <- (heatmap_width_cm + legend_width_cm) / 2.54
if (ENFORCE_MCP_PAGE_LIMITS) {
  fits <- STANDARD_COLUMN_WIDTHS_IN >= natural_width_in
  FIGURE_WIDTH_IN <- if (any(fits)) min(STANDARD_COLUMN_WIDTHS_IN[fits]) else max(STANDARD_COLUMN_WIDTHS_IN)
} else {
  FIGURE_WIDTH_IN <- natural_width_in
}


pdf(OUTPUT_FILE, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN,
    family = FONT_FAMILY, useDingbats = FALSE)
draw(
  ht_list,
  ht_gap = unit(rep(3, length(AMP_PREDICTOR_COLS) - 1), "mm"),
  column_title = PLOT_TITLE,
  column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY),
  heatmap_legend_side = "right"
)
dev.off()

cat("\nSaved heatmap to", OUTPUT_FILE, "\n")


EXPORT_PNG_PREVIEW <- TRUE
if (EXPORT_PNG_PREVIEW) {
  png(sub("\\.pdf$", "_preview.png", OUTPUT_FILE),
      width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in",
      res = 1000, type = "cairo", family = FONT_FAMILY)
  draw(
    ht_list,
    ht_gap = unit(rep(3, length(AMP_PREDICTOR_COLS) - 1), "mm"),
    column_title = PLOT_TITLE,
    column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY),
    heatmap_legend_side = "right"
  )
  dev.off()
}
