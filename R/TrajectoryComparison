#Library calls
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

INPUT_FILE <- file.path("data", "merged_data.xlsx")
SHEET_NAME <- "Sheet1"
ID_COL     <- "gene_id"

# AMP-specific predictor scores used to decide which genes are plotted
AMP_PREDICTOR_COLS  <- c("ampir prob_AMP", "multipep antimicrobial", "multipep antibacterial")
ROW_FILTER_THRESHOLD <- ROW_FILTER_THRESHOLD
STAGE_ORDER  <- c("L3", "xL3", "L4", "aMale", "aFemale")
Z_SCORE_ROWS <- TRUE

PLOT_TITLE  <- "Relative expression trajectory of\n AMP-predicted peptides across life stages"

FONT_FAMILY <- "Helvetica"                                       
STANDARD_COLUMN_WIDTHS_IN <- c("1 column" = 3.5, "1.5 column" = 5, "2 column" = 7.2)
MAX_FIGURE_HEIGHT_IN <- 6            
CELL_WIDTH_CM    <- 0.9              
ROW_HEIGHT_IN    <- 0.15           
header_margin_in <- 1.3 + ifelse(grepl("\n", PLOT_TITLE), 0.3, 0)    


#Loading data
df <- read_excel(INPUT_FILE, sheet = SHEET_NAME)
cat("Loaded", nrow(df), "rows.\n")
cat("Columns:", paste(colnames(df), collapse = ", "), "\n\n")

df[AMP_PREDICTOR_COLS] <- lapply(df[AMP_PREDICTOR_COLS], function(x) ifelse(is.na(x), 0, x))

n_before <- nrow(df)
keep <- apply(df[AMP_PREDICTOR_COLS] > ROW_FILTER_THRESHOLD, 1, any)
df <- df[keep, ]
cat("Filtered genes: kept", nrow(df), "of", n_before,
    "(> ", ROW_FILTER_THRESHOLD, " in at least one of: ",
    paste(AMP_PREDICTOR_COLS, collapse = ", "), ")\n", sep = "")

gene_ids <- df[[ID_COL]]

#Build the chained relative-expression matrix: rows = genes, columns = stages
expr <- matrix(NA_real_, nrow = nrow(df), ncol = length(STAGE_ORDER),
               dimnames = list(gene_ids, STAGE_ORDER))
expr[, "L3"] <- 0  

fc <- function(comparison) {
  col <- paste0("log2FoldChange_", comparison)
  df[[col]] * SIGN_FLIP[[comparison]]
}

expr[, "xL3"]     <- expr[, "L3"] - fc("L3_v_xL3")
expr[, "L4"]      <- expr[, "xL3"] - fc("xL3_v_L4")
expr[, "aMale"]   <- expr[, "L4"] - fc("L4_v_aMale")
expr[, "aFemale"] <- expr[, "L4"] - fc("L4_v_aFemale")

n_incomplete <- sum(!complete.cases(expr))
if (n_incomplete > 0) {
  cat(n_incomplete, "of", nrow(expr),
      "genes are missing a value somewhere along the chain and are dropped ",
      "(", nrow(expr) - n_incomplete, " remain with a complete trajectory).\n", sep = "")
  df   <- df[complete.cases(expr), , drop = FALSE]
  expr <- expr[complete.cases(expr), , drop = FALSE]
}

if (nrow(expr) == 0) {
  stop("No genes have a complete trajectory across all 4 comparisons - ",
       "the pairwise sheets may not share enough genes in common for this chaining approach.")
}

#Order genes by their strongest AMP evidence, highest first
max_amp_score <- apply(df[AMP_PREDICTOR_COLS], 1, max, na.rm = TRUE)
expr <- expr[order(-max_amp_score), , drop = FALSE]

#Optionally z-score each row
plot_mat <- if (Z_SCORE_ROWS) t(scale(t(expr))) else expr

# ComplexHeatmap title
legend_title <- if (Z_SCORE_ROWS) "Relative change in\nlog2FC (z-scored)" else "Relative change in\nlog2FC"

# A gene with zero variance 
zero_var_rows <- apply(plot_mat, 1, function(r) all(is.na(r)))
if (any(zero_var_rows)) {
  cat(sum(zero_var_rows), "gene(s) show no variance across stages and will be shown as blank cells.\n")
}

#Colour: diverging, centred at 0
val_range <- max(abs(plot_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-val_range, 0, val_range), c("#2166AC", "white", "#B2182B"))

# Figure Sizing
if (ENFORCE_MCP_PAGE_LIMITS) {
  natural_height_in <- nrow(plot_mat) * ROW_HEIGHT_IN + header_margin_in
  if (natural_height_in > MAX_FIGURE_HEIGHT_IN) {
    effective_row_height_in <- (MAX_FIGURE_HEIGHT_IN - header_margin_in) / nrow(plot_mat)
    cat(sprintf(
      "Note: %d genes at %.3f in/row would need %.2f in of height (over MCP's %.0f in max). ",
      nrow(plot_mat), ROW_HEIGHT_IN, natural_height_in, MAX_FIGURE_HEIGHT_IN
    ))
    cat(sprintf("Row height compressed to %.4f in to fit.\n", effective_row_height_in))
    if (effective_row_height_in * 72 < 7) {  # ~7pt is about the smallest a row can be and still hold 7pt row labels
      cat("Warning: compressed row height is tight for 7pt row labels - consider raising",
          "ROW_FILTER_THRESHOLD to plot fewer genes.\n")
    }
  } else {
    effective_row_height_in <- ROW_HEIGHT_IN
  }
} else {
  effective_row_height_in <- ROW_HEIGHT_IN
}
FIGURE_HEIGHT_IN <- max(nrow(plot_mat) * effective_row_height_in + header_margin_in, 3)

# Width
heatmap_width_cm <- length(STAGE_ORDER) * CELL_WIDTH_CM + 2   
legend_width_cm  <- 4
natural_width_in <- (heatmap_width_cm + legend_width_cm) / 2.54
if (ENFORCE_MCP_PAGE_LIMITS) {
  fits <- STANDARD_COLUMN_WIDTHS_IN >= natural_width_in
  if (any(fits)) {
    FIGURE_WIDTH_IN <- min(STANDARD_COLUMN_WIDTHS_IN[fits])
  } else {
    FIGURE_WIDTH_IN <- max(STANDARD_COLUMN_WIDTHS_IN)
    cat("Note: natural width exceeds the 2-column maximum (7.2 in) - using 7.2 in; the figure will be tight.\n")
  }
} else {
  FIGURE_WIDTH_IN <- natural_width_in
}

#Legend styling - matches Heatmap_Protocol.r's legend_params block
legend_params <- list(
  title     = legend_title,
  title_gp  = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT_FAMILY),
  labels_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY)
)

ht <- Heatmap(
  plot_mat,
  name = "expr",
  col = col_fun,
  na_col = "grey95",
  rect_gp = gpar(col = "grey50", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  height = unit(nrow(plot_mat) * effective_row_height_in, "in"),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
  column_names_rot = 45,                                       # was 90 - matches Heatmap_Protocol.r
  column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
  heatmap_legend_param = legend_params,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(plot_mat[i, j])) {
      grid.text("-", x, y, gp = gpar(fontsize = 8, col = "grey60", fontfamily = FONT_FAMILY))
    }
  },
  width = unit(length(STAGE_ORDER) * CELL_WIDTH_CM, "cm")
)

pdf(OUTPUT_FILE, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN,
    family = FONT_FAMILY, useDingbats = FALSE)
draw(
  ht,
  column_title = PLOT_TITLE,
  column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY),
  column_title_side = "top",
  heatmap_legend_side = "right",
  padding = unit(c(2, 2, 2, 2), "mm")   # tight margins - MCP asks for no unnecessary white space
)
dev.off()

cat("\nSaved heatmap to", OUTPUT_FILE, "\n")

EXPORT_PNG_PREVIEW <- TRUE
if (EXPORT_PNG_PREVIEW) {
  png(sub("\\.pdf$", "_preview.png", OUTPUT_FILE),
      width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in",
      res = 1000, type = "cairo", family = FONT_FAMILY)
  draw(
    ht,
    column_title = PLOT_TITLE,
    column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY),
    column_title_side = "top",
    heatmap_legend_side = "right",
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  dev.off()
}
