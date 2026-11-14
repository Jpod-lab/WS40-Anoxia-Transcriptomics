repo <- normalizePath(Sys.getenv("WS40_REPO", unset="."))
# =============================================================================
# 12dpd 4-PANEL VOLCANO PLOT REMAKE (Panels A-D)
# Built directly from the DESeq2 results CSVs:
#   Panel A = 12dpd_SA_t0_res_sig_fullset.csv         (4hA vs t0)
#   Panel B = 12dpd_LA_t0_res_sig_fullset.csv         (24hA vs t0)
#   Panel C = 12dpd_SR_t0_res_sig_fullset.csv         (2hR vs t0)
#   Panel D = 12dpd_LR_t0_res_sig_fullset.csv         (24hR vs t0)
# Gene symbol labels are pulled from Orthofinder_blast_merge_Hs_230628_MapIDS.csv
# (joined on Locus_Al == the DESeq2 row name; falls back to the row name
# itself — i.e. the Locus/LOC id — when no human ortholog symbol exists).
#
# =============================================================================
# LATEST FIX (this round): LABELS OVERLAPPING THE BLACK NS POINT CLOUD
# =============================================================================
# WHAT HAPPENED: labels no longer overlap EACH OTHER (previous round's
# fix), but several still land directly on top of the dense black
# non-significant (NS) point cloud, making them hard to read against all
# the black dots (most visible in Panel D around MMP9/ABCB11).
#
# WHY: ggrepel only knows to avoid overlap between the items you actually
# hand it. Since only the ~8 labeled genes per panel were passed to
# geom_text_repel(), it had no awareness that a huge cloud of OTHER points
# (the NS genes, plotted via a separate geom_point() layer) existed at
# all — so it happily placed labels right on top of that cloud whenever
# doing so kept the labels away from each other.
#
# THE FIX: build_obstacle_df() constructs a set of invisible "obstacle"
# points — every other point in the panel (NS and unlabeled significant
# genes alike) that falls within the visible plot region — with an EMPTY
# label (""). These are combined with the real label_df and passed to
# ggrepel TOGETHER. Because ggrepel treats every row it receives as a
# label-placement candidate, the blank-label obstacle rows still occupy
# space that ggrepel must route real labels around, even though nothing
# visibly prints for them. This is the standard trick for forcing ggrepel
# to avoid arbitrary background points, not just other visible labels.
#
# PERFORMANCE NOTE: panels B/C/D have on the order of ~20,000-26,000 total
# points, which is too many obstacles to hand ggrepel directly without a
# major slowdown (collision-avoidance cost scales with how many points it
# has to check against). MAX_OBSTACLES caps how many obstacle points are
# used per panel (default 3000, random-sampled with a fixed seed for
# reproducibility) — enough to make the dense NS region something ggrepel
# actively routes around, without ballooning runtime. If labels are still
# landing on the cloud after re-running, raise MAX_OBSTACLES (e.g. to
# 6000-8000); if the script runs too slowly, lower it.
#
# =============================================================================
# EARLIER FIXES (still in effect):
# =============================================================================
#   - Each panel gets its OWN x-axis range from its own data (+ padding),
#     so Panel C's genuine outlier (iws1, log2FC ~ -19.8) doesn't force
#     every panel onto an unnecessarily wide, mostly-empty axis.
#   - Stronger ggrepel settings (force, max.iter, box.padding,
#     point.padding, fixed seed, min.segment.length = 0) resolve
#     label-vs-label overlap; Panel D gets extra panel-specific tuning
#     (higher force, lower force_pull, more y-axis headroom) for its
#     especially tight low-value cluster.
#   - Bold A/B/C/D tags via patchwork; every figure saves both a vector
#     PDF and a 600 DPI *_FOR_WORD.png via save_for_word().
#
# REQUIRED PACKAGES: readr, dplyr, stringr, ggplot2, ggrepel, patchwork,
#                     scales, (optional) ragg
#   install.packages(c("readr","dplyr","stringr","ggplot2","ggrepel",
#                       "patchwork","scales","ragg"))
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

# -----------------------------------------------------------------------
# 0. INPUT FILES — update these paths to match your local system.
# -----------------------------------------------------------------------
sa_path <- file.path(repo, "data/input/historical_threshold", "12dpd_SA_t0_res_sig_fullset.csv")
la_path <- file.path(repo, "data/input/historical_threshold", "12dpd_LA_t0_res_sig_fullset.csv")
sr_path <- file.path(repo, "data/input/historical_threshold", "12dpd_SR_t0_res_sig_fullset.csv")
lr_path <- file.path(repo, "data/input/historical_threshold", "12dpd_LR_t0_res_sig_fullset.csv")
ortho_path <- file.path(repo, "data/input/Orthofinder_blast_merge_Hs_230628_MapIDS.csv")

out_dir <- file.path(repo, "results/reproduced/figures/volcano")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------
# 0b. WORD-QUALITY EXPORT SETTINGS + HELPER (same pattern as the other
#     figure scripts in this project, for consistency).
# -----------------------------------------------------------------------
WORD_PNG_DPI <- 600
.have_ragg <- base::requireNamespace("ragg", quietly = TRUE)

save_for_word <- function(plot, path_no_ext, width, height, dpi = WORD_PNG_DPI) {
  png_path <- base::paste0(path_no_ext, "_FOR_WORD.png")
  if (.have_ragg) {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = dpi, device = ragg::agg_png, units = "in", bg = "white")
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = dpi, device = "png", units = "in", bg = "white")
  }
  message("Word-ready PNG saved (", dpi, " DPI, ", width, "x", height, "in): ", png_path)
}

# -----------------------------------------------------------------------
# 1. SIGNIFICANCE THRESHOLDS (adjust here if your criteria differ)
# -----------------------------------------------------------------------
PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1     # i.e. 2-fold change; the two near-zero vertical
# reference lines in the plot mark +/- this value

# How many top (smallest padj) significant genes to label per panel.
N_LABELS_PER_PANEL <- 8

# Max number of "invisible obstacle" background points fed to ggrepel per
# panel (see LATEST FIX note above). Higher = more thorough avoidance of
# the dense point cloud, but slower. Lower = faster, but labels may more
# easily still land on an unrepresented part of the cloud.
MAX_OBSTACLES <- 3000

# -----------------------------------------------------------------------
# 2. LOAD DESeq2 RESULT TABLES
# -----------------------------------------------------------------------
read_deseq_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  names(df)[1] <- "GeneID"
  df
}

sa <- read_deseq_csv(sa_path)
la <- read_deseq_csv(la_path)
sr <- read_deseq_csv(sr_path)
lr <- read_deseq_csv(lr_path)

# -----------------------------------------------------------------------
# 3. LOAD ORTHOFINDER MAP AND BUILD A CLEAN Locus_Al -> Gene_Hs LOOKUP
# -----------------------------------------------------------------------
ortho <- readr::read_csv(ortho_path, show_col_types = FALSE)

gene_symbol_lookup <- ortho %>%
  dplyr::filter(!base::is.na(Locus_Al), Locus_Al != "") %>%
  dplyr::filter(!base::is.na(Gene_Hs), Gene_Hs != "") %>%
  dplyr::distinct(Locus_Al, .keep_all = TRUE) %>%
  dplyr::select(Locus_Al, Gene_Hs)

add_gene_symbol <- function(df) {
  df %>%
    dplyr::left_join(gene_symbol_lookup, by = c("GeneID" = "Locus_Al")) %>%
    dplyr::mutate(display_label = dplyr::if_else(base::is.na(Gene_Hs) | Gene_Hs == "",
                                                 GeneID, Gene_Hs))
}

sa <- add_gene_symbol(sa)
la <- add_gene_symbol(la)
sr <- add_gene_symbol(sr)
lr <- add_gene_symbol(lr)

# -----------------------------------------------------------------------
# 4. CLASSIFY SIGNIFICANCE (Up / Down / NS) AND DROP NA-padj ROWS
# -----------------------------------------------------------------------
classify_volcano <- function(df, padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LOG2FC_CUTOFF) {
  df %>%
    dplyr::filter(!base::is.na(padj)) %>%
    dplyr::mutate(
      neg_log10_padj = -base::log10(padj),
      category = dplyr::case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      ),
      category = factor(category, levels = c("Up", "Down", "NS"))
    )
}

sa_v <- classify_volcano(sa)
la_v <- classify_volcano(la)
sr_v <- classify_volcano(sr)
lr_v <- classify_volcano(lr)

volcano_colors <- c(Up = "red", Down = "blue", NS = "black")

# -----------------------------------------------------------------------
# 4b. OBSTACLE BUILDER (this round's fix)
# -----------------------------------------------------------------------
# Builds a set of invisible (blank-label) points representing the
# surrounding point cloud, restricted to the panel's actual visible
# x/y range, and downsampled to MAX_OBSTACLES if needed. Passing these
# alongside the real labels forces ggrepel to steer labels around the
# dense point cloud instead of placing them directly on top of it.
build_obstacle_df <- function(df, label_df, x_min, x_max, y_max,
                              max_obstacles = MAX_OBSTACLES, seed = 1) {
  obstacles <- df %>%
    dplyr::filter(!(GeneID %in% label_df$GeneID)) %>%
    dplyr::filter(log2FoldChange >= x_min, log2FoldChange <= x_max,
                  neg_log10_padj <= y_max) %>%
    dplyr::transmute(log2FoldChange, neg_log10_padj, category,
                     display_label = "")
  
  if (base::nrow(obstacles) > max_obstacles) {
    base::set.seed(seed)
    obstacles <- obstacles %>% dplyr::slice_sample(n = max_obstacles)
  }
  obstacles
}

# -----------------------------------------------------------------------
# 5. VOLCANO PLOT BUILDER — one call per panel, EACH gets its own x-axis
#    range computed from its own data (with padding), and now also
#    obstacle-aware label placement (see 4b above).
# -----------------------------------------------------------------------
plot_volcano_panel <- function(df, n_labels = N_LABELS_PER_PANEL,
                               point_size = 1.1, label_size = 3.2,
                               x_pad_frac = 0.08,
                               repel_seed = 42,
                               force = 12, force_pull = 1,
                               box_padding = 0.6, point_padding = 0.3,
                               y_pad_top_frac = 0.08,
                               max_obstacles = MAX_OBSTACLES) {
  
  x_min_raw <- base::min(df$log2FoldChange)
  x_max_raw <- base::max(df$log2FoldChange)
  x_range   <- x_max_raw - x_min_raw
  x_min <- base::floor(x_min_raw - x_range * x_pad_frac)
  x_max <- base::ceiling(x_max_raw + x_range * x_pad_frac)
  
  y_max_raw <- base::max(df$neg_log10_padj, na.rm = TRUE)
  y_max <- y_max_raw * (1 + y_pad_top_frac)
  
  label_df <- df %>%
    dplyr::filter(category != "NS") %>%
    dplyr::arrange(padj) %>%
    utils::head(n_labels)
  
  # Invisible background obstacles (this round's fix) — combined with the
  # real labels so ggrepel actively avoids the dense point cloud, not just
  # other visible labels.
  obstacle_df <- build_obstacle_df(df, label_df, x_min, x_max, y_max,
                                   max_obstacles = max_obstacles)
  repel_data <- dplyr::bind_rows(
    label_df %>% dplyr::select(log2FoldChange, neg_log10_padj, category, display_label),
    obstacle_df
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = category), size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = volcano_colors, guide = "none") +
    ggplot2::geom_hline(yintercept = -base::log10(PADJ_CUTOFF), color = "red", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), color = "black", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = repel_data,
                             ggplot2::aes(label = display_label, color = category),
                             size = label_size, max.overlaps = Inf,
                             force = force, force_pull = force_pull, max.iter = 20000,
                             box.padding = box_padding, point.padding = point_padding,
                             segment.size = 0.3, min.segment.length = 0,
                             seed = repel_seed,
                             show.legend = FALSE) +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max)) +
    ggplot2::scale_y_continuous(limits = c(NA, y_max)) +
    ggplot2::labs(x = "log2FoldChange", y = "-log10(padj)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

panel_A <- plot_volcano_panel(sa_v)
panel_B <- plot_volcano_panel(la_v)
panel_C <- plot_volcano_panel(sr_v)

# Panel D-SPECIFIC OVERRIDE: this panel's tight low-value cluster
# (MMP9/ABCB11/MIOX all crowded near the bottom of a much smaller y-axis
# range than the other 3 panels) needed extra help beyond the global
# defaults. force raised further (12->18), force_pull lowered (1->0.5) so
# labels can drift further from crowded neighbors, box_padding increased
# (0.6->0.9), and extra y-axis headroom (y_pad_top_frac->0.15) gives
# ggrepel more free vertical space to spread this cluster into. Obstacle
# avoidance (new this round) applies here too via the shared function.
panel_D <- plot_volcano_panel(lr_v, force = 18, force_pull = 0.5,
                              box_padding = 0.9, y_pad_top_frac = 0.15)

# -----------------------------------------------------------------------
# 6. COMBINE INTO A 2x2 GRID WITH BOLD A/B/C/D TAGS (patchwork)
# -----------------------------------------------------------------------
combined_volcano <- (panel_A | panel_B) / (panel_C | panel_D) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))

VOLCANO_WIDTH  <- 9
VOLCANO_HEIGHT <- 8
volcano_path <- base::file.path(out_dir, "12dpd_volcano_4panel_ABCD")
ggplot2::ggsave(base::paste0(volcano_path, ".pdf"), combined_volcano,
                width = VOLCANO_WIDTH, height = VOLCANO_HEIGHT)
save_for_word(combined_volcano, volcano_path,
              width = VOLCANO_WIDTH, height = VOLCANO_HEIGHT)
# -----------------------------------------------------------------------
# 7. EXPORT THE FOUR VOLCANO PANELS INDIVIDUALLY
# -----------------------------------------------------------------------
# These exports intentionally omit panel letters. The final combined
# A-H figure can add one consistent set of panel tags during assembly.

INDIVIDUAL_VOLCANO_WIDTH  <- 5
INDIVIDUAL_VOLCANO_HEIGHT <- 4.5

individual_volcano_plots <- list(
  "12dpd_volcano_4hA_vs_t0"  = panel_A,
  "12dpd_volcano_24hA_vs_t0" = panel_B,
  "12dpd_volcano_2hR_vs_t0"  = panel_C,
  "12dpd_volcano_24hR_vs_t0" = panel_D
)

for (plot_name in names(individual_volcano_plots)) {
  
  current_plot <- individual_volcano_plots[[plot_name]]
  
  output_path <- base::file.path(
    out_dir,
    plot_name
  )
  
  # Vector PDF
  ggplot2::ggsave(
    filename = base::paste0(output_path, ".pdf"),
    plot = current_plot,
    width = INDIVIDUAL_VOLCANO_WIDTH,
    height = INDIVIDUAL_VOLCANO_HEIGHT,
    units = "in",
    bg = "white"
  )
  
  # 600-DPI PNG using your existing save_for_word() helper
  save_for_word(
    plot = current_plot,
    path_no_ext = output_path,
    width = INDIVIDUAL_VOLCANO_WIDTH,
    height = INDIVIDUAL_VOLCANO_HEIGHT,
    dpi = WORD_PNG_DPI
  )
}

message(
  "Finished exporting all four individual volcano plots to: ",
  out_dir
)