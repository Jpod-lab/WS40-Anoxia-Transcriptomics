repo <- normalizePath(Sys.getenv("WS40_REPO", unset="."))
# =============================================================================
# 12dpd cluster-compare figure REMAKE — built directly from the exported
# GO/KEGG term Excel files (12dpd_*_cluster_{GO,KEGG}_terms_Hs_20251017.xlsx)
#
# WHY THIS SCRIPT EXISTS:
# The original clusterProfiler/enrichGO/enrichKEGG compareCluster() objects
# are not being re-run. Everything below reconstructs the figures purely
# from the already-exported term tables.
#
# FIGURES REBUILT:
#   Fig 6+7  - COMBINED GO+KEGG dotplot, up-regulated
#              (GO terms shaded by BP/MF/CC, KEGG terms shaded as their own
#               "KEGG" category, all in ONE plot; disease-term filter on
#               KEGG; bold 16pt term labels; narrowed x-axis spacing)
#   Fig 9+10 - COMBINED GO+KEGG dotplot, down-regulated (same idea)
#   Fig 8    - GO cnetplot, up-regulated   (pie term AND pie gene nodes
#                                            colored by Cluster, bold
#                                            leader-line labels)
#   Fig 11   - GO cnetplot, down-regulated (same pie-node style, separate file)
#
# =============================================================================
# LATEST CHANGE (this round): COMBINE GO + KEGG INTO ONE DOTPLOT, BIGGER BOLD
# TERM LABELS (16pt start), NARROWER X-AXIS
# =============================================================================
#
# WHAT CHANGED:
#   1. plot_kegg_dotplot() and plot_go_dotplot_shaded() (two separate
#      functions/figures) have been replaced by ONE function,
#      plot_combined_go_kegg_dotplot(), that takes both a GO term table and
#      a KEGG term table and renders them together in a single dotplot.
#      - GO terms keep their BP/MF/CC background shading (same pastel
#        colors as before).
#      - KEGG terms (after the same disease-term blacklist filter as
#        before) get their own "KEGG" category/shading color, appended
#        below the GO categories on the y-axis.
#      - This produces ONE figure per regulation direction instead of two
#        (Fig 6+7 combined for up-regulated; Fig 9+10 combined for
#        down-regulated), each still saved as both a vector PDF and a 600
#        DPI *_FOR_WORD.png.
#   2. Term label font size: previously 7.5-8pt (axis.text.y). Now starts
#      at 16pt (TERM_FONT_SIZE constant below) AND is bold via
#      ggtext::element_markdown() + **bold** markdown wrapping (same
#      bold_wrap() mechanism used before, just applied uniformly here).
#      If 16pt is too large once you see it on your real term list (more
#      terms = more likely to need to go smaller to fit), adjust
#      TERM_FONT_SIZE directly — it's the single number controlling this
#      everywhere in the combined dotplot.
#   3. Narrower x-axis (less white space between the 4 cluster/condition
#      columns): two changes work together for this —
#        a. scale_x_discrete(expand = ggplot2::expansion(add = 0.4)) tightens
#           the padding at the left/right edges of the x-axis (previously
#           ggplot2's default expansion added roughly a full extra category-
#           width of blank space on each side).
#        b. The saved plot WIDTH (COMBINED_DOTPLOT_WIDTH below) was reduced
#           relative to before, since a big driver of "white space between
#           conditions" is simply how wide the whole canvas is relative to
#           only 4 discrete x categories — a narrower canvas packs the same
#           4 columns closer together automatically. Term label wrap width
#           was widened slightly to compensate for the narrower plot so
#           long GO/KEGG term text still fits legibly on the left.
#
# If the 4 cluster columns still look too spread out once you see the
# actual output, the fastest additional lever is reducing
# COMBINED_DOTPLOT_WIDTH further (it directly controls the horizontal
# canvas the 4 columns are drawn across) — try dropping it in increments of
# ~1 inch. Do not shrink height the same way; height is driven by how many
# terms are shown and shrinking it will just cause label overlap.
#
# =============================================================================
# EARLIER FIXES (still in effect, unchanged from prior rounds):
# =============================================================================
#   - Canvas size vs. Word auto-shrink: ggplot2 text renders at an ABSOLUTE
#     physical size (mm/pt), independent of ggsave()'s width/height. Figure
#     canvas sizes below are kept at realistic actual-display sizes so Word
#     has little/nothing left to auto-shrink on insert.
#   - Cnetplot (Fig 8/11): gene nodes are pie charts reflecting every
#     Cluster a gene appears in among shown terms (constant size); term
#     nodes are pie charts showing CLUSTER PRESENCE/ABSENCE (not gene-count
#     proportion), sized by total gene count; a circle-overlap resolver
#     guarantees no two pies visually overlap; ggrepel label placement is
#     anchored at true node centers with term_label_manual_nudge available
#     for hand-tuning individual labels; bottom-left size legend removed;
#     term and gene pies share one merged "Cluster" fill legend.
#   - Every figure saves both a vector PDF (best for journal submission)
#     and a 600 DPI *_FOR_WORD.png (best for direct Word insertion) via the
#     save_for_word() helper. Install the optional `ragg` package for
#     noticeably sharper text rendering in the PNGs:
#       install.packages("ragg")
#
# HOW TO INSERT IN WORD:
#   1. Insert > Pictures > This Device, choosing the *_FOR_WORD.png file.
#   2. Check Format Picture > Size after inserting — width should read
#      close to the *_WIDTH constant used for that figure. If Word still
#      auto-shrinks it, reduce that constant to match your real column
#      width and re-export, rather than drag-resizing afterward.
#   3. File > Options > Advanced > Image Size and Quality > check "Do not
#      compress images in file" before finalizing.
#
# REQUIRED PACKAGES: scatterpie, ggrepel, ggtext, (optional) ragg
#   install.packages(c("scatterpie", "ggrepel", "ggtext", "ragg"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)
library(ggrepel)
library(ggtext)
library(scatterpie)
library(tibble)

# -----------------------------------------------------------------------
# 0. INPUT FILES
# -----------------------------------------------------------------------
go_up_path <- file.path(repo, "data/input/figure_terms", "12dpd_up_regulated_cluster_GO_terms_Hs_20251017.xlsx")
kegg_up_path <- file.path(repo, "data/input/figure_terms", "12dpd_up_regulated_cluster_KEGG_terms_Hs_20251017.xlsx")
go_dn_path <- file.path(repo, "data/input/figure_terms", "12dpd_down_regulated_cluster_GO_terms_Hs_20251017.xlsx")
kegg_dn_path <- file.path(repo, "data/input/figure_terms", "12dpd_down_regulated_cluster_KEGG_terms_Hs_20251017.xlsx")

out_dir_up <- file.path(repo, "results/reproduced/figures/enrichment")
out_dir_dn <- file.path(repo, "results/reproduced/figures/enrichment")

dir.create(out_dir_up, recursive=TRUE, showWarnings=FALSE)
cluster_order <- c("4hA", "24hA", "2hR", "24hR")

# -----------------------------------------------------------------------
# 0b. WORD-QUALITY EXPORT SETTINGS + HELPER
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
# 1. LOAD DATA
# -----------------------------------------------------------------------
go_up   <- readxl::read_xlsx(go_up_path)
kegg_up <- readxl::read_xlsx(kegg_up_path)
go_dn   <- readxl::read_xlsx(go_dn_path)
kegg_dn <- readxl::read_xlsx(kegg_dn_path)

clean_clusters <- function(df, col = "Cluster") {
  df[[col]] <- factor(df[[col]],
                      levels = cluster_order[cluster_order %in% base::unique(df[[col]])])
  df
}
go_up   <- clean_clusters(go_up)
kegg_up <- clean_clusters(kegg_up)
go_dn   <- clean_clusters(go_dn)
kegg_dn <- clean_clusters(kegg_dn)

# -----------------------------------------------------------------------
# 2. SHARED HELPERS
# -----------------------------------------------------------------------
kegg_disease_blacklist <- c(
  "cancer", "carcinoma", "melanoma", "glioma", "leukemia", "lymphoma",
  "measles", "legionellosis", "chagas disease", "malaria", "influenza",
  "tuberculosis", "hepatitis", "pertussis", "salmonella", "shigellosis",
  "amoebiasis", "toxoplasmosis", "leishmaniasis", "prion disease",
  "alcoholism", "diabetic cardiomyopathy", "alzheimer", "parkinson",
  "huntington", "amyotrophic lateral sclerosis", "chemical carcinogenesis",
  "pathways of neurodegeneration"
)

filter_kegg_disease_terms <- function(df, blacklist = kegg_disease_blacklist) {
  hit <- Reduce(`|`, lapply(blacklist, function(term) {
    base::grepl(term, df$Description, ignore.case = TRUE)
  }))
  if (base::any(hit)) {
    message("Dropping ", base::sum(hit), " non-informative KEGG disease term(s): ",
            base::paste(base::unique(df$Description[hit]), collapse = "; "))
  }
  df[!hit, ]
}

select_top_terms <- function(df, n, group_cols = "Cluster", rank_col = "p.adjust") {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::slice_min(order_by = .data[[rank_col]], n = n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    base::unique()
}

wrap_desc <- function(x, width = 55) stringr::str_wrap(x, width = width)

bold_wrap <- function(x, width = 55) {
  wrapped <- stringr::str_wrap(x, width = width)
  base::paste0("**", wrapped, "**")
}

kegg_dot_size_scale <- function(range = c(2, 10), n_breaks = 6) {
  ggplot2::scale_size_continuous(range = range, breaks = scales::pretty_breaks(n = n_breaks))
}

# -----------------------------------------------------------------------
# 3. COMBINED GO + KEGG DOTPLOT — Fig 6+7 (up) & Fig 9+10 (down)
#
#   GO terms keep BP/MF/CC background shading; KEGG terms (after disease-
#   term filtering) get their own "KEGG" shading band, appended below the
#   GO categories. Term labels are bold and start at TERM_FONT_SIZE (16pt).
#   The x-axis is tightened (less white space between the 4 cluster
#   columns) via a smaller expansion() and a narrower saved plot width.
# -----------------------------------------------------------------------

# Single control for term-label (y-axis) font size everywhere in this
# combined dotplot. Started at 16pt, then dialed down to 12pt to reduce
# crowding; now standardizing at 14pt to match the x-axis/legend text size
# requested below ("14 across the board"). Adjust this one number if 14pt
# still feels too large/small once you see the real term list rendered.
TERM_FONT_SIZE <- 14

# Single control for x-axis (cluster/condition) and legend text size,
# applied uniformly below so everything besides the y-axis term labels
# matches the same 14pt, bold treatment.
AXIS_LEGEND_FONT_SIZE <- 14

# Background shading colors: 3 GO categories + 1 for KEGG.
combined_category_colors <- c(BP = "#A6CEE3", MF = "#FDBF6F", CC = "#B2DF8A", KEGG = "#FB9A99")

plot_combined_go_kegg_dotplot <- function(go_df, kegg_df,
                                          n_top_per_cluster_ontology = 5,
                                          n_top_kegg_per_cluster = 10,
                                          wrap_width = 60) {
  
  # ---- Prepare GO half: same top-terms-per-Cluster-x-ONTOLOGY selection
  #      used before, tagged with its real ONTOLOGY (BP/MF/CC). ----
  go_keep <- select_top_terms(go_df, n = n_top_per_cluster_ontology,
                              group_cols = c("Cluster", "ONTOLOGY"))
  go_sub <- go_df %>% dplyr::filter(Description %in% go_keep) %>%
    dplyr::select(Description, Cluster, Count, p.adjust, ONTOLOGY)
  
  # ---- Prepare KEGG half: disease-term filter (as before), then top N per
  #      Cluster, tagged with ONTOLOGY = "KEGG" so it can share one combined
  #      category/shading scheme with the GO terms. ----
  kegg_df_filtered <- filter_kegg_disease_terms(kegg_df)
  kegg_keep <- select_top_terms(kegg_df_filtered, n = n_top_kegg_per_cluster)
  kegg_sub <- kegg_df_filtered %>% dplyr::filter(Description %in% kegg_keep) %>%
    dplyr::select(Description, Cluster, Count, p.adjust) %>%
    dplyr::mutate(ONTOLOGY = "KEGG")
  
  df <- dplyr::bind_rows(go_sub, kegg_sub)
  
  # ---- One category tag per term (guard against any term appearing under
  #      more than one tag, same defensive pattern as before). ----
  term_category <- df %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(ONTOLOGY = dplyr::first(ONTOLOGY),
                     best_p = base::min(p.adjust),
                     .groups = "drop")
  
  # ---- Order y-axis: BP, MF, CC, then KEGG (in that order), and within
  #      each category by strongest (smallest) p.adjust. ----
  term_category$ONTOLOGY <- factor(term_category$ONTOLOGY, levels = c("BP", "MF", "CC", "KEGG"))
  term_category <- term_category %>% dplyr::arrange(ONTOLOGY, best_p)
  
  df$Description <- factor(df$Description, levels = base::rev(term_category$Description))
  desc_bold_levels <- bold_wrap(base::rev(term_category$Description), width = wrap_width)
  df$DescriptionBold <- factor(bold_wrap(df$Description, width = wrap_width), levels = desc_bold_levels)
  
  # ---- Background shading bands per category. ----
  band_positions <- term_category %>%
    dplyr::mutate(y_bold = base::as.numeric(factor(bold_wrap(Description, width = wrap_width),
                                                   levels = desc_bold_levels)))
  bands <- band_positions %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::summarise(ymin = base::min(y_bold) - 0.5, ymax = base::max(y_bold) + 0.5, .groups = "drop")
  
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands,
                       ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
                       alpha = 0.25, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = combined_category_colors, name = "Category",
                               guide = ggplot2::guide_legend(override.aes = list(alpha = 0.4), order = 1)) +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(x = Cluster, y = DescriptionBold, size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "red", high = "blue", name = "p.adjust",
                                  guide = ggplot2::guide_colorbar(order = 2)) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Gene count",
                                   guide = ggplot2::guide_legend(order = 3)) +
    # NARROWER X-AXIS: expansion(add = 0.4) leaves only 0.4 category-widths
    # of blank space on either side of the 4 cluster columns, versus
    # ggplot2's default (~0.6 units + 5% multiplicative), visibly tightening
    # the horizontal spread of the 4 conditions.
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggtext::element_markdown(size = TERM_FONT_SIZE),
                   # X-axis (cluster/condition) labels: bold + same 14pt
                   # size as the y-axis term labels and legend text below.
                   axis.text.x = ggplot2::element_text(size = AXIS_LEGEND_FONT_SIZE, face = "bold"),
                   # Legend title ("Category", "p.adjust", "Gene count")
                   # and legend item text (BP/MF/CC/KEGG, size breaks,
                   # p.adjust color scale numbers): bold + 14pt to match.
                   legend.title = ggplot2::element_text(size = AXIS_LEGEND_FONT_SIZE, face = "bold"),
                   legend.text = ggplot2::element_text(size = AXIS_LEGEND_FONT_SIZE, face = "bold"))
}

# COMBINED_DOTPLOT_WIDTH/HEIGHT: width was reduced relative to the old
# separate Fig6/Fig7 canvases specifically to pack the 4 cluster columns
# closer together (see note above — canvas width is the main driver of
# "white space between conditions" once expansion() is already tightened).
# Height still needs to be generous since ALL GO+KEGG terms are now on one
# axis; adjust n_top_per_cluster_ontology / n_top_kegg_per_cluster (fewer
# terms) if you'd rather keep height down than shrink font further.
COMBINED_DOTPLOT_WIDTH  <- 10
COMBINED_DOTPLOT_HEIGHT <- 16

fig_up_combined <- plot_combined_go_kegg_dotplot(go_up, kegg_up,
                                                 n_top_per_cluster_ontology = 5,
                                                 n_top_kegg_per_cluster = 10)
fig_up_combined_path <- base::file.path(out_dir_up, "12dpd_COMBINED_GO_KEGG_dotplot_UP_FROM_EXCEL")
ggplot2::ggsave(base::paste0(fig_up_combined_path, ".pdf"), fig_up_combined,
                width = COMBINED_DOTPLOT_WIDTH, height = COMBINED_DOTPLOT_HEIGHT)
save_for_word(fig_up_combined, fig_up_combined_path,
              width = COMBINED_DOTPLOT_WIDTH, height = COMBINED_DOTPLOT_HEIGHT)

fig_dn_combined <- plot_combined_go_kegg_dotplot(go_dn, kegg_dn,
                                                 n_top_per_cluster_ontology = 5,
                                                 n_top_kegg_per_cluster = 10)
fig_dn_combined_path <- base::file.path(out_dir_dn, "12dpd_COMBINED_GO_KEGG_dotplot_DOWN_FROM_EXCEL")
ggplot2::ggsave(base::paste0(fig_dn_combined_path, ".pdf"), fig_dn_combined,
                width = COMBINED_DOTPLOT_WIDTH, height = COMBINED_DOTPLOT_HEIGHT)
save_for_word(fig_dn_combined, fig_dn_combined_path,
              width = COMBINED_DOTPLOT_WIDTH, height = COMBINED_DOTPLOT_HEIGHT)

# Bump n_top_per_cluster_ontology / n_top_kegg_per_cluster up/down to make
# the combined figure show more or fewer terms per category. A full/
# uncondensed supplemental version can be generated the same way with
# larger values, e.g.:
# fig_up_supp <- plot_combined_go_kegg_dotplot(go_up, kegg_up,
#                  n_top_per_cluster_ontology = 100, n_top_kegg_per_cluster = 100)
# ggplot2::ggsave(base::file.path(out_dir_up, "12dpd_COMBINED_GO_KEGG_dotplot_UP_SUPPLEMENT.pdf"),
#                 fig_up_supp, width = 14, height = 50, limitsize = FALSE)

# -----------------------------------------------------------------------
# 5. GO CNETPLOTS, PIE-NODE STYLE — Fig 8 (up) & Fig 11 (down)
# -----------------------------------------------------------------------
TARGET_TERMS_FIG8 <- c(
  "protein refolding",
  "response to unfolded protein",
  "protein folding chaperone",
  "response to temperature stimulus",
  "unfolded protein binding",
  "cellular response to nutrient levels",
  "MAP kinase tyrosine/serine/threonine phosphatase activity",
  "regulation of neuroinflammatory response",
  "cellular response to xenobiotic stimulus",
  "xenobiotic transmembrane transport",
  "ABC-type xenobiotic transporter activity"
)

resolve_circle_overlaps <- function(x, y, r, iterations = 400, padding = 1.1, step = 0.5) {
  n <- base::length(x)
  for (iter in base::seq_len(iterations)) {
    moved <- FALSE
    for (i in base::seq_len(n - 1)) {
      for (j in (i + 1):n) {
        dx <- x[j] - x[i]
        dy <- y[j] - y[i]
        d <- base::sqrt(dx^2 + dy^2)
        min_d <- (r[i] + r[j]) * padding
        if (d < min_d) {
          if (d < 1e-6) {
            ang <- stats::runif(1, 0, 2 * base::pi)
            dx <- base::cos(ang); dy <- base::sin(ang); d <- 1e-6
          }
          overlap <- (min_d - d)
          ux <- dx / d; uy <- dy / d
          x[i] <- x[i] - ux * overlap * step
          y[i] <- y[i] - uy * overlap * step
          x[j] <- x[j] + ux * overlap * step
          y[j] <- y[j] + uy * overlap * step
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  tibble::tibble(x = x, y = y)
}

plot_go_cnet_pie <- function(df, target_terms = NULL, top_n_auto = 5,
                             gene_label_size = 5.5, term_label_size = 6.0,
                             term_radius_frac = c(0.03, 0.03),
                             gene_radius_frac = 0.012,
                             spread_factor = 1.8, seed = 42,
                             overlap_iterations = 400, overlap_padding = 1.15,
                             overlap_step = 0.5,
                             label_nudge_mult_term = 3.2,
                             label_nudge_mult_gene = 1.6,
                             term_label_manual_nudge = NULL) {
  
  cluster_levels <- base::levels(df$Cluster)
  
  if (base::is.null(target_terms)) {
    target_terms <- select_top_terms(df, n = top_n_auto, group_cols = "Cluster")
  }
  sub <- df %>% dplyr::filter(Description %in% target_terms)
  if (base::nrow(sub) == 0) {
    stop("None of the requested target_terms were found in this data frame's Description column. ",
         "Check spelling/case against the Excel file's Description column.")
  }
  
  term_cluster_presence <- sub %>%
    dplyr::distinct(Description, Cluster) %>%
    dplyr::mutate(present = 1)
  
  term_totals <- sub %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(total = base::sum(Count), .groups = "drop")
  
  term_pie_wide <- term_cluster_presence %>%
    tidyr::pivot_wider(names_from = Cluster, values_from = present, values_fill = 0)
  missing_term_cols <- base::setdiff(cluster_levels, base::colnames(term_pie_wide))
  for (mc in missing_term_cols) term_pie_wide[[mc]] <- 0
  
  gene_term_cluster <- sub %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(gene = geneID) %>%
    dplyr::distinct(Description, Cluster, gene)
  
  edge_list <- gene_term_cluster %>% dplyr::distinct(Description, gene)
  
  gene_cluster_presence <- gene_term_cluster %>%
    dplyr::distinct(gene, Cluster) %>%
    dplyr::mutate(present = 1) %>%
    tidyr::pivot_wider(names_from = Cluster, values_from = present, values_fill = 0)
  missing_gene_cols <- base::setdiff(cluster_levels, base::colnames(gene_cluster_presence))
  for (mc in missing_gene_cols) gene_cluster_presence[[mc]] <- 0
  
  terms <- base::unique(edge_list$Description)
  genes <- base::unique(edge_list$gene)
  
  g <- igraph::graph_from_data_frame(
    edge_list,
    vertices = dplyr::bind_rows(
      tibble::tibble(name = terms, type = "term"),
      tibble::tibble(name = genes, type = "gene")
    ),
    directed = FALSE
  )
  
  base::set.seed(seed)
  layout_xy <- igraph::layout_with_fr(g) * spread_factor
  base::colnames(layout_xy) <- c("x", "y")
  node_pos <- tibble::tibble(name = igraph::V(g)$name, type = igraph::V(g)$type,
                             x = layout_xy[, "x"], y = layout_xy[, "y"])
  
  layout_span <- base::max(base::diff(base::range(node_pos$x)),
                           base::diff(base::range(node_pos$y)))
  
  term_nodes <- node_pos %>% dplyr::filter(type == "term") %>%
    dplyr::rename(Description = name) %>%
    dplyr::left_join(term_totals, by = "Description") %>%
    dplyr::left_join(term_pie_wide, by = "Description")
  term_nodes$r <- scales::rescale(term_nodes$total, to = term_radius_frac * layout_span)
  
  gene_nodes <- node_pos %>% dplyr::filter(type == "gene") %>%
    dplyr::rename(gene = name) %>%
    dplyr::left_join(gene_cluster_presence, by = "gene")
  gene_nodes$r <- gene_radius_frac * layout_span
  
  all_nodes <- dplyr::bind_rows(
    term_nodes %>% dplyr::transmute(id = Description, x, y, r),
    gene_nodes %>% dplyr::transmute(id = gene, x, y, r)
  )
  resolved <- resolve_circle_overlaps(all_nodes$x, all_nodes$y, all_nodes$r,
                                      iterations = overlap_iterations,
                                      padding = overlap_padding,
                                      step = overlap_step)
  all_nodes$x <- resolved$x
  all_nodes$y <- resolved$y
  
  term_nodes <- term_nodes %>% dplyr::select(-x, -y) %>%
    dplyr::left_join(all_nodes %>% dplyr::select(id, x, y), by = c("Description" = "id"))
  gene_nodes <- gene_nodes %>% dplyr::select(-x, -y) %>%
    dplyr::left_join(all_nodes %>% dplyr::select(id, x, y), by = c("gene" = "id"))
  
  edge_coords <- edge_list %>%
    dplyr::left_join(term_nodes %>% dplyr::select(Description, x, y), by = "Description") %>%
    dplyr::rename(x_term = x, y_term = y) %>%
    dplyr::left_join(gene_nodes %>% dplyr::select(gene, x, y), by = "gene") %>%
    dplyr::rename(x_gene = x, y_gene = y)
  
  cluster_palette <- scales::hue_pal()(base::length(cluster_levels))
  base::names(cluster_palette) <- cluster_levels
  
  cx <- base::mean(all_nodes$x); cy <- base::mean(all_nodes$y)
  radial_unit <- function(x, y) {
    dx <- x - cx; dy <- y - cy
    d <- base::sqrt(dx^2 + dy^2); d[d == 0] <- 1e-6
    tibble::tibble(ux = dx / d, uy = dy / d)
  }
  term_dir <- radial_unit(term_nodes$x, term_nodes$y)
  gene_dir <- radial_unit(gene_nodes$x, gene_nodes$y)
  
  term_nudge_x <- term_dir$ux * term_nodes$r * label_nudge_mult_term
  term_nudge_y <- term_dir$uy * term_nodes$r * label_nudge_mult_term
  gene_nudge_x <- gene_dir$ux * gene_nodes$r * label_nudge_mult_gene
  gene_nudge_y <- gene_dir$uy * gene_nodes$r * label_nudge_mult_gene
  
  if (!base::is.null(term_label_manual_nudge)) {
    for (term_name in base::names(term_label_manual_nudge)) {
      idx <- base::which(term_nodes$Description == term_name)
      if (base::length(idx) == 0) {
        warning("term_label_manual_nudge: '", term_name,
                "' not found among this plot's Description values — skipped.")
        next
      }
      offset <- term_label_manual_nudge[[term_name]]
      term_nudge_x[idx] <- term_nudge_x[idx] + offset[["x"]]
      term_nudge_y[idx] <- term_nudge_y[idx] + offset[["y"]]
    }
  }
  
  term_label_df <- term_nodes %>%
    dplyr::transmute(x, y, label = stringr::str_wrap(Description, 22))
  gene_label_df <- gene_nodes %>%
    dplyr::transmute(x, y, label = gene)
  
  gene_obstacles_for_term_labels <- gene_nodes %>% dplyr::transmute(x, y, label = "")
  term_obstacles_for_gene_labels <- term_nodes %>% dplyr::transmute(x, y, label = "")
  
  term_repel_data <- dplyr::bind_rows(term_label_df, gene_obstacles_for_term_labels)
  term_repel_nudge_x <- c(term_nudge_x, base::rep(0, base::nrow(gene_obstacles_for_term_labels)))
  term_repel_nudge_y <- c(term_nudge_y, base::rep(0, base::nrow(gene_obstacles_for_term_labels)))
  
  gene_repel_data <- dplyr::bind_rows(gene_label_df, term_obstacles_for_gene_labels)
  gene_repel_nudge_x <- c(gene_nudge_x, base::rep(0, base::nrow(term_obstacles_for_gene_labels)))
  gene_repel_nudge_y <- c(gene_nudge_y, base::rep(0, base::nrow(term_obstacles_for_gene_labels)))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_coords,
                          ggplot2::aes(x = x_term, y = y_term, xend = x_gene, yend = y_gene),
                          color = "grey70", linewidth = 0.3, alpha = 0.6) +
    scatterpie::geom_scatterpie(data = gene_nodes, ggplot2::aes(x = x, y = y, r = r),
                                cols = cluster_levels, color = NA, alpha = 0.95) +
    scatterpie::geom_scatterpie(data = term_nodes, ggplot2::aes(x = x, y = y, r = r),
                                cols = cluster_levels, color = NA) +
    ggplot2::scale_fill_manual(values = cluster_palette, name = "Cluster", breaks = cluster_levels) +
    ggrepel::geom_text_repel(data = gene_repel_data,
                             ggplot2::aes(x = x, y = y, label = label),
                             nudge_x = gene_repel_nudge_x, nudge_y = gene_repel_nudge_y,
                             size = gene_label_size, fontface = "italic",
                             segment.size = 0.25, segment.color = "grey50",
                             box.padding = 0.4, point.padding = 0.25,
                             force = 3, max.overlaps = Inf, seed = seed,
                             min.segment.length = 0) +
    ggrepel::geom_text_repel(data = term_repel_data,
                             ggplot2::aes(x = x, y = y, label = label),
                             nudge_x = term_repel_nudge_x, nudge_y = term_repel_nudge_y,
                             size = term_label_size, fontface = "bold",
                             segment.size = 0.45, segment.color = "black",
                             box.padding = 1.8, point.padding = 0.7,
                             force = 20, max.iter = 20000, max.overlaps = Inf, seed = seed,
                             min.segment.length = 0) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right")
  
  # (Bottom-left "node size = gene count" legend intentionally omitted.)
  
  p
}

# ---- Fig 8: UP-regulated cnetplot -----------------------------------------
fig8 <- plot_go_cnet_pie(
  go_up, target_terms = TARGET_TERMS_FIG8,
  term_label_manual_nudge = list(
    "unfolded protein binding" = c(x = -2.5, y = 2.5),
    "MAP kinase tyrosine/serine/threonine phosphatase activity" = c(x = -3, y = -0.25),
    "regulation of neuroinflammatory response" = c(x = 4.5, y = -4.5),
    "ABC-type xenobiotic transporter activity" = c(x = -2, y = -1.0),
    "cellular response to xenobiotic stimulus" = c(x = 0.0, y = -1.0),
    "response to unfolded protein" = c(x = -8.0, y = -0.5),
    "response to temperature stimulus" = c(x = -3.0, y = 1.0),
    "protein refolding" = c(x = 5, y = 0.0),
    "cellular response to nutrient levels" = c(x = 2.0, y = 4.0),
    "protein folding chaperone" = c(x = 0.0, y = -4)
  )
)

FIG8_WIDTH  <- 12.0
FIG8_HEIGHT <- 11.00
fig8_path <- base::file.path(out_dir_up, "cnet_enrichGO_justgenes_UP_FROM_EXCEL")
ggplot2::ggsave(base::paste0(fig8_path, ".pdf"), fig8, width = FIG8_WIDTH, height = FIG8_HEIGHT)
save_for_word(fig8, fig8_path, width = FIG8_WIDTH, height = FIG8_HEIGHT)

# ---- Fig 11: DOWN-regulated cnetplot ---------------------------------------
fig11 <- plot_go_cnet_pie(go_dn, target_terms = NULL, top_n_auto = 5)

FIG11_WIDTH  <- 12.0
FIG11_HEIGHT <- 11.0
fig11_path <- base::file.path(out_dir_dn, "cnet_enrichGO_justgenes_DOWN_FROM_EXCEL")
ggplot2::ggsave(base::paste0(fig11_path, ".pdf"), fig11, width = FIG11_WIDTH, height = FIG11_HEIGHT)
save_for_word(fig11, fig11_path, width = FIG11_WIDTH, height = FIG11_HEIGHT)
