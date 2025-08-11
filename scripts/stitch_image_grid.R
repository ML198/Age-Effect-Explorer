library(magick)
library(stringr)

label_from_filename <- function(fp) {
  nm <- tools::file_path_sans_ext(basename(fp))
  m  <- stringr::str_match(nm, "^DAS_(.*?)_GSEA(?:.*)?$")
  lab <- ifelse(is.na(m[,2]), nm, m[,2])
  trimws(gsub("_"," ", lab))
}

make_image_grids <- function(input_dir,
                             pattern = "^DAS_.*_GSEA\\.png$",
                             out_dir = file.path(input_dir, "grids"),
                             file_prefix = "GSEA_grid",
                             nrow = 2, ncol = 3,
                             cell_w = 1200, cell_h = 900,
                             bg = "white",
                             label = TRUE,
                             pointsize = 28,
                             add_index = TRUE,
                             index_size = 40,
                             index_offset = c(24, 14)) {   
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  files <- stringr::str_sort(files, numeric = TRUE)
  if (length(files) == 0) stop("No files matched pattern in input_dir.")
  
  k <- nrow * ncol
  pages <- split(seq_along(files), ceiling(seq_along(files) / k))
  
  for (pi in seq_along(pages)) {
    idxs  <- pages[[pi]]
    chunk <- files[idxs]
    
    imgs <- lapply(seq_along(chunk), function(j) {   # ← use j
      fp <- chunk[j]
      im <- image_read(fp)
      im <- image_resize(im, sprintf("%dx%d!", cell_w, cell_h))
      
      if (add_index) {
        num <- idxs[j]                               # ← now defined
        im <- image_annotate(
          im, as.character(num),
          gravity  = "northwest",
          location = sprintf("+%d+%d", index_offset[1], index_offset[2]),
          size     = index_size,
          color    = "black",
          boxcolor = "rgba(255,255,255,0.6)"
        )
      }
      
      im
    })
    
    rows     <- split(imgs, ceiling(seq_along(imgs) / ncol))
    row_imgs <- lapply(rows, function(r) image_append(image_join(r), stack = FALSE))
    grid     <- image_append(image_join(row_imgs), stack = TRUE)
    
    outfile <- file.path(out_dir, sprintf("%s_%02d.png", file_prefix, pi))
    image_write(grid, path = outfile, format = "png")
    message("Saved: ", outfile)
  }
}


app_dir <- here::here()
data_dir <- file.path(app_dir, "data")
tissue_file <- file.path(data_dir, "tissue_names.txt")
tissues <- readLines(tissue_file)

input_dir <- file.path(app_dir, "figures")

# 6 图一页（2x3）
make_image_grids(input_dir,
                 file_prefix = "GSEA_2x3",
                 nrow = 3, ncol = 2,
                 cell_w = 1200, cell_h = 900)

