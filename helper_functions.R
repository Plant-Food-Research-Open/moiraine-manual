## Have the intro code ready
read_targets_intro <- function(){
  targets_file <- file("_targets.R", open = "r")
  targets_lines <- readLines(targets_file)
  close(targets_file)

  intro_end <- stringr::str_which(targets_lines, "^\\)$")[[1]]

  c(targets_lines[seq_len(intro_end)], "")
}

## Read in list of chapters in order
get_chapters <- function(){
  quarto_yml <- file("_quarto.yml", open = "r")
  yml_lines <- readLines(quarto_yml)
  close(quarto_yml)

  ## Selecting section containing chapter list
  yaml_start <- stringr::str_which(yml_lines, "chapters")[2]
  yml_lines <- yml_lines[-(seq_len(yaml_start))]
  yaml_end <- stringr::str_which(yml_lines, "references.qmd")[1] - 1
  yml_lines <- yml_lines[seq_len(yaml_end)]

  yml_lines[stringr::str_which(yml_lines, "\\.qmd")] |>
    stringr::str_remove("\\s+- ")
}

## Read last chunk from each chapter
parse_chapter <- function(chapter){
  chapter_content <- parsermd::parse_rmd(chapter)

  chapter_title <- parsermd::as_tibble(chapter_content) |>
    dplyr::filter(!is.na(sec_h1)) |>
    dplyr::pull(sec_h1) |>
    unique() |>
    stringr::str_remove(" \\{.+")

  chapter_title_line <- paste0("# ", chapter_title, " ") |>
    stringr::str_pad(width = 80, side = "right", pad = "-")

  lines <- chapter_content |>
    parsermd::rmd_select(parsermd::has_label("*recap-targets-list")) |>
    parsermd::as_document()

  to_remove <- c(
    stringr::str_which(lines, "```"),
    stringr::str_which(lines, "^list\\("),
    stringr::str_which(lines, "^\\)$")
  )

  c(chapter_title_line, "", lines[-to_remove])
}

get_previous_content <- function(current_chapter) {
  current_chapter <- stringr::str_replace(current_chapter, "rmarkdown", "qmd")

  intro <- read_targets_intro()

  chapters <- get_chapters()
  curr_chap <- which(chapters == current_chapter)
  chapters <- chapters[seq_len(curr_chap - 1)]

  chapters_content <- purrr::map(chapters, parse_chapter)

  chapters_content <- unlist(chapters_content)

  if (chapters_content[length(chapters_content)] == "") {
    chapters_content <- chapters_content[-length(chapters_content)]
  }

  c(
    "::: {.cell}",
    "",
    "\`\`\`{.r .cell-code}",
    intro,
    "list(",
    chapters_content,
    ")",
    "\`\`\`",
    ":::"
  ) |>
    paste0(collapse = "\n")
}





## for each chapter, combine previous chapters' code
