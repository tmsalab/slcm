custom_user_spelling_dictionary = function(custom_words, pkg) {
  
  # Check if package has aspell dictionary
  aspell_location = ".aspell"
  if(!dir.exists(aspell_location)) dir.create(aspell_location)
  
  # For newer versions of R, we create an R-level dictionary using an RDS file.
  # This is read using
  saveRDS(custom_words, file = file.path(aspell_location, paste0(pkg, ".rds")), version = 2)
  
  defaults_file_text = paste0('Rd_files <- vignettes <- R_files <- description <-
  list(encoding = "UTF-8",
       language = "en",
       dictionaries = c("en_stats", "',pkg,'"))
  ')
  writeLines(defaults_file_text, file.path(aspell_location, "defaults.R"))
  
  # For the {spelling} package, write out a simple WORDLIST.
  writeLines(custom_words, file.path("inst", "WORDLIST"))
  
  # Provide backward compatibility support for base R
  # Note, must be run from package root.
  # For more details, see: ?"aspell-utils"
  utils::aspell_write_personal_dictionary_file(
    custom_words,
    out = file.path(aspell_location, "words.pws")
  )
  
  invisible(custom_words)
}

author = c(
  "Balamuta",
  "Culpepper"
)
other = c(
  "Liang"
)

custom_words = c(author, other)

custom_user_spelling_dictionary(custom_words, "slcm")