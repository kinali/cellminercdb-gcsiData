# ============================================================
# Add PubChem IDs to a drug synonym table
#
# Requirements:
#   install.packages("remotes")
#   remotes::install_github("cannin/slugify")
#   remotes::install_github("cannin/rcellminerPubchem")
#
# Expected input columns:
#   NAME_SET
#   nci60
#   PUBCHEM_ID
#
# You can check this code with a small data
# test_table <- drugSynonymTab_pubchem[1:5, ]
# test_result <- add_pubchem_ids(test_table)

# ============================================================

library(rcellminerPubchem)

cacheDir <- "cacheDir"

# Make sure the cache directory exists
if (!dir.exists(cacheDir)) {
  dir.create(cacheDir, recursive = TRUE)
}


# extract the first value from the NAME_SET column
# NAME_SET contains multiple drug synonyms, so return one
get_first_value <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  
  as.character(x[[1]])
}


# Give drug_table = drugSynonymTab_pubchem and it will return the same table with PubChem IDs added

add_pubchem_ids <- function(drug_table) {
  
  
  # Check required columns
    # NAME_SET: used to search PubChem by drug name
    # nci60: used to search PubChem by NSC number
    # PUBCHEM_ID: existing PubChem ID, if one is already available
  required_cols <- c(
    "NAME_SET",
    "nci60",
    "PUBCHEM_ID"
  )
  
  # If at least one required column is missing stop the script and tell the user which column is missing.
  missing_cols <- setdiff(
    required_cols,
    colnames(drug_table)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # 1. Extract primary drug name
  # ----------------------------------------------------------
  
  drug_table$primary_name <- vapply(
    drug_table$NAME_SET,
    get_first_value,
    character(1)
  )
  
  
  # ----------------------------------------------------------
  # 2. Extract NCI-60 NSC identifier
  # ----------------------------------------------------------
  
  # This goes through every entry in NAME_SET and applies get_first_value()
  
  drug_table$nsc_id <- vapply(
    drug_table$nci60,
    get_first_value,
    character(1)
  )
  
  
  # ----------------------------------------------------------
  # 3. Search PubChem using drug names [1] and [2] from NAME_SET
  # ----------------------------------------------------------

  
  drug_table$PUBCHEM_from_name <- vapply(
    drug_table$NAME_SET,
    function(x) {
      
      # Check the first two names
      for (name in head(x, 2)) {
        
        if (is.na(name) || name == "") {
          next
        }
        
        result <- as.character(getPubchemFromName(name))
        
        # Return the first valid PubChem ID
        if (length(result) == 1 &&
            !is.na(result) &&
            result != "") {
          return(result)
        }
      }
      
      # If neither name returns a valid PubChem ID
      return(NA_character_)
    },
    character(1)
  )
  
  
  # ----------------------------------------------------------
  # 4. Search PubChem using NSC number
  # ----------------------------------------------------------
  
  drug_table$PUBCHEM_from_nsc <- vapply(
    drug_table$nsc_id,
    function(x) {
      
      if (is.na(x) || x == "") {
        return(NA_character_)
      }
      
      as.character(getPubchemFromNsc(x))
      
    },
    character(1)
  )
  
  
  
  # ----------------------------------------------------------
  # 5. Combine PubChem IDs
  #
  # Priority:
  #   1. PubChem ID obtained from drug name
  #   2. Existing PUBCHEM_ID
  #   3. PubChem ID obtained from NSC
  # ----------------------------------------------------------
  
  drug_table$PUBCHEM_ID_updated <- drug_table$PUBCHEM_from_name
  
  missing <- is.na(drug_table$PUBCHEM_ID_updated)
  
  drug_table$PUBCHEM_ID_updated[missing] <-
    as.character(drug_table$PUBCHEM_ID[missing])
  
  
  missing <- is.na(drug_table$PUBCHEM_ID_updated)
  
  drug_table$PUBCHEM_ID_updated[missing] <-
    drug_table$PUBCHEM_from_nsc[missing]
  
  
  return(drug_table)
}


test_table <- drugSynonymTab_pubchem[1:5, ]
test_result <- add_pubchem_ids(test_table)

