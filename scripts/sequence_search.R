# ----------------------------------
#  Search for ITS2 sequences
# ----------------------------------

library(rentrez)
library(taxize)
library(taxizedb)
library(glue)
library(purrr)
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(Biostrings)

organism = snakemake@config$organism
marker = snakemake@config$marker

seq_outfile = snakemake@output$seqs
taxa_outfile = snakemake@output$taxa

download_taxdb = snakemake@config$download_taxdb
delete_taxdb = snakemake@config$delete_taxdb
taxdb_cache = snakemake@params$taxdb_cache

# make sure we've got a supported marker
if (!marker %in% c("ITS2", "18S")) {
	stop(marker, " not supported")
}

search_text = list(
	"ITS2" = "(ITS2 OR internal transcribed spacer 2)",
	"18S" = "(18S OR SSU rRNA OR SSU OR rRNA)"
)

term = glue("{organism}[ORGN] {search_text[[marker]]}")

# search nucleotide using provided organism -------------------------------

run_search = function(term) {
	message("Search NCBI with the following search term ", term)
	search = entrez_search(db = "nucleotide", term = term, use_history = TRUE)
	message("Found ", search$count, " hits")
	return(search)
}

# Get summaries -----------------------------------------------------------
# needed so we can add taxonomy
get_summaries = function(search) {
	message("Getting summaries...")
	summaries = seq(0, search$count, by = 200) %>%
	  map(~entrez_summary(db = "nucleotide",
	                      web_history = search$web_history,
	                      retmax = 200,
	                      retstart = .x)) %>%
	  flatten()
}

# turn search summaries into a dataframe
summaries_to_dataframe = function(summaries) {
	map_dfr(summaries,
					~tibble(
						uid = as.character(.x[["uid"]]),
						taxid = as.character(.x[["taxid"]])))
}

# add accession numbers
add_acc = function(sum_df) {
	ids = split(sum_df$uid, ceiling(seq_along(sum_df$uid) / 200))
	accns = map(ids, ~entrez_fetch("nucleotide", .x, rettype = "acc"))
	accns = str_c(accns, collapse = "\n")
	accns = str_split(accns, "\n")[[1]]
	accns = discard(accns, ~.x == "")

	if (length(accns) != nrow(sum_df)) {
	  stop("Number of accession ids retrieved does not match the number of search hits")
	}

	return(mutate(sum_df, accn = accns))
}


# Get taxonomic classifications -------------------------------------------

add_taxonomy = function(summaries, sum_df, download_taxdb, delete_taxdb, taxdb_cache) {
	message("Retrieving taxonomic classifications...")

	taxids = map_int(summaries, "taxid") %>% unique() %>% as.character()
	ranks = c("superkingdom", "kingdom", "phylum", "class", "order",
	          "family", "genus", "species" )

	taxize::taxize_options(ncbi_sleep = 0.350) # prevent reaching rate limit

	if (download_taxdb) {
	  taxizedb::tdb_cache$cache_path_set(full_path = taxdb_cache)
	  taxizedb::db_download_ncbi() # downloads to "taxdb_cache" dir
	  taxizedb::src_ncbi()
	  tax = map(setNames(taxids, nm = taxids), ~ retry_classification(.x, download_taxdb))

	  # delete taxonomic database cache
	  if (delete_taxdb) {
	    message("Deleting downloaded tax DB...")
	    taxizedb::tdb_cache$delete_all()
    }
	} else {
	  tax = map(setNames(taxids, nm = taxids), ~ retry_classification(.x, download_taxdb))
	}

	tax = tax %>%
	  bind_rows(.id = "taxid") %>%
	  select(-id) %>%
	  filter(rank %in% ranks) %>%
	  spread(rank, name) %>%
	  select(taxid, all_of(ranks))

	# add in the taxonomies and write it out
	sum_df = sum_df %>% left_join(tax, by = "taxid")


	# filter out bad taxonomies, species without proper identification,
	# contains sp. and numbers.  Also filters out taxa that are NA
	sum_df = sum_df %>%
		filter(!str_detect(sum_df$species, "sp\\.|\\d+"),
					 !is.na(genus),
					 !is.na(family),
					 !is.na(order),
					 !is.na(class),
					 !is.na(phylum),
					 !is.na(kingdom),
					 !is.na(superkingdom))


	return(sum_df)

}


# give NULL value instead of stopping when HTTP error occurs
safe_classification = possibly(
  ~ taxize::classification(.x, db = "ncbi", max_tries = 10)[[1]], # [[1]] to extract dataframe only without attributes
  otherwise = NULL
  )

# retry fetching taxonomy if error occurs
retry_classification = function(taxid, download_taxdb) {

  if (download_taxdb) {
    tax = taxizedb::classification(taxid, db = "ncbi")[[1]]
    # check if taxonomy is fetched, i.e. returns a dataframe
    tax_exist = is.data.frame(tax)
    if (tax_exist) {
      remote_taxdb = FALSE
    } else {
      # fetch missing tax from NCBI remotely
      # NCBI taxonomy dump file doesn't contain all tax id
      # e.g. species with sp. and numbers
      remote_taxdb = TRUE
      message(cat("Tax ID", taxid, "not found in downloaded tax DB, fetching from NCBI remotely..."))
    }
  } else {
    remote_taxdb = TRUE
  }

  while (remote_taxdb == TRUE) {
    tax = safe_classification(taxid)
    tax_exist = is.data.frame(tax)
    if (tax_exist) {remote_taxdb = FALSE}
  }

  tax
}


# Get fasta sequences -----------------------------------------------------

retrieve_and_write_sequences = function(ids, outfile, marker) {
	message("Retrieving sequences...")
	ids = split(ids, ceiling(seq_along(ids) / 50))
	seqs = map(ids, ~entrez_fetch("nucleotide",
																id = .x,
																rettype = "fasta",
																retmode = "text")) %>%
		flatten()

	message("Writing fasta of hits")
	if (file.exists(outfile)) {
		warning(outfile, " exists, removing..")
		file.remove(outfile)
	}

	seqs %>%
  	walk(~write_file(.x, outfile, append = TRUE))


  if (marker == "ITS2") {

    # add additional filtering to remove false positive hits to ITS1 - not sure why these are there
    raw = readDNAStringSet(outfile)
    reg = "(?:ITS2)|(?:ITS-2)|([Ii]nternal transcribed spacer 2)|(second internal transcribed spacer)|(internal transcribed spacers 1 and 2)"
    keeps = str_detect(names(raw), reg)
    raw = raw[keeps]

    if (length(raw) == 0) {
      stop("Problem filtering out mismatched ITS2 - bad problem - please contact developers.")
    }

    writeXStringSet(raw, outfile)
  }
}



# Run it ------------------------------------------------------------------

search = run_search(term)
summaries = get_summaries(search)
sum_df = summaries_to_dataframe(summaries)
sum_df = add_acc(sum_df)
sum_df = add_taxonomy(summaries, sum_df, download_taxdb, delete_taxdb, taxdb_cache)
retrieve_and_write_sequences(sum_df$uid, seq_outfile, marker)

message("Writing taxonomy info")
write_tsv(sum_df, taxa_outfile)



