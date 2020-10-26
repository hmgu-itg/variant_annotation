searchAtlasBaselineExperiments <- function (species = NULL)
{
    if (missing(species)) {
        message("No species was provided. Will search for data from all available species.")
    }
    else if (typeof(species) != "character") {
        stop("Please provide species as a character vector.")
    }
    else if (length(species) > 1) {
        stop("More than one species found. You may only specify one species at a time.")
    }
    aeAPIbase <- "http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?experimentType=baseline"
    queryURL <- paste(aeAPIbase,"&gxa=TRUE", sep = "")
    if (!missing(species)) {
        species <- URLencode(species)
        queryURL <- paste(queryURL, "&species=", species, sep = "")
    }
    message("Searching for Expression Atlas experiments matching your query ...")
    response <- httr::GET(queryURL)
    if (httr::status_code(response) != 200) {
        stop(paste("Error running query. Received HTTP error code",
            status_code(response), "from server. Please try again later. If you continue to experience problems please email atlas-feedback\
@ebi.ac.uk"))
    }
    else {
        message("Query successful.")
    }
    parsedXML <- XML::xmlParse(httr::content(response))
    allExpsNode <- XML::xmlRoot(parsedXML)
    numExps <- XML::xmlAttrs(allExpsNode)["total"]
    if (numExps == 0) {
        return(message("No results found. Cannot continue."))
    }
    else {
        message(paste("Found", numExps, "experiments matching your query."))
    }
    allExperiments <- XML::xmlElementsByTagName(allExpsNode, "experiment")
    resultsList <- lapply(allExperiments, function(experimentNode) {
        expAcc <- XML::xmlValue(XML::xmlElementsByTagName(experimentNode,
            "accession")$accession)
        expTitle <- XML::xmlValue(XML::xmlElementsByTagName(experimentNode,
            "name")$name)
        species <- XML::xmlValue(XML::xmlElementsByTagName(experimentNode,
            "organism")$organism)
        expType <- XML::xmlValue(XML::xmlElementsByTagName(experimentNode,
            "experimenttype")$experimenttype)
        list(accession = expAcc, title = expTitle, species = species,
            expType = expType)
    })
    allAccessions <- sapply(resultsList, function(x) {
        x$accession
    })
    allExpTypes <- sapply(resultsList, function(x) {
        x$expType
    })
    allSpecies <- sapply(resultsList, function(x) {
        x$species
    })
    allTitles <- sapply(resultsList, function(x) {
        x$title
    })
    names(allAccessions) <- NULL
    names(allExpTypes) <- NULL
    names(allSpecies) <- NULL
    names(allTitles) <- NULL
    resultsSummary <- DataFrame(Accession = allAccessions, Species = allSpecies,
        Type = allExpTypes, Title = allTitles)
    resultsSummary <- resultsSummary[order(resultsSummary$Species,
        resultsSummary$Type, resultsSummary$Accession), ]
    return(resultsSummary)
}
