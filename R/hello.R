#' RevComp: Reverse Complement function
#'
#' This function takes the reverse, complement, or reverse complement of 1 or more DNA sequence strings.
#'
#' @param  x a DNA string, vector of DNA strings, or an object representing a vector of DNA strings.
#' @param  designates the type of output.  "comp", "rev", or "revcomp" designate taking the complement, reverse,
#' or reverse complement of the sequence.
#' @return a list of DNA strings from the desired DNA strings computed according to the selected functionn type as well as the original string.
#' @author Brendan Gongol
#' @export
#' @examples
#'hi <- c("TTGGGGGCGTTAG", "CATTGGTTCTAGTT", "GCATTGGTGTGCATTAG")
#'RevComp(myseq=hi, type= "rev")
#'RevComp(myseq=hi, type= "comp")
#'RevComp(myseq=hi, type= "revcomp")
#'RevComp(myseq="ATGCATTGGACGTTAG", type="rev")
#'RevComp(myseq="ATGCATTGGACGTTAG", type="comp")
#'RevComp(myseq="ATGCATTGGACGTTAG", type="revcomp")

RevComp <- function(myseq, type){

  RevComp2 <- function(myseq, type){ ### if you call RevComp() without calling any additional arguments, these values will be the defaults.
    if(type=="rev") {
      y <- substring(myseq, 1:nchar(myseq), 1:nchar(myseq)) ## nchar() returns the number of elements assigned to X.
      y <- rev(y)
      y <- paste(y, collapse="")
    }
    if(type=="comp") {
      y <- chartr("ATGC", "TACG", myseq) ##### chartr() is Translate character
    }
    if(type=="revcomp") {
      y <- substring(myseq, 1:nchar(myseq), 1:nchar(myseq))
      y <- (rev(y))
      y <- paste(y, collapse="")
      y <- chartr("ATGC", "TACG", y)
    }
    return(y)
  }


  if(type=="rev"){
    ret <- sapply(myseq, RevComp2, type="rev")
  }
  if(type=="comp"){
    ret <- sapply(myseq, RevComp2, type="comp")
  }
  if(type=="revcomp"){
    ret <-  sapply(myseq, RevComp2, type="revcomp")
  }
  return(ret)
}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' RevComp2: Reverse Complement function
#'
#' This function takes the reverse, complement, or reverse complement of 1 or more DNA sequence strings.
#'
#' @param  x a DNA string, vector of DNA strings, or an object representing a vector of DNA strings.
#' @param  designates the type of output.  "comp", "rev", or "revcomp" designate taking the complement, reverse,
#' or reverse complement of the sequence.
#' @return a list of DNA strings from the desired DNA strings computed according to the selected functionn type.
#' @author Brendan Gongol
#' @export
#' @examples
#'hi <- c("TTGGGGGCGTTAG", "CATTGGTTCTAGTT", "GCATTGGTGTGCATTAG")
#'RevComp2(myseq=hi, type= "rev")
#'RevComp2(myseq=hi, type= "comp")
#'RevComp2(myseq=hi, type= "revcomp")
#'RevComp2(myseq="ATGCATTGGACGTTAG", type="rev")
#'RevComp2(myseq="ATGCATTGGACGTTAG", type="comp")
#'RevComp2(myseq="ATGCATTGGACGTTAG", type="revcomp")

RevComp2 <- function(myseq, type){

  RevComp3 <- function(myseq2, type2){ ### if you call RevComp() without calling any additional arguments, these values will be the defaults.
    if(type2=="rev") {
      y <- substring(myseq2, 1:nchar(myseq2), 1:nchar(myseq2)) ## nchar() returns the number of elements assigned to X.
      y <- rev(y)
      y <- paste(y, collapse="")
    }
    if(type2=="comp") {
      y <- chartr("ATGC", "TACG", myseq2) ##### chartr() is Translate character
    }
    if(type2=="revcomp") {
      y <- substring(myseq2, 1:nchar(myseq2), 1:nchar(myseq2))
      y <- (rev(y))
      y <- paste(y, collapse="")
      y <- chartr("ATGC", "TACG", y)
    }
    return(y)
  }

  ret <- NULL
  if(type=="rev"){
    pb <- txtProgressBar(min = 0, max = length(myseq), style = 3)
    for(i in 1:length(myseq)){
      ret[i] <- RevComp3(myseq[i], type="rev")
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  if(type=="comp"){
    pb <- txtProgressBar(min = 0, max = length(myseq), style = 3)
    for(i in 1:length(myseq)){
      ret[i] <- RevComp3(myseq[i], type="comp")
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  if(type=="revcomp"){
    pb <- txtProgressBar(min = 0, max = length(myseq), style = 3)
    for(i in 1:length(myseq)){
      ret[i] <- RevComp3(myseq[i], type="revcomp")
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  return(ret)
}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################



#' Translate: Translate function.
#'
#' This function will translate one or many DNA sequences in all three reading frames into proteins.
#'
#' @param seq a vector of DNA strings.
#' @author Brendan Gongol
#' @export
#' @examples
#'hi <- c("TTGGGGGCGTTAG", "CATTGGTTCTAGTT", "GCATTGGTGTGCATTAG")
#'Translate(seq=hi)
#'seq <- c(s1="CTATATAGAGAGAGGAGAAATCTCT", s2="ACACACTGACGATTGTGACACAGATC")
#'Translate(seq=seq)

Translate <- function(seq=c(s1="CTATATAGAGAGA", s2="ACACACTGACGAT") ){

  AAv <- c("S", "S", "S", "S", "F", "F", "L", "L", "Y", "Y", "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R",
           "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
  names(AAv) <- c("TCA", "TCG", "TCC", "TCT", "TTT", "TTC", "TTA", "TTG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTA", "CTG", "CTC", "CTT", "CCA", "CCG", "CCC", "CCT", "CAT", "CAC", "CAA", "CAG", "CGA", "CGG", "CGC", "CGT",
                  "ATT", "ATC", "ATA", "ATG", "ACA", "ACG", "ACC", "ACT", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTA", "GTG", "GTC", "GTT", "GCA", "GCG", "GCC", "GCT", "GAT", "GAC", "GAA", "GAG", "GGA", "GGG", "GGC", "GGT")

  Translate2 <- function(x="ATGCATTGGACGTTAG"){
    seq <- NULL

    a <- gsub("(...)", "\\1_", x)
    a <- unlist(strsplit(a, "_"))
    a <- a[grep("^...$", a)]

    b <- gsub("^.", "", x) #^ anchors to the start of the string.  the dot represents one character.
    b <- gsub("(...)", "\\1_", b)
    b <- unlist(strsplit(b, "_"))
    b <- b[grep("^...$", b)]

    c <- gsub("^..", "", x)
    c <- gsub("(...)", "\\1_", c)
    c <- unlist(strsplit(c, "_"))
    c <- c[grep("^...$", c)]

    resultlist <- list(FR1=AAv[a], FR2=AAv[b], FR3=AAv[c])
    return(resultlist)
  }

  lapply(seq, Translate2)

}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' GOsorteR
#'
#' Requires a data.table with a column labeled "mergecol" and a vector of query terms
#'
#' @param DT a data table containing a column labeled "mergecol" housing the characters to query.
#' @param terms a vector of character strings indicating the terms to query "mergecol" for.
#' @return The data table that is subsetted according to the search terms listed in the "terms" argument.
#' @author Brendan Gongol
#' @importFrom data.table data.table
#' @export
#' @examples
#' library(data.table)
#' gene_symbol <- c("AKT", "PI3K", "SREBP", "FOXO", "ABCA1", "Caspase-1", "PIGPEN", "SNAIL")
#' chromosome_name <- c("1", "5", "10", "X", "Y", "2", "3", "4")
#' go_id <- c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004", "GO:0000005", "GO:0000006", "GO:0000007", "GO:0000008")
#' name_1006 <- c("biological_process", "cholesterol", "apoptosis", "differentiation", "Mitochondria Function", "cellular differentiation", "inflammation", "cell cycle")
#' definition_1006 <- c("Any process specifically pertinent to the functioning of integrated living units",
#'                      "Any process specifically pertinent to the functioning of cholesterol",
#'                      "Any process specifically pertinent to the functioning of apoptosis",
#'                      "Any process specifically pertinent to the functioning of differentiation",
#'                      "Any process specifically pertinent to the functioning of mitochondrial function",
#'                      "Any process specifically pertinent to the functioning of cellular differentiation",
#'                      "Any process specifically pertinent to the functioning of inflammation",
#'                      "Any process specifically pertinent to the functioning of cell cycle")
#' data <- data.table(cbind(gene_symbol, chromosome_name, go_id, name_1006, definition_1006))
#' data[, mergecol := paste(go_id,name_1006,definition_1006,sep="@")]
#'
#' words <- c("cholesterol", "differentiation", "functioning of mitochondria")
#'
#' GOsorteR(DT = data, terms = words)

GOsorteR <- function(DT, terms){

  Sorted_file <- NULL
  for(i in 1:length(terms)){
    sor <- DT[grepl(terms[i], mergecol)]
    Sorted_file <- rbind(Sorted_file, sor)
  }
  return(Sorted_file)
}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' BarcodescoreR
#'
#' Barcodes are used for multiplexing deep sequencing reactions.  However, in order to ensure the barcodes can be demultiplexed
#' efficiently, illumina requirea as close to a AC:TG ratio of 1:1 as possible for each position of the barcode.  Therefore, the
#' AC percentage of the barcode at a given position should be close to .5 across all barcodes.
#'
#' @param barcodes a character vector containing barcodes either 6 or 8 nucleotides in length
#' @param barlength Either 6 or 8 indicating the length of the barcodes to compute the AC percentage at each position for.
#' @return a vector containing the percentages at each position as well as a message indicating whether a barcode has been duplicated or not.
#' @author Brendan Gongol
#' @export
#' @examples
#' barcodes1 <- c("AACGCC","AAGGTA","AATTCC","ACACAG","ACACTC","ACACTG","ACAGGA","ACCTGT","ACGAAG","ACGACT","ACGTCA","CTGTTC",
#'                "GCCTAA","GATCTG","GTAGCC","GGAACT","GGACGG","GCGGAC","GGCCAC","GCTACC","GCTCAT","GTATAG","GTCGTC","GAATGA")
#' BarcodescoreR(barcodes1, barlength = 6)
#'
#' barcodes2 <- c("AACGCCAT","AAGGTACG","AATTCCGG","ACACAGAG","ACACTCAG","ACACTGTG","ACAGGACA","ACCTGTAG","ACGAAGGT","ACGACTTG",
#'                "ACGTCAAC","CTGTTCAC","GCCTAAGT","GATCTGGT","GTAGCCGT","GGAACTGT","GGACGGGT","GCGGACGT","GGCCACGT","GCTACCGT",
#'                "GCTCATGT","GTATAGGT","GTCGTCGT","GAATGAGT")
#' BarcodescoreR(barcodes2, barlength = 8)
#'
#' barcodes3 <- c("AACGCC","AAGGTA","AATTCC","ACACAG","ACACTC","ACACTG","ACAGGA","ACCTGT","ACGAAG","ACGACT","ACGTCA","CTGTTC",
#'                "GCCTAA","GATCTG","GTAGCC","GGAACT","GGACGG","GCGGAC","GGCCAC","GCTACC","GCTCAT","GTATAG","GTCGTC","AACGCC")
#' BarcodescoreR(barcodes3, barlength = 6)

BarcodescoreR <- function(barcodes, barlength){
  if(sum(duplicated(barcodes))>0){
    print("You have duplicated barcodes")
  }else{

    SP1 <- strsplit(barcodes, split = "")

    if(barlength == 8){
      AT1 <- NULL; AT2 <- NULL; AT3 <- NULL; AT4 <- NULL; AT5 <- NULL; AT6 <- NULL; AT7 <- NULL; AT8 <- NULL
      for(i in 1:length(SP1)){
        if(SP1[[i]][1] == "A" | SP1[[i]][1] == "C"){
          AT1[i] <- 1
        }else{
          AT1[i] <- 0
        }
        if(SP1[[i]][2] == "A" | SP1[[i]][2] == "C"){
          AT2[i] <- 1
        }else{
          AT2[i] <- 0
        }
        if(SP1[[i]][3] == "A" | SP1[[i]][3] == "C"){
          AT3[i] <- 1
        }else{
          AT3[i] <- 0
        }
        if(SP1[[i]][4] == "A" | SP1[[i]][4] == "C"){
          AT4[i] <- 1
        }else{
          AT4[i] <- 0
        }
        if(SP1[[i]][5] == "A" | SP1[[i]][5] == "C"){
          AT5[i] <- 1
        }else{
          AT5[i] <- 0
        }
        if(SP1[[i]][6] == "A" | SP1[[i]][6] == "C"){
          AT6[i] <- 1
        }else{
          AT6[i] <- 0
        }
        if(SP1[[i]][7] == "A" | SP1[[i]][7] == "C"){
          AT7[i] <- 1
        }else{
          AT7[i] <- 0
        }
        if(SP1[[i]][8] == "A" | SP1[[i]][8] == "C"){
          AT8[i] <- 1
        }else{
          AT8[i] <- 0
        }
      }
      print("No duplications found")
      return(c(sum(AT1)/length(SP1), sum(AT2)/length(SP1), sum(AT3)/length(SP1), sum(AT4)/length(SP1), sum(AT5)/length(SP1), sum(AT6)/length(SP1),
               sum(AT7)/length(SP1), sum(AT8)/length(SP1)))
    }

    if(barlength == 6){
      AT1 <- NULL; AT2 <- NULL; AT3 <- NULL; AT4 <- NULL; AT5 <- NULL; AT6 <- NULL
      for(i in 1:length(SP1)){
        if(SP1[[i]][1] == "A" | SP1[[i]][1] == "C"){
          AT1[i] <- 1
        }else{
          AT1[i] <- 0
        }
        if(SP1[[i]][2] == "A" | SP1[[i]][2] == "C"){
          AT2[i] <- 1
        }else{
          AT2[i] <- 0
        }
        if(SP1[[i]][3] == "A" | SP1[[i]][3] == "C"){
          AT3[i] <- 1
        }else{
          AT3[i] <- 0
        }
        if(SP1[[i]][4] == "A" | SP1[[i]][4] == "C"){
          AT4[i] <- 1
        }else{
          AT4[i] <- 0
        }
        if(SP1[[i]][5] == "A" | SP1[[i]][5] == "C"){
          AT5[i] <- 1
        }else{
          AT5[i] <- 0
        }
        if(SP1[[i]][6] == "A" | SP1[[i]][6] == "C"){
          AT6[i] <- 1
        }else{
          AT6[i] <- 0
        }
      }
      print("No duplications found")
      return(c(sum(AT1)/length(SP1), sum(AT2)/length(SP1), sum(AT3)/length(SP1), sum(AT4)/length(SP1), sum(AT5)/length(SP1), sum(AT6)/length(SP1)))
    }

  }
}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' BarcodeselectR
#'
#' Barcodes are used for multiplexing deep sequencing reactions.  However, in order to ensure the barcodes can be demultiplexed
#' efficiently, illumina requirea as close to a AC:TG ratio of 1:1 as possible for each position of the barcode.  Therefore, the
#' AC percentage of the barcode at a given position should be close to .5 across all barcodes.  This function returns the AC percentage
#' for all combinations of specified barcodes.
#'
#' @param barcodes A character vector containing barcodes either 6 or 8 nucleotides in length
#' @param number A numerical value indicating the number of barcodes the user wishes the function to return.
#' @param bar_length A numerical value, either 6 or 8, indicating the length of the barcodes to compute.
#' @return A data table.  Each row contains the combination of barcodes analyzed, the AC percentage at each position, and a numerical score of the overal barcode AC:TG diversity.
#' @author Brendan Gongol
#' @importFrom data.table data.table
#' @export
#' @examples
#' practice8 <- c("GTGTTCTA","AAACATCG","CATCAAGT","GTCTGTCA","AACAACCA","CATACCAA",
#'                "TGGAACAA","ACCTCCAA","ACACGACC","CTGAGCCA","ACACAGAA","TATCAGCA",
#'                "CTGTAGCC","TGGCTTCA","CCAGTTCA","AGATGTAC","ACAGCAGA","ACAAGCTA",
#'                "ACGCTCGA","CACTTCGA","AAGACGGA","AACGCTTA","GAATCTGA", "ATCCTGTA")
#' DT <- BarcodeselectR(barcodes = practice8, number = 23, barcode_length = 8)
#' DT <- DT[order(DT$score),]
#' DT
#'
#' practice6 <- c("GTGTTC","AAACAT","CATCAA","GTCTGT","AACAAC","CATACC",
#'                "TGGAAC","ACCTCC","ACACGA","CTGAGC","ACACAG","TATCAG",
#'                "CTGTAG","TGGCTT","CCAGTT","AGATGT","ACAGCA","ACAAGC",
#'                "ACGCTC","CACTTC","AAGACG","AACGCT","GAATCT", "ATCCTG")
#' DT <- BarcodeselectR(barcodes = practice6, number = 23, barcode_length = 6)
#' DT <- DT[order(DT$score),]
#' DT

BarcodeselectR <- function(barcodes, number, barcode_length){

  BarcodeR2 <- function(barcodes, barlength){
    if(sum(duplicated(barcodes))>0){
      print("You have duplicated barcodes")
    }else{

      SP1 <- strsplit(barcodes, split = "")

      if(barlength == 8){
        AT1 <- NULL; AT2 <- NULL; AT3 <- NULL; AT4 <- NULL; AT5 <- NULL; AT6 <- NULL; AT7 <- NULL; AT8 <- NULL
        for(i in 1:length(SP1)){
          if(SP1[[i]][1] == "A" | SP1[[i]][1] == "C"){
            AT1[i] <- 1
          }else{
            AT1[i] <- 0
          }
          if(SP1[[i]][2] == "A" | SP1[[i]][2] == "C"){
            AT2[i] <- 1
          }else{
            AT2[i] <- 0
          }
          if(SP1[[i]][3] == "A" | SP1[[i]][3] == "C"){
            AT3[i] <- 1
          }else{
            AT3[i] <- 0
          }
          if(SP1[[i]][4] == "A" | SP1[[i]][4] == "C"){
            AT4[i] <- 1
          }else{
            AT4[i] <- 0
          }
          if(SP1[[i]][5] == "A" | SP1[[i]][5] == "C"){
            AT5[i] <- 1
          }else{
            AT5[i] <- 0
          }
          if(SP1[[i]][6] == "A" | SP1[[i]][6] == "C"){
            AT6[i] <- 1
          }else{
            AT6[i] <- 0
          }
          if(SP1[[i]][7] == "A" | SP1[[i]][7] == "C"){
            AT7[i] <- 1
          }else{
            AT7[i] <- 0
          }
          if(SP1[[i]][8] == "A" | SP1[[i]][8] == "C"){
            AT8[i] <- 1
          }else{
            AT8[i] <- 0
          }
        }

        return(c(sum(AT1)/length(SP1), sum(AT2)/length(SP1), sum(AT3)/length(SP1), sum(AT4)/length(SP1), sum(AT5)/length(SP1), sum(AT6)/length(SP1),
                 sum(AT7)/length(SP1), sum(AT8)/length(SP1)))
      }

      if(barlength == 6){
        AT1 <- NULL; AT2 <- NULL; AT3 <- NULL; AT4 <- NULL; AT5 <- NULL; AT6 <- NULL
        for(i in 1:length(SP1)){
          if(SP1[[i]][1] == "A" | SP1[[i]][1] == "C"){
            AT1[i] <- 1
          }else{
            AT1[i] <- 0
          }
          if(SP1[[i]][2] == "A" | SP1[[i]][2] == "C"){
            AT2[i] <- 1
          }else{
            AT2[i] <- 0
          }
          if(SP1[[i]][3] == "A" | SP1[[i]][3] == "C"){
            AT3[i] <- 1
          }else{
            AT3[i] <- 0
          }
          if(SP1[[i]][4] == "A" | SP1[[i]][4] == "C"){
            AT4[i] <- 1
          }else{
            AT4[i] <- 0
          }
          if(SP1[[i]][5] == "A" | SP1[[i]][5] == "C"){
            AT5[i] <- 1
          }else{
            AT5[i] <- 0
          }
          if(SP1[[i]][6] == "A" | SP1[[i]][6] == "C"){
            AT6[i] <- 1
          }else{
            AT6[i] <- 0
          }
        }
        return(c(sum(AT1)/length(SP1), sum(AT2)/length(SP1), sum(AT3)/length(SP1), sum(AT4)/length(SP1), sum(AT5)/length(SP1), sum(AT6)/length(SP1)))
      }

    }
  }


  grid <- combn(barcodes, number)
  prac <- grid

  if(barcode_length == 6){
    pb <- txtProgressBar(min = 0, max = ncol(prac), style = 3)
    com <- NULL
    numb <- NULL
    for(i in 1:ncol(prac)){
      cal <- BarcodeR2(prac[,i], barlength = barcode_length)
      com <- rbind(com, cal)
      numb[i] <- sum(abs((com[i,]/.5) -1))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    rownames(com) <- 1:nrow(com)
    com <- data.table(com)
    com$score <- numb
    seq <- data.table(t(prac))
    setnames(com, c("V1", "V2", "V3", "V4", "V5", "V6"), c("AC1", "AC2", "AC3", "AC4", "AC5", "AC6"))
    DT <- cbind(seq, com)
    return(DT)
  }

  if(barcode_length == 8){
    pb <- txtProgressBar(min = 0, max = ncol(prac), style = 3)
    com <- NULL
    numb <- NULL
    for(i in 1:ncol(prac)){
      cal <- BarcodeR2(prac[,i], barlength = barcode_length)
      com <- rbind(com, cal)
      numb[i] <- sum(abs((com[i,]/.5) -1))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    rownames(com) <- 1:nrow(com)
    com <- data.table(com)
    com$score <- numb
    seq <- data.table(t(prac))
    setnames(com, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"), c("AC1", "AC2", "AC3", "AC4", "AC5", "AC6", "AC7", "AC8"))
    DT <- cbind(seq, com)
    return(DT)
  }

}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' file_sortR a file sorting function
#'
#' This function takes a directory of your choosing that contains everal files with similar names and
#' creates new directory containing subdirectories containing files named according to the subsetted file name specified,
#' and copies the files into the new directories.
#'
#'
#' @param directory path to the directory storing the files of interest. Must not have any spaces in the file path.
#' @param newdirectory a path containing the name of the new directory to make. Must not have any spaces in the file path.
#' @param files a numerical string indicating what files in the directory to move
#' @param Name_start a number indicating the start position of the file name to create a directory for.
#' @param Name_end a numer indicating the end position of the file name to create a directory for.
#' @param remove_old a TRUE/FALSE logical.  If TRUE, will delete the old files.
#' @return a new directory containing sub directories housing original files.
#' @author Brendan Gongol
#' @export
#' @examples
#' file_sorteR(directory = "C:/Users/Brendan/Desktop/AMPKKO2_1", newdirectory = "C:/Users/Brendan/Desktop/AMPK_practice",
#' files = 2:343, Name_start = 14, Name_end = 19, remove_old = TRUE)

file_sorteR <- function(directory, newdirectory, files, Name_start, Name_end, remove_old = FALSE){

  # directory <- "C:/Users/Brendan/Desktop/AMPKKO2_1"
  # newdirectory <- "C:/Users/Brendan/Desktop/AMPK_practice"
  # files = 2:343
  # Name_start = 14
  # Name_end = 19
  # remove_old = TRUE

  setwd(directory)
  fil <- dir()[files]
  # fil

  #### Obtain directory names ####
  dir_names <- NULL
  for(i in 1:length(fil)){
    spl <- strsplit(fil[i], split = "")[[1]][Name_start:Name_end]
    spl2 <- paste(spl, collapse = "")
    dir_names <- c(dir_names, spl2)
  }
  # dir_names

  #### Create a new directory and create all the files in the new directory ####
  dir.create(newdirectory)

  dupnames <- dir_names[!duplicated(dir_names)]
  for(i in 1:length(dupnames)){
    Dir <- paste(newdirectory, "/", dupnames[i], collapse = "")
    Dir <- gsub(" ", "", Dir)
    dir.create(Dir)
  }

  #### Identify the number of each file that will be move in to the new corresponding file.
  numbers <- NULL
  for(i in 1:length(dupnames)){
    numbers[i] <- sum(dir_names == dupnames[i])
  }
  # numbers

  #### move the files from one director to the new corresponding directory ####

  pb <- txtProgressBar(min = 0, max = length(numbers), style = 3)
  for(a in 1:length(numbers)){

    for(i in (sum(numbers[1:a])- numbers[a] + 1): sum(numbers[1:a])){
      fro <- paste(directory, "/", fil[i], collapse = "")
      fro <- gsub(" ", "", fro)
      TO <- paste(newdirectory, "/", dir_names[sum(numbers[1:a])- numbers[a] + 1], collapse = "")
      TO <- gsub(" ", "", TO)
      file.copy(from = fro, to = TO)
    }
    setTxtProgressBar(pb, a)
  }
  close(pb)

  if(remove_old == TRUE){
    file.remove(dir()[files])
  }

}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' FCcalculator
#'
#' This function takes two vectors and calculates the fold changes between them.
#'
#' @param vector1 a vector containing the first set of values
#' @param vector2 a vector containing the second set of values
#' @return a new vector containing the fold changes of the original values.
#' @author Brendan Gongol
#' @export
#' @examples
#'

FCcalculator <- function(vectorTX, vectorCTRL){
  FC <- NULL
  for(i in 1:length(vectorTX)){
    if(vectorTX[i] > vectorCTRL[i]){
      FC[i] <- vectorTX[i]/vectorCTRL[i]
    }else{
      FC[i] <- -(vectorTX[i]/vectorCTRL[i])
    }
  }
  return(FC)

}


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' demultiplex
#'
#' This function demultiplexes a fastq file into individual fastq files
#'
#' @param x fastq file to demultiplex
#' @param barcode vector containing barcodes to demultiplex
#' @param nreads number of records in fastq file to stream through in a given block
#' @param cutoff minimum read score cutoff
#' @return a new vector containing the fold changes of the original values.
#' @author Brendan Gongol
#' @export
#' @examples
#' demultiplex(x=fastq[1], barcode=c("TT", "AA", "GG"), nreads=50, cutoff=20)


demultiplex <- function(x, barcode, nreads, cutoff) {
  f <- FastqStreamer(x, nreads)
  while(length(fq <- yield(f))) {
    for(i in barcode) {
      pattern <- paste("^", i, sep="")
      fq <- trimTails(fq, k=2, a=rawToChar(as.raw(cutoff+33)), successive=FALSE)
      fqsub <- fq[grepl(pattern, sread(fq))]
      if(length(fqsub) > 0) {
        writeFastq(fqsub, paste(x, i, sep="_"), mode="a", compress=FALSE)
      }
    }
  }
  close(f)
}





##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################



#' ChiPseqPeakPlotter
#'
#' This function plotts the peak height of a CHiPseq experiment for a given chromosome location.
#'
#' @param treatment_bw a character vector indicating the location of the treatment bigwig file.
#' @param control_bw a character vector indicating the location of the treatment bigwig file.
#' @param chr a character vector indicating the chromosome number.
#' @param start a number indicating the start site.
#' @param end a number indicating the end site.
#' @param average_dist a number indicating the number of bases to average the peak height.
#' @param fill_dist a number indicating the number of average replicates to fill.
#' @param type a character either "all" or "called" indicating whether to include all peaks for a given area or only those called by peak caller.
#' @param peaks_tx a character vector indicating the location of the treatment MACS2 peak file.
#' @param peaks_ctrl a character vector indicating the location of the control MACS2 peak file.
#' @return a plot indicating the peak score for a specified chromosome position.
#' @author Brendan Gongol
#' @importFrom data.table data.table
#' @importFrom rtracklayer import
#' @export
#' @examples
#' ChiPseqPeakPlotter(treatment_bw <- "./example_datasets/8hr-PS1_S16_L008_R1_001.fastq.gz.hisat_PS_8hr.bw",
#'                    control_bw <- "./example_datasets/bigwig_files/8hr-OS1_S7_L008_R1_001.fastq.gz.hisat_OS_8hr.bw",
#'                    chr <- "10",
#'                    start <- 69266000,
#'                    end <- 69405881,
#'                    average_dist <- 50,
#'                    fill_dist <- 50,
#'                    type = "all")
#'
#' ChiPseqPeakPlotter(treatment_bw <- "./example_datasets/8hr-PS1_S16_L008_R1_001.fastq.gz.hisat_PS_8hr.bw",
#'                    control_bw <- "./example_datasets/8hr-OS1_S7_L008_R1_001.fastq.gz.hisat_OS_8hr.bw",
#'                    chr <- "10",
#'                    start <- 69266000,
#'                    end <- 69405881,
#'                    average_dist <- 50,
#'                    fill_dist <- 50,
#'                    type = "called",
#'                    peaks_tx <- "./example_datasets/8hr-PS1_S16_L008_R1_001.fastq.gz.hisat_PS_8hr.bam_macs2_peaks.annotatedpeakAnno.xls",
#'                    peaks_ctrl <- "./example_datasets/8hr-OS1_S7_L008_R1_001.fastq.gz.hisat_OS_8hr.bam_macs2_peaks.annotatedpeakAnno.xls")


ChiPseqPeakPlotter <- function(treatment_bw, control_bw, chr, start, end, average_dist, fill_dist, type="all", peaks_tx, peaks_ctrl){

  if(type == "all"){
    ##############################################
    ######### Code to get score from bw files ####
    ##############################################
    which <- GRanges(c(chr), IRanges(c(start), c(end)))
    bw_tr <- import(treatment_bw, which = which)
    bw_ct <- import(control_bw, which = which)

    #### format treatment data ####
    ###############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_tr)))){
      num <- start(ranges(bw_tr))[i]:end(ranges(bw_tr))[i]
      sco <- score(bw_tr)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    len <- data.table(start:end); setnames(len, colnames(len), "num")
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)

    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))######################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_tx <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_tx <- rbind(data_tx, d)
    }
    data_tx$num <- 1:nrow(data_tx)

    #### format control data ####
    #############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_ct)))){
      num <- start(ranges(bw_ct))[i]:end(ranges(bw_ct))[i]
      sco <- score(bw_ct)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)

    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))####################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_ctrl <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_ctrl <- rbind(data_ctrl, d)
    }
    data_ctrl$num <- 1:nrow(data_ctrl)

    #### combine TX and CTRL data together ####
    ###########################################
    data_tx$fill <- "tx"
    data_ctrl$fill <- "ctrl"
    dt <- rbind(data_tx, data_ctrl)

    #### plot data ####
    ###################
    p <- ggplot(aes(x=num, y=mea), data = dt)+
      geom_line(aes(color = fill))+   #color = c("#1B9E77"))+
      scale_color_brewer(palette = "Dark2")+ # "Set1"

      #### add fill below line and color ####
    geom_area(aes(fill=fill))+
      scale_fill_brewer(palette="Dark2")+ # "Dark2"  "Set1"
      # scale_fill_manual(values=c("#1B9E77", "#7570B3"))+ # Manually fill colors
      # scale_color_manual(values=c("#1B9E77", "#7570B3")) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(p)
  }

  else if(type == "called"){

    #########################################
    ### plot only peaks that were called ####
    #########################################

    #### get peaks for treatment and control files ####
    ###################################################
    df_tx <- read.delim(peaks_tx, comment = "#")
    df_tx$seqnames <- gsub("chr", "", df_tx$seqnames, ignore.case = TRUE)
    peaks_tx <- as(df_tx, "GRanges")

    df_ctrl <- read.delim(peaks_ctrl, comment = "#")
    df_ctrl$seqnames <- gsub("chr", "", df_ctrl$seqnames, ignore.case = TRUE)
    peaks_ctrl <- as(df_ctrl, "GRanges")

    #### subset out peaks for region of interest ####
    #################################################
    which <- GRanges(c(chr), IRanges(c(start), c(end)))

    which_tx <- subsetByOverlaps(peaks_tx, which)
    which_ctrl <- subsetByOverlaps(peaks_ctrl, which)

    #### get score for selected regions ####
    ########################################
    bw_tr <- suppressWarnings(import(treatment_bw, which = which_tx))########################
    bw_ct <- suppressWarnings(import(control_bw, which = which_ctrl))########################

    #### format treatment data ####
    ###############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_tr)))){
      num <- start(ranges(bw_tr))[i]:end(ranges(bw_tr))[i]
      sco <- score(bw_tr)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    len <- data.table(start:end); setnames(len, colnames(len), "num")
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)

    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))#####################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_tx <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_tx <- rbind(data_tx, d)
    }
    data_tx$num <- 1:nrow(data_tx)

    #### format control data ####
    #############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_ct)))){
      num <- start(ranges(bw_ct))[i]:end(ranges(bw_ct))[i]
      sco <- score(bw_ct)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)

    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))#####################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_ctrl <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_ctrl <- rbind(data_ctrl, d)
    }
    data_ctrl$num <- 1:nrow(data_ctrl)

    #### combine TX and CTRL data together ####
    ###########################################
    data_tx$fill <- "tx"
    data_ctrl$fill <- "ctrl"
    dt <- rbind(data_tx, data_ctrl)

    #### plot data ####
    ###################
    p <- ggplot(aes(x=num, y=mea), data = dt)+
      geom_line(color = "#1B9E77")+
      geom_area(aes(fill=fill))+
      scale_fill_brewer(palette="Dark2")+
      # scale_fill_manual(values=c("#1B9E77", "#7570B3"))+ # Manually fill colors
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(p)
  }

}



##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' ChiPseqregionscoreR
#'
#' This function returns the average peak height of a CHiPseq experiment for a given chromosome location.
#'
#' @param treatment_bw a character vector indicating the location of the treatment bigwig file.
#' @param control_bw a character vector indicating the location of the treatment bigwig file.
#' @param chr a character vector indicating the chromosome number.
#' @param start a number indicating the start site.
#' @param end a number indicating the end site.
#' @param average_dist a number indicating the number of bases to average the peak height.
#' @return a data.table of the average peak score for a specified chromosome region.
#' @author Brendan Gongol
#' @importFrom data.table data.table
#' @importFrom rtracklayer import
#' @export
#' @examples
#' dt <- ChiPseqregionscoreR(treatment_bw <- "./results/bigwig_files/2hr-PS1_S13_L008_R1_001.fastq.gz.hisat_PS_2hr.bw",
#'                           control_bw <- "./results/bigwig_files/2hr-OS1_S4_L008_R1_001.fastq.gz.hisat_OS_2hr.bw",
#'                           chr <- "12",
#'                           start <- 48105253 - 5000, # 48101253
#'                           end <- 48105253 + 5000, # 48150404
#'                           average_dist <- 50)
#' hr2 <- data.table(dt[dt$fill == "tx",]$mea/dt[dt$fill == "ctrl",]$mea)
#' hr2$time = "hr2"
#'
#' dt <- ChiPseqregionscoreR(treatment_bw <- "./results/bigwig_files/8hr-PS1_S16_L008_R1_001.fastq.gz.hisat_PS_8hr.bw",
#'                    control_bw <- "./results/bigwig_files/8hr-OS1_S7_L008_R1_001.fastq.gz.hisat_OS_8hr.bw",
#'                    chr <- "12",
#'                    start <- 48105253 - 5000, # 48101253
#'                    end <- 48105253 + 5000, # 48150404
#'                    average_dist <- 50)
#' hr8 <- data.table(dt[dt$fill == "tx",]$mea/dt[dt$fill == "ctrl",]$mea)
#' hr8$time = "hr8"
#'
#' dt <- ChiPseqregionscoreR(treatment_bw <- "./results/bigwig_files/16hr-OS1_S10_L008_R1_001.fastq.gz.hisat_OS_16hr.bw",
#'                           control_bw <- "./results/bigwig_files/16hr-PS1_S19_L008_R1_001.fastq.gz.hisat_PS_16hr.bw",
#'                           chr <- "12",
#'                           start <- 48105253 - 5000, # 48101253
#'                           end <- 48105253 + 5000, # 48150404
#'                           average_dist <- 50)
#' hr16 <- data.table(dt[dt$fill == "tx",]$mea/dt[dt$fill == "ctrl",]$mea)
#' hr16$time = "hr16"
#'
#' dt2 <- rbind(hr2, hr8)
#' dt2 <- rbind(dt2, hr16); setnames(dt2, colnames(dt2), c("mea", "fill"))
#' dt2 <- dt2[dt2$mea > 0,]
#' dt2$mea <- log2(dt2$mea)
#' dt2 <- dt2[!dt2$mea == "-Inf",]
#' dt2 <- dt2[!dt2$mea == "Inf",]
#' dt2


ChiPseqregionscoreR <- function(treatment_bw, control_bw, chr, start, end, average_dist){

  ##############################################
  ######### Code to get score from bw files ####
  ##############################################
  which <- GRanges(c(chr), IRanges(c(start), c(end)))
  bw_tr <- import(treatment_bw, which = which)
  bw_ct <- import(control_bw, which = which)

  #### format treatment data ####
  ###############################
  scores <- data.table()
  for(i in 1:length(start(ranges(bw_tr)))){
    num <- start(ranges(bw_tr))[i]:end(ranges(bw_tr))[i]
    sco <- score(bw_tr)[i]
    dt <- data.table(num, sco)
    scores <- rbind(scores, dt)
  }
  scores
  len <- data.table(start:end); setnames(len, colnames(len), "num")
  zz <- merge(len, scores, all = TRUE)
  zz[is.na(zz)] <- 0
  zz$key <- 1:nrow(zz)

  #### create averages across specified windows ####
  ##################################################
  mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))######################
  data_tx <- data.table(mea)
  # #### fill in lost rows ####
  # data_tx <- data.table()
  # for(i in 1:nrow(dt)){
  #   d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
  #   data_tx <- rbind(data_tx, d)
  # }
  data_tx$num <- 1:nrow(data_tx)

  #### format control data ####
  #############################
  scores <- data.table()
  for(i in 1:length(start(ranges(bw_ct)))){
    num <- start(ranges(bw_ct))[i]:end(ranges(bw_ct))[i]
    sco <- score(bw_ct)[i]
    dt <- data.table(num, sco)
    scores <- rbind(scores, dt)
  }
  scores
  zz <- merge(len, scores, all = TRUE)
  zz[is.na(zz)] <- 0
  zz$key <- 1:nrow(zz)

  #### create averages across specified windows ####
  ##################################################
  mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))####################
  data_ctrl <- data.table(mea)
  # #### fill in lost rows ####
  # data_ctrl <- data.table()
  # for(i in 1:nrow(dt)){
  #   d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
  #   data_ctrl <- rbind(data_ctrl, d)
  # }
  data_ctrl$num <- 1:nrow(data_ctrl)

  #### combine TX and CTRL data together ####
  ###########################################
  data_tx$fill <- "tx"
  data_ctrl$fill <- "ctrl"
  dt <- rbind(data_tx, data_ctrl)
  return(dt)

}



##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


#' multiplot
#'
#' Plot multiple plots in a single window
#'
#' @param file
#' @param cols
#' @return a plot containing all figures
#' @author Brendan Gongol
#' @export
#' @examples
#' multiplot(p1, p2, p3, cols = 1)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}












