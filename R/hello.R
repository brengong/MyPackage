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








