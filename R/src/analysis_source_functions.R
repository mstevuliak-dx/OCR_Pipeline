# ======================================= FUNCTIONS USED IN ANALYSIS
# Author: Matej Stevuliak


# ------------------- BACKGROUND CORRECTION

correct_background <- function(d){
  # correct for outliars in background measurements and remove background
  # takes wells marked as Backfround and for each plate and time substracts median of those from OCR measurements

  background <- d %>%
    filter(Protocol == "Background") %>%
    group_by(Measurement, plate_id) %>%
    mutate( Z_mod = abs(0.6745*(OCR - median(OCR))/mad(OCR)) ) %>%
    mutate(Mean_to_substr = mean(OCR)) # take mean to substract

  out        <- background %>% filter(Z_mod > 3.5)
  background <- background %>% filter(Z_mod < 3.5)
  background <- background[!duplicated(background[,"Mean_to_substr"]),] # reduce the data before joining them

  d <- d %>%
    left_join(background, by = c("Measurement", "plate_id" ), suffix  = c("",".y") ) %>%
    mutate(OCR = OCR - Mean_to_substr ) %>%
    select(-c(contains("y"),"Z_mod", "Mean_to_substr")) %>%
    filter(Protocol != "Background" )

  return(list(data = d, removed = out))
  rm(background,d)
}

# --------------------------------------------------------------------------------------CONSTRUCT DATA

read_xlsx_set <- function(path_, pattern_){
  # returns a combined dataframe from all .xlsx files in folder (path_).
  # can specify pattern to distinguish "pre" "post"  treatment as: pattern_ = "*pre/post.xlsx"

  get_intervals <- function(x, y){
    # Identifies Interval form measurement time  and group info Used with mapply
    return(
      ifelse(grepl("Background",y), "Background",
             ifelse(grepl(as.character(x), "123"), "Int1",
                    ifelse((grepl(as.character(x), "456") & (grepl("Oligo", y)|grepl("Glyco", y))), "Int2",
                           ifelse((grepl(as.character(x), "456") & grepl("FCCP", y)), "Int3",
                                  ifelse(grepl(as.character(x), "789") & grepl("Glyco", y) , "Int5",
                                         ifelse(grepl(as.character(x), "789") & ! grepl("Glyco", y), "Int4", "Other" )))))))
  }

  # list of filenames
  files    <- list.files(path=path_, pattern=pattern_, full.names=TRUE, recursive=FALSE)
  merged_d <- data_frame()

  # for every file create a data frame and merge into one
  i <- 0   # to tag sample_id
  list.data_frames <- lapply(files, function(x){
    # Read rates
    d <- read_excel(x,  sheet = "Rate")

    # read assay configuration
    head       <- read_excel(x,  sheet = "Assay Configuration")
    assay_name <- as.character(head[3,2])
    rm(head)

    # # Read mmHg from raw data
    # raw <- read_excel(x,  sheet = "Raw")
    # raw <- raw %>% filter(Tick == "0") %>% select(c("O2 (mmHg)","Well", "Measurement") )
    #
    # d <- d %>% left_join(raw, by = c("Measurement", "Well"))
    # rm(raw)
    #
    # add plateID column
    d$plate_id <- rep(assay_name, times = nrow(d))

    # add Protocol column
    d <- d %>%
      arrange(Well) %>%
      mutate(Protocol = ifelse(! grepl("\\+", as.character(x)) & !Well %in% c( "A12", "H12", "D07") & 6<as.numeric(substr(Well,2,3)), NA,
                          ifelse(grepl("Oligo", Group), "Olig",
                            ifelse(grepl("FCCP",Group), "FCCP",
                              ifelse(grepl("Glyco", Group), "Glyco",
                                ifelse(grepl("Backgr",Group), "Background", "Other") )))))

    # remove Background from OCR data
    if (mean(filter(d, Protocol == "Background")$OCR) == 0 ) {
      d <- d %>% filter(Protocol != "Background")
    } else {
      corr <- correct_background(d)
      d    <- corr$data
      out  <- corr$removed
      cat("There were ", nrow(out), " background measurements removed on plate: ", assay_name, "\n")
    }


    # add interval column
    d$Interval <- mapply(get_intervals, d$Measurement, d$Group)

    # add column with time on experiment scale
    d <- d %>% mutate(Time = ifelse(Interval == "Int1" | Interval == "Int2", Measurement,
                              ifelse(Interval == "Int3" | Interval == "Int4", Measurement + 3, Measurement + 6)))

    # the Enries with blanks should be filtered, Must be here !
    d <- d %>%filter(!is.na(Protocol))

    # add log OCR
    d$LOCR <- sapply(d$OCR, log)

    # add sample ID   If "+" in the filename recognizes two samples on plate
    if (grepl("\\+", as.character(x))) {
      d$sample_id <- as.character(sapply(d$Group, function(x){ substr(x,(nchar(x)+1)-7,nchar(x)) }))
    }else{
      d$sample_id <- rep(file_path_sans_ext(basename(x)), nrow(d))  # Takes sample name from filename
    }

    d <- d %>% mutate(sample_id = paste0(sample_id, "-", i))
    i <<- i+1

    return(d)
  })
  merged_d <- do.call(rbind, list.data_frames)

  return(merged_d)
}

# -------------------------------------------------------------------------------------- IDENTITY OUTLIARS

idfy_outliar <- function(DT, x, cut.well, cut.point ){

  # if statements ...
  # if (x == "ECAR") {
  #   x == "LECAR"
  # }


  dm <- DT
  dm$x <- dm[[x]] # Create column with regressed variable
  i <- 1
  keep <- T
  size <- nrow(dm)
  dm <- dm %>% drop_na()
  while (keep) {
    # --- Fitting linear regresion, separate wells in separate intervals (and sample)
    fit_all <- lm(x ~ -1 + sample_id + Interval + Well, data = dm) # can be one sample_ID
    dm$fitted <- fitted(fit_all)
    dm$residuals = residuals(fit_all)

    dm <- dm %>%  group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
             n= n(),
             sq_err   = (x-int_mean)^2) %>%
      ungroup()
    # get mean error in well
    dm <- dm %>%  group_by(sample_id, Well) %>%
      mutate(mean_err_sq = mean(sq_err, na.rm = T) ) %>%
      ungroup()
    # Done for Well
    dm <- dm %>% group_by(sample_id) %>%
      mutate(median_mean_sq_err = median(mean_err_sq, na.rm = T ),
             mad_mean_sqE       = mad(mean_err_sq, na.rm = T),
             is.out.w           = median_mean_sq_err + cut.well * mad_mean_sqE < mean_err_sq) %>%
      ungroup()
    # Print
    n_out_w <-  nrow(dm %>% filter(is.out.w == TRUE))
    cat(i, " Well outliares: ", n_out_w, "--", n_out_w/size*100 ,  "% \n")

    dm <- dm %>% filter(is.out.w == F)
    if( n_out_w == 0) keep <-FALSE # Stop when no more outliars found
    i <- i+1
  }
  # complete data
  dm_r <- DT %>%  drop_na() %>%
    left_join(dm, by = c("Measurement", "Well", "Group", "Time", "plate_id", "Protocol", "Interval"), suffix  = c("",".y") ) %>%
    mutate(is.out.w = replace_na(is.out.w, T)) %>%
    select(-c(contains("y"), "mad_mean_sqE", "sq_err", "int_mean"))

  # ------------ Single Point Outliars

  i <- 1
  keep <- T
  size <- nrow(dm)
  dm <- dm %>% drop_na()
  while (keep) {

    fit_all <- lm(x ~ -1 + sample_id + Interval + Well, data = dm)
    dm$fitted <- fitted(fit_all)
    dm$residuals = residuals(fit_all)

    dm <- dm %>% group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
             sq_err   = (x-int_mean)^2) %>%
      ungroup()
    # get mean error in well
    dm <- dm %>% group_by(sample_id, Well) %>%
      mutate(mean_err_sq = mean(sq_err, na.rm = T) ) %>%
      ungroup()
    # SINGLE
    dm <- dm %>%  group_by(sample_id, Interval) %>%
      mutate(median_sqE = median(sq_err),
             mad_sqE    = mad(sq_err, na.rm=T),
             is.out.p    = median_sqE + cut.point * mad_sqE < sq_err) %>%
      ungroup()

    n_out_p <-  nrow(dm %>% filter(is.out.p == T))
    cat(i, " Point outliares: ", n_out_p, "--", n_out_p/size*100 ,  "% \n")

    dm <- dm %>% filter(is.out.p == F)
    if( n_out_p == 0) keep <-F # Stop when no more outliars found

    i <- i+1

  }

  # complete data
  dm_r <- dm_r %>%
    left_join(dm, by = c("Measurement", "Well", "Group", "Time", "plate_id", "Protocol", "Interval"), suffix = c("",".y") ) %>%
    mutate(is.out.p = ifelse(is.out.w == T, NA, replace_na(is.out.p, T)),
           out      = ifelse(is.out.p == F & is.out.w == F, "NO",ifelse(is.out.w == T, "WELL", "SINGLE"))) %>%
    select(-c(contains("y"), "median_sqE", "mad_sqE", "sq_err", "int_mean"))

  # Print summary
  cat("Tolat well outliars: ", nrow(filter(dm_r, is.out.w == T))/nrow(dm_r)*100, "% \n" )
  cat("Tolat single point outliars: ", nrow(filter(dm_r, is.out.p == T))/nrow(dm_r)*100, "% \n" )

  return(dm_r)
}

# -------------------------------------------------------------------------------------- COMPUTE BIOENERGETICS

compute_bioenergetics <- function(dm_r){
  dt_rem <- dm_r %>%
    filter(out == "NO")

  estimates   <- data.frame()
  errors      <- data.frame()
  coef_melted <- data.frame()

  for (sample in unique(dt_rem$sample_id)) {
    d <- dt_rem %>% filter(sample_id == sample)
    fit  <- lm(LOCR ~ -1  + Interval, data = d)
    coef <- summary(fit)$coefficients
    # Trim the row names to extract only factor
    rownames(coef) <- sapply(rownames(coef), function(x){ substr(x,(nchar(x)+1)-4,nchar(x)) })


    # ---------- Interval Estimates
    ints_est <- coef[grepl("Int", rownames(coef)),1 ] # take only rows with "Int"
    #ints <- cbind(Interval = rownames(ints_est), data.frame(ints, row.names=NULL))
    ints_est <- t(ints_est)
    ints_est <- as.data.frame(ints_est) %>% mutate(Sample = sample )
    # ---------- Estimates Standar Error
    ints_sde <- coef[grepl("Int", rownames(coef)),2 ] # take only rows with "Int"
    ints_sde <- t(ints_sde)
    ints_sde <- as.data.frame(ints_sde) %>% mutate(Sample = sample )
    #est_sde <- cbind(est = ints_est, Std_e = ints_sde)

    # have combined data melted for ploting
    est_melted  <- cbind(est = melt(ints_est),melt( ints_sde))
    est_melted  <- est_melted %>% select(c("Sample", Interval = "variable", Estimate = "est.value", SdEr = "value"))
    coef_melted <- rbind(coef_melted, est_melted)


    estimates   <- rbind(estimates, ints_est)
    errors      <- rbind(errors, ints_sde)


  }

  bio_e <- estimates %>%
    mutate(Basal.Resp       = (Int1 - Int4) / Int1,
           ATP.linked.Resp  = (Int1 - Int2) / Int1,
           Proton.Leak      = (Int2 - Int4) / Int2,
           Spare.Resp.Cpcty = Int3 / Int1,
           Maximal.Resp     = Int3 / Int4) %>%
    select(c("Basal.Resp", "ATP.linked.Resp", "Proton.Leak", "Spare.Resp.Cpcty", "Maximal.Resp", "Sample" ))


  return(list(bio_e = bio_e, estim = estimates, err = errors, coef_melted = coef_melted))
}
