#' Resample Daily Dates for Stochastic Weather Generation
#'
#' @description
#' Generates a simulated daily weather sequence using a Markov Chain model and KNN-based resampling,
#' conditioned on annual precipitation and climate statistics. This is designed for use in stochastic
#' weather generators where daily sequences are resampled to match statistical properties of observed data.
#'
#' @param PRCP_FINAL_ANNUAL_SIM Numeric vector. Simulated annual precipitation totals for each simulation year.
#' @param ANNUAL_PRCP Numeric vector. Observed annual precipitation for historical years.
#' @param PRCP Numeric vector. Observed daily precipitation series, length must match `dates.d`.
#' @param TEMP Numeric vector. Observed daily temperature series (e.g., mean temp), length must match `dates.d`.
#' @param TMAX Numeric vector. Observed daily maximum temperature.
#' @param TMIN Numeric vector. Observed daily minimum temperature.
#' @param START_YEAR_SIM Integer. The starting year of the simulation period.
#' @param k1 Integer. Simulation trace index (seed modifier).
#' @param ymax Integer. Number of years to simulate.
#' @param dates.d Data frame. Columns: `date`, `month`, `day`, `wyear` for the observed period.
#' @param sim.dates.d Data frame. Columns: `month`, `day`, `wyear` for the simulation period.
#' @param YEAR_D Integer vector. Water year for each observed day (must match `dates.d`).
#' @param month.start Integer (1-12). Month at which the water year starts.
#' @param knn.annual.sample.num Integer. Number of nearest years to sample in annual KNN (default 50).
#' @param wet.quantile Numeric (0-1). Precip threshold quantile for defining wet days (default 0.2).
#' @param extreme.quantile Numeric (0-1). Precip threshold quantile for extreme wet days (default 0.8).
#' @param dry.spell.change Numeric vector (length 12). Adjustment factors for dry spells per month.
#' @param wet.spell.change Numeric vector (length 12). Adjustment factors for wet spells per month.
#' @param seed Integer or NULL. Random seed for reproducibility (optional).
#'
#' @return
#' Returns a vector of `Date` objects representing the resampled daily dates for the simulation period.
#'
#' @details
#' Uses a Markov Chain to simulate daily precipitation state transitions (dry, wet, extreme)
#' and KNN-based selection for matching temperature/precipitation values.
#'
#' @examples
#' \dontrun{
#' sim_dates <- resampleDates(
#'   PRCP_FINAL_ANNUAL_SIM = runif(3, 900, 1100),
#'   ANNUAL_PRCP = runif(20, 900, 1100),
#'   PRCP = rnorm(7300, 3, 5),
#'   TEMP = rnorm(7300, 18, 5),
#'   TMAX = rnorm(7300, 25, 6),
#'   TMIN = rnorm(7300, 11, 4),
#'   START_YEAR_SIM = 2001,
#'   k1 = 1,
#'   ymax = 3,
#'   dates.d = data.frame(date = seq.Date(as.Date("2000-01-01"), as.Date("2019-12-31"), by = "day"),
#'                        month = rep(1:12, each = 365, length.out = 7300),
#'                        day = rep(1:31, length.out = 7300),
#'                        wyear = rep(2000:2019, each = 365)),
#'   sim.dates.d = data.frame(month = rep(1:12, each = 365, length.out = 1095),
#'                            day = rep(1:31, length.out = 1095),
#'                            wyear = rep(2001:2003, each = 365)),
#'   YEAR_D = rep(2000:2019, each = 365),
#'   month.start = 1,
#'   seed = 123
#' )
#' head(sim_dates)
#' }
#'
#' @export
resampleDates <- function(
  PRCP_FINAL_ANNUAL_SIM,
  ANNUAL_PRCP,
  PRCP,
  TEMP,
  TMAX,
  TMIN,
  START_YEAR_SIM,
  k1,
  ymax,
  dates.d,
  sim.dates.d,
  YEAR_D,
  month.start = 1,
  knn.annual.sample.num = 50,
  wet.quantile = 0.2,
  extreme.quantile = 0.8,
  dry.spell.change = rep(1,12),
  wet.spell.change = rep(1,12),
  seed = NULL

) {

  ######################## --- Input Validation --- ############################

  stopifnot(
    !is.null(PRCP_FINAL_ANNUAL_SIM),
    !is.null(ANNUAL_PRCP),
    !is.null(PRCP),
    !is.null(TEMP),
    !is.null(TMAX),
    !is.null(TMIN),
    !is.null(START_YEAR_SIM),
    !is.null(dates.d),
    !is.null(sim.dates.d),
    !is.null(YEAR_D),
    !is.null(month.start)
  )

  # Numeric checks
  if (!is.numeric(PRCP_FINAL_ANNUAL_SIM) || !is.numeric(ANNUAL_PRCP))
    stop("PRCP_FINAL_ANNUAL_SIM and ANNUAL_PRCP must be numeric vectors.")
  if (!is.numeric(PRCP) || !is.numeric(TEMP) || !is.numeric(TMAX) || !is.numeric(TMIN))
    stop("PRCP, TEMP, TMAX, and TMIN must be numeric vectors.")
  if (!is.numeric(wet.quantile) || wet.quantile < 0 || wet.quantile > 1)
    stop("wet.quantile must be between 0 and 1.")
  if (!is.numeric(extreme.quantile) || extreme.quantile < 0 || extreme.quantile > 1)
    stop("extreme.quantile must be between 0 and 1.")
  if (!is.numeric(knn.annual.sample.num) || knn.annual.sample.num < 1)
    stop("knn.annual.sample.num must be a positive integer.")
  if (!is.numeric(seed) && !is.null(seed))
    stop("seed must be numeric or NULL.")

  # Data frame checks
  if (!is.data.frame(dates.d) || !is.data.frame(sim.dates.d))
    stop("dates.d and sim.dates.d must be data.frames.")

  # Check required columns
  required_cols <- c("date", "month", "day", "wyear")
  missing_cols_dates <- setdiff(required_cols, names(dates.d))
  missing_cols_sim   <- setdiff(c("month", "day", "wyear"), names(sim.dates.d))
  if (length(missing_cols_dates) > 0)
    stop(paste("dates.d is missing columns:", paste(missing_cols_dates, collapse = ", ")))
  if (length(missing_cols_sim) > 0)
    stop(paste("sim.dates.d is missing columns:", paste(missing_cols_sim, collapse = ", ")))

  # Check lengths (for main time series vectors)
  expected_length <- nrow(dates.d)
  if (length(PRCP) != expected_length)
    stop("Length of PRCP does not match dates.d.")
  if (length(TEMP) != expected_length)
    stop("Length of TEMP does not match dates.d.")
  if (length(TMAX) != expected_length)
    stop("Length of TMAX does not match dates.d.")
  if (length(TMIN) != expected_length)
    stop("Length of TMIN does not match dates.d.")
  if (length(YEAR_D) != expected_length)
    stop("Length of YEAR_D does not match dates.d.")

  ######################## --- LOCAL RANDOM SEED --- ##########################

  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({.Random.seed <<- old_seed}, add = TRUE)
    set.seed(seed + k1)
  }

  # Workaround for rlang warning
  month <- day <- wyear <- 0

  # Prepare month list and reference vectors
  month_list <- if(month.start == 1) 1:12 else c(month.start:12, 1:(month.start-1))
  WATER_YEAR_A <- dplyr::filter(dates.d, month == month.start & day == 1) %>% dplyr::pull(wyear)
  DATE_D <- dates.d$date
  MONTH_D <- dates.d$month
  WATER_YEAR_D <- dates.d$wyear
  MONTH_DAY_D <- dates.d[,c("month","day")]
  MONTH_SIM = sim.dates.d$month
  DAY_SIM = sim.dates.d$day
  WATER_YEAR_SIM = sim.dates.d$wyear
  SIM_LENGTH <- length(MONTH_SIM)
  water_year_start = dates.d$wyear[1]

  # Initialize Markov transition arrays and simulated states
	p00_final <- array(NA,SIM_LENGTH)
	p01_final <- array(NA,SIM_LENGTH)
	p02_final <- array(NA,SIM_LENGTH)
	p10_final <- array(NA,SIM_LENGTH)
	p11_final <- array(NA,SIM_LENGTH)
	p12_final <- array(NA,SIM_LENGTH)
	p20_final <- array(NA,SIM_LENGTH)
	p21_final <- array(NA,SIM_LENGTH)
	p22_final <- array(NA,SIM_LENGTH)
	OCCURENCES <- array(0,c(SIM_LENGTH))
	SIM_PRCP <- array(0, c(SIM_LENGTH))
	SIM_TEMP <- array(25, c(SIM_LENGTH))
  SIM_TMAX <- array(30, c(SIM_LENGTH))
	SIM_TMIN <- array(20, c(SIM_LENGTH))
	SIM_DATE <- array(as.Date(paste(water_year_start+1, month.start,"01",sep="-")), SIM_LENGTH)
	count <- 1
	set.seed(seed + k1)
	rn_all <- stats::runif(SIM_LENGTH,0,1)
  kk <- max(round(sqrt(length(ANNUAL_PRCP)),0),round(length(ANNUAL_PRCP),0)*.5)


  # ----- 3. Simulate years -----
	for (y in seq_len(ymax)) {

	  # - For each simulated year, select similar years from the observed record via KNN
	  # - Calculate monthly precipitation thresholds (wet, extreme)
	  # - Calculate Markov transition probabilities
	  # - For each simulated day, use Markov Chain to determine state, then sample a day from historical record

	  # --- 3a. Sample annual years ---
		sim_annual_prcp <- PRCP_FINAL_ANNUAL_SIM[y]
		CUR_YEARS <- knnAnnual(
		  sim_annual_prcp = sim_annual_prcp,
		  ANNUAL_PRCP = ANNUAL_PRCP,
		  WATER_YEAR_A = WATER_YEAR_A,
		  kk = kk,
		  k1 = k1,
		  y = y,
		  seed = seed,
		  y_sample_size = knn.annual.sample.num
		)

		# # Find indices of days in all sampled years in CUR_YEARS
		# conditional_selection <- NULL
		#
		# #Indices of days in the selected 100 years....
		# for (yy in 1:length(CUR_YEARS)) {
		#   conditional_selection <- c(conditional_selection, which(WATER_YEAR_D==CUR_YEARS[yy]))
		# }

		# Find all variables and date indices in the conditional selection
		conditional_selection <- unlist(lapply(CUR_YEARS, function(cur_y) which(WATER_YEAR_D == cur_y)))
		PRCP_CURRENT <- PRCP[conditional_selection]
		TEMP_CURRENT <- TEMP[conditional_selection]
		TMAX_CURRENT <- TMAX[conditional_selection]
		TMIN_CURRENT <- TMIN[conditional_selection]
		DATE_D_CURRENT <- DATE_D[conditional_selection]
		MONTH_D_CURRENT <- MONTH_D[conditional_selection]
		YEAR_D_CURRENT <- YEAR_D[conditional_selection]
		MONTH_DAY_D_CURRENT <- MONTH_DAY_D[conditional_selection,]

		# --- 3b. Thresholds ---
		wet_threshold <- rep(stats::quantile(PRCP_CURRENT, wet.quantile, names = FALSE), 12)
    extreme_threshold <- sapply(1:12, function(m)
      stats::quantile(PRCP_CURRENT[which(MONTH_D_CURRENT==month_list[m] & PRCP_CURRENT > wet_threshold[m])],
        extreme.quantile, names = F))

    ############################################################################

		#Define lagged variables on daily time-series (for current year set)
		PRCP_LAG0  <- PRCP_CURRENT[2:length(PRCP_CURRENT)]
		PRCP_LAG1  <- PRCP_CURRENT[1:(length(PRCP_CURRENT)-1)]
		MONTH_LAG0 <- MONTH_D_CURRENT[2:length(PRCP_CURRENT)]
		MONTH_LAG1 <- MONTH_D_CURRENT[1:(length(PRCP_CURRENT)-1)]

		# --- 3d. Vectorized Markov transition probability calculation ---
		probs <- calculateMarkovProbs(
		  PRCP_LAG0, PRCP_LAG1, MONTH_LAG0, MONTH_LAG1,
		  wet_threshold, extreme_threshold, month_list,
		  MONTH_SIM, WATER_YEAR_SIM, y, START_YEAR_SIM,
		  dry.spell.change, wet.spell.change, SIM_LENGTH
		)

		p00_final <- probs$p00_final; p01_final <- probs$p01_final; p02_final <- probs$p02_final
		p10_final <- probs$p10_final; p11_final <- probs$p11_final; p12_final <- probs$p12_final
		p20_final <- probs$p20_final; p21_final <- probs$p21_final; p22_final <- probs$p22_final

		############################################################################
		############################################################################

		# MARKOV-CHAIN AND DAILY KNN SAMPLING.......................................
		for (j in 1:365) {


		  count <- count + 1

			if (count <= SIM_LENGTH) {

  			  # Random number generated for the current day
    			rn <- rn_all[(count-1)]

    			# If day is in state 0
    			if (OCCURENCES[(count-1)]==0) {
    				pp1 <- p00_final[(count-1)]
    				pp2 <- p00_final[(count-1)] + p01_final[(count-1)]
    			}

    			# If day is in state 1
    			if (OCCURENCES[(count-1)]==1) {
    				pp1 <- p10_final[(count-1)]
    				pp2 <- p10_final[(count-1)] + p11_final[(count-1)]
    			}

    			# If day is in state 2
    			if (OCCURENCES[(count-1)]==2) {
    				pp1 <- p20_final[(count-1)]
    				pp2 <- p20_final[(count-1)] + p21_final[(count-1)]
    			}

    			# Set state for the next day
    			if(is.na(pp1)) pp1 <- 0
    			if(is.na(pp2)) pp2 <- 0

    			if(rn < pp1) {
    				OCCURENCES[count] <- 0
    			} else if (rn >= pp1 & rn < pp2) {
    				OCCURENCES[count] <- 1
    			} else {
    				OCCURENCES[count] <- 2
    			}

    			# Current day's month of the year and day of the year indices
    			m <- MONTH_SIM[(count-1)]
    			mmm <- which(month_list==m)
    			d <- DAY_SIM[(count-1)]

    			# Current occurrance
    			cur_OCCERENCE <- OCCURENCES[(count-1)]

    			# Next day's occurrance
    			next_OCCURENCE <- OCCURENCES[(count)]

    			# Subset non-zero days?
    			cur_day <- which(MONTH_DAY_D_CURRENT[,1]==m & MONTH_DAY_D_CURRENT[,2]==d)
    			cur_day <- c((cur_day-3),(cur_day-2),(cur_day-1),cur_day,(cur_day+1),(cur_day+2),(cur_day+3))
    			cur_day <- subset(cur_day,cur_day > 0)

    			# Set current day's state.............
    			if (cur_OCCERENCE==0 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    			if (cur_OCCERENCE==0 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    			if (cur_OCCERENCE==0 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}
    			if (cur_OCCERENCE==1 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    			if (cur_OCCERENCE==1 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    			if (cur_OCCERENCE==1 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}
    			if (cur_OCCERENCE==2 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    			if (cur_OCCERENCE==2 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    			if (cur_OCCERENCE==2 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}

    			# If length of the state index is zero......
    			if (length(cur_day_cur_state)==0) {
    				cur_day <- which(MONTH_DAY_D_CURRENT[,1]==m & MONTH_DAY_D_CURRENT[,2]==d)
    				cur_day_final <- array(NA,length(cur_day)*61)
    				cur_day_window <- seq(-30,30)
    				for (cc in 1:61) {
    					cur_day_final[(1 + length(cur_day)*(cc-1)):(length(cur_day) + length(cur_day)*(cc-1))] <- (cur_day+cur_day_window[cc])
    				}
    				cur_day <- subset(cur_day_final, cur_day_final > 0)
    				if (cur_OCCERENCE==0 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    				if (cur_OCCERENCE==0 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    				if (cur_OCCERENCE==0 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}
    				if (cur_OCCERENCE==1 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    				if (cur_OCCERENCE==1 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    				if (cur_OCCERENCE==1 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>wet_threshold[m] & PRCP_CURRENT[cur_day]<=extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}
    				if (cur_OCCERENCE==2 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]<=wet_threshold[m])}
    				if (cur_OCCERENCE==2 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>wet_threshold[m] & PRCP_CURRENT[(cur_day+1)]<=extreme_threshold[mmm])}
    				if (cur_OCCERENCE==2 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>extreme_threshold[mmm] & PRCP_CURRENT[(cur_day+1)]>extreme_threshold[mmm])}
    			}

  			#	Define Possible days to choose from and set their values
  			possible_days <- cur_day[cur_day_cur_state]
  			PRCP_TODAY <- PRCP_CURRENT[possible_days]
  			TEMP_TODAY <- TEMP_CURRENT[possible_days]
  			PRCP_TOMORROW <- PRCP_CURRENT[possible_days+1]
				TEMP_TOMORROW <- TEMP_CURRENT[possible_days+1]
				TMAX_TOMORROW <- TMAX_CURRENT[possible_days+1]
				TMIN_TOMORROW <- TMIN_CURRENT[possible_days+1]
  			DATE_TOMORROW <- DATE_D_CURRENT[possible_days+1]

  			cur_sim_PRCP <- SIM_PRCP[(count-1)]
  			cur_sim_TEMP <- SIM_TEMP[(count-1)]

  			mm <- which(MONTH_D_CURRENT==m)
  			mm_p <- which(MONTH_D_CURRENT==m & PRCP_CURRENT > 0)

  			sd_monthly_TEMP <- stats::sd(TEMP_CURRENT[mm], na.rm=TRUE)
  			sd_monthly_PRCP <- stats::sd(PRCP_CURRENT[mm_p], na.rm=TRUE)

  			mean_monthly_TEMP <- mean(TEMP_CURRENT[mm], na.rm=TRUE)
  			mean_monthly_PRCP <- mean(PRCP_CURRENT[mm_p], na.rm=TRUE)

  			k <- round(sqrt(length(possible_days)))

  			RESULT <- knnDaily(
  			  cur_sim_PRCP = cur_sim_PRCP,
  			  cur_sim_TEMP = cur_sim_TEMP,
  			  PRCP_TODAY = PRCP_TODAY,
  			  TEMP_TODAY = TEMP_TODAY,
  			  k = k,
  			  sd_monthly_PRCP = sd_monthly_PRCP,
  			  sd_monthly_TEMP = sd_monthly_TEMP,
  			  mean_monthly_PRCP = mean_monthly_PRCP,
  			  mean_monthly_TEMP = mean_monthly_TEMP,
  			  k1 = k1,
  			  count = count,
  			  seed = seed
  			)

  			# Fallback if no valid RESULT
  			if (is.na(RESULT) || length(RESULT) == 0 ||
  			    RESULT < 1 || RESULT > length(PRCP_TOMORROW)) {
  			  # Fallback: sample randomly, or assign NA, or use 1st available
  			  if (length(PRCP_TOMORROW) > 0) {
  			    fallback_index <- sample(seq_along(PRCP_TOMORROW), 1)
  			  } else {
  			    fallback_index <- NA_integer_
  			  }
  			  RESULT <- fallback_index
  			}

  			# Now safely assign, even if fallback_index is NA
  			SIM_PRCP[count] <- if (!is.na(RESULT)) PRCP_TOMORROW[RESULT] else NA_real_
  			SIM_TEMP[count] <- if (!is.na(RESULT)) TEMP_TOMORROW[RESULT] else NA_real_
  			SIM_TMAX[count] <- if (!is.na(RESULT)) TMAX_TOMORROW[RESULT] else NA_real_
  			SIM_TMIN[count] <- if (!is.na(RESULT)) TMIN_TOMORROW[RESULT] else NA_real_
  			SIM_DATE[count] <- if (!is.na(RESULT)) DATE_TOMORROW[RESULT] else as.Date(NA)

      	#SIM_PRCP[count] <- PRCP_TOMORROW[RESULT]
    		#SIM_TEMP[count] <- TEMP_TOMORROW[RESULT]
    		#SIM_TMAX[count] <- TMAX_TOMORROW[RESULT]
    		#SIM_TMIN[count] <- TMIN_TOMORROW[RESULT]
    		#SIM_DATE[count] <- DATE_D_CURRENT[which(as.numeric(DATE_D_CURRENT)==DATE_TOMORROW[RESULT])][1]

			} # if-condition count close

		} # daily-counter close

	} # year-counter close

  class(SIM_DATE) <- "Date"

	return(SIM_DATE)

}

