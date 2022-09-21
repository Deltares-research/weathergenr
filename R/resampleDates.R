
#' Dissaggregation function
#'
#' @param k1 placeholder
#' @param ymax placeholder
#' @param PRCP_FINAL_ANNUAL_SIM placeholder
#' @param ANNUAL_PRCP placeholder
#' @param PRCP placeholder
#' @param TEMP placeholder
#' @param YEAR_D placeholder
#' @param START_YEAR_SIM placeholder
#' @param dates.d placeholder
#' @param sim.dates.d placeholder
#' @param month.start placeholder
#' @param wet.quantile Precipitation threshold (mm) to distinguish between wet and dry days
#' @param extreme.quantile quantile threshold to distinguish between extreme wet days
#' @param knn.annual.sample.num placeholder
#' @param seed random seed value
#'
#' @return
#' @export
resampleDates <- function(
  PRCP_FINAL_ANNUAL_SIM = NULL,
  ANNUAL_PRCP = NULL,
  PRCP = NULL,
  TEMP = NULL,
  START_YEAR_SIM = NULL,
  k1 = NULL,
  ymax = NULL,
  dates.d = NULL,
  sim.dates.d = NULL,
  YEAR_D = NULL,
  month.start = NULL,
  knn.annual.sample.num = 50,
  wet.quantile = 0.2,
  extreme.quantile = 0.8,
  seed = NULL)

  {

  # If ne seed provided, generate a random number
  if(is.null(seed)) seed = sample.int(1e10,1)

  # Workaround for rlang warning
  month <- day <- wyear <- 0

  if(month.start == 1) {
    month_list <- 1:12
  } else {
    month_list <- c(month.start:12,1:(month.start-1))
  }

  WATER_YEAR_A <- dates.d %>%
    dplyr::filter(month == month.start & day == 1) %>% pull(wyear)

  DATE_D <- dates.d$date
  MONTH_D <- dates.d$month
  WATER_YEAR_D <- dates.d$wyear
  MONTH_DAY_D <- dates.d[,c("month","day")]

  MONTH_SIM = sim.dates.d$month
  DAY_SIM = sim.dates.d$day

  WATER_YEAR_SIM = sim.dates.d$wyear
  SIM_LENGTH <- length(MONTH_SIM)

  #***
  water_year_start = dates.d$wyear[1]
  water_year_end = dates.d$wyear[nrow(dates.d)]
  #****

	# Vectors to store transition probabilities
	p00_final <- array(NA,SIM_LENGTH)
	p01_final <- array(NA,SIM_LENGTH)
	p02_final <- array(NA,SIM_LENGTH)
	p10_final <- array(NA,SIM_LENGTH)
	p11_final <- array(NA,SIM_LENGTH)
	p12_final <- array(NA,SIM_LENGTH)
	p20_final <- array(NA,SIM_LENGTH)
	p21_final <- array(NA,SIM_LENGTH)
	p22_final <- array(NA,SIM_LENGTH)

	# Vetor to store simulated results???
	OCCURENCES <- array(0,c(SIM_LENGTH))
	SIM_PRCP <- array(0, c(SIM_LENGTH))
	SIM_TEMP <- array(25, c(SIM_LENGTH))

	SIM_DATE <- array(as.Date(paste(water_year_start+1, month.start,"01",sep="-")), SIM_LENGTH)

	# Current Stochastic trace....
	count <- 1

	# Generate random numbers btw 0 and 1 for each day of the simulation period
	set.seed(seed+k1)
	rn_all <- stats::runif(SIM_LENGTH,0,1)

  # knn sample size
  kk <- max(round(sqrt(length(ANNUAL_PRCP)),0),round(length(ANNUAL_PRCP),0)*.5)

	# For each year start sampling....
	for (y in 1:ymax) {

	  #print(y)

	  # Current simulated annual precip at y
		sim_annual_prcp <- PRCP_FINAL_ANNUAL_SIM[y]

		# Sample 100 similar years from the historical record for the current year
		CUR_YEARS <- knnAnnual(
		  sim_annual_prcp = sim_annual_prcp,
		  ANNUAL_PRCP = ANNUAL_PRCP,
		  WATER_YEAR_A = WATER_YEAR_A,
		  kk = kk,
		  k1 = k1,
		  y = y,
		  seed = seed,
		  y_sample_size = knn.annual.sample.num)

		# Find indices of days in all sampled years in CUR_YEARS
		conditional_selection <- NULL

		#Indices of days in the selected 100 years....
		for (yy in 1:length(CUR_YEARS)) {
		  conditional_selection <- c(conditional_selection, which(WATER_YEAR_D==CUR_YEARS[yy]))
		}

		# Find all variables and date indices in the current selection years
		PRCP_CURRENT <- PRCP[conditional_selection]
		TEMP_CURRENT <- TEMP[conditional_selection]
		DATE_D_CURRENT <- DATE_D[conditional_selection]
		MONTH_D_CURRENT <- MONTH_D[conditional_selection]
		YEAR_D_CURRENT <- YEAR_D[conditional_selection]
		MONTH_DAY_D_CURRENT <- MONTH_DAY_D[conditional_selection,]

		# Calculate thresholds for each month
    wet_threshold <- sapply(1:12, function(m)
      stats::quantile(PRCP_CURRENT[which(MONTH_D_CURRENT==month_list[m])],
        wet.quantile, names = F))

    extreme_threshold <- sapply(1:12, function(m)
      stats::quantile(PRCP_CURRENT[which(MONTH_D_CURRENT==month_list[m] & PRCP_CURRENT > wet_threshold[m])],
        extreme.quantile, names = F))

		#Define lagged variables on daily time-series (for current year set)
		PRCP_LAG0  <- PRCP_CURRENT[2:length(PRCP_CURRENT)]
		PRCP_LAG1  <- PRCP_CURRENT[1:(length(PRCP_CURRENT)-1)]
		MONTH_LAG0 <- MONTH_D_CURRENT[2:length(PRCP_CURRENT)]
		MONTH_LAG1 <- MONTH_D_CURRENT[1:(length(PRCP_CURRENT)-1)]
		YEAR_LAG0  <- YEAR_D_CURRENT[2:length(PRCP_CURRENT)]
		YEAR_LAG1  <- YEAR_D_CURRENT[1:(length(PRCP_CURRENT)-1)]

		# Calculate monthly trahsition probabilities
		for (m in 1:12) {

		  # day of the year index in each month
		  x <- which(MONTH_LAG1==month_list[m])

			r <- which(MONTH_SIM==month_list[m] & WATER_YEAR_SIM == (y+START_YEAR_SIM))

			CUR_PRCP0 <- PRCP_LAG0[x]
			CUR_PRCP1 <- PRCP_LAG1[x]

			# Transition probabilities
			p00_final[r] <- length(which(PRCP_LAG1[x]<=wet_threshold[m] & PRCP_LAG0[x]<=wet_threshold[m])) / length(which(PRCP_LAG1[x]<=wet_threshold[m]))
			p01_final[r] <- length(which(PRCP_LAG1[x]<=wet_threshold[m] & PRCP_LAG0[x]>wet_threshold[m] & PRCP_LAG0[x]<=extreme_threshold[m])) / length(which(PRCP_LAG1[x]<=wet_threshold[m]))
			p02_final[r] <- length(which(PRCP_LAG1[x]<=wet_threshold[m] & PRCP_LAG0[x]>extreme_threshold[m])) / length(which(PRCP_LAG1[x]<=wet_threshold[m]))
			p10_final[r] <- length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m] & PRCP_LAG0[x]<=wet_threshold[m])) / length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m]))
			p11_final[r] <- length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m] & PRCP_LAG0[x]>wet_threshold[m] & PRCP_LAG0[x]<=extreme_threshold[m])) / length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m]))
			p12_final[r] <- length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m] & PRCP_LAG0[x]>extreme_threshold[m])) / length(which(PRCP_LAG1[x]>wet_threshold[m] & PRCP_LAG1[x]<=extreme_threshold[m]))
			p20_final[r] <- length(which(PRCP_LAG1[x]>extreme_threshold[m] & PRCP_LAG0[x]<=wet_threshold[m])) / length(which(PRCP_LAG1[x]>extreme_threshold[m]))
			p21_final[r] <- length(which(PRCP_LAG1[x]>extreme_threshold[m] & PRCP_LAG0[x]>wet_threshold[m] & PRCP_LAG0[x]<=extreme_threshold[m])) / length(which(PRCP_LAG1[x]>extreme_threshold[m]))
			p22_final[r] <- length(which(PRCP_LAG1[x]>extreme_threshold[m] & PRCP_LAG0[x]>extreme_threshold[m])) / length(which(PRCP_LAG1[x]>extreme_threshold[m]))

		} #month-counter close

		############################################################################
		############################################################################

		# MARKOV-CHAIN AND DAILY KNN SAMPLING.......................................
		for (j in 1:365) {

		  #print(j)

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
  			  DATE_TOMORROW = DATE_TOMORROW,
  			  k = k,
  			  sd_monthly_PRCP = sd_monthly_PRCP,
  			  sd_monthly_TEMP = sd_monthly_TEMP,
  			  mean_monthly_PRCP = mean_monthly_PRCP,
  			  mean_monthly_TEMP = mean_monthly_TEMP,
  			  k1 = k1,
  			  count = count,
  			  seed = seed)

  			SIM_DATE[count] <- DATE_D_CURRENT[which(as.numeric(DATE_D_CURRENT)==RESULT)][1]

			} # if-condition count close

		} # daily-counter close

	} # year-counter close

  class(SIM_DATE) <- "Date"

	return(SIM_DATE)

}

