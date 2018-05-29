##' Convert NONMEM/Monolix data to RxODE/nlmixr style data
##'
##' @param nonmem.data NONMEM dataset
##' @return RxODE compatible data.
##' @author Yuan Xiong and Matthew L. Fidler
##' @export
nmDataConvert <- function(nonmem.data)
{
    if (!is(nonmem.data, "data.frame")){
        stop("Input must be a data frame.")
    }
    d <- as.data.frame(nonmem.data)
    col.names <- colnames(d)
    col.names <- toupper(col.names)
    row.names(d) <- NULL;  ## Make sure row.names are named by R's internals
    ## Fix a few possible problems
    colnames(d) <- col.names
    diff <- setdiff(c("TIME","DV"), col.names);
    if (length(diff) > 0){
        stop(sprintf("The conversion is missing the following required column(s): %s", paste(diff, collapse=", ")));
    }

    if (!any(col.names == "EVID") && any(col.names == "MDV")){
        d$EVID <- d$MDV;
        d <- d[, names(d) != "MDV"];
    }

    if (any(d$EVID > 4) | any(d$EVID < 0)){
        stop("EVIDs do not match NONMEM-style EVIDs.");
    }

    if (any(d$EVID >= 3)){
        stop("A reset of the ODE system is not currently supported.")
    }

    if (any(col.names == "SS" & any(d$SS > 0))){
        stop("Steady State dosing events are not currently supported.")
    }

    ## Handle DATE TIME; DAT1 TIME; DAT2 TIME and DAT3 TIME
    ## Note NONMEM handles dates of the format DAY-MONTH and DAY as well for the DATE class of objects.
    ## It is too complex to handle, and not very common so it will throw an error
    do.date <- FALSE;
    bad.msg <- "Dates formated as MONTH-DAY or DAY alone are not supported in this conversion."
    dup.date <- "Dates can only be specified by one of DATE, DAT1, DAT2, DAT3 / TIME;  This dataset has multiple DATE columns."
    check.bad <- function(d){
        d <- paste(d);
        if (any(unlist(lapply(strsplit(d, "[^0-9]+"), length)) != 3)){
            stop(bad.msg)
        }
        return(d)
    }
    if (any(col.names == "DATE")){
        ##  Month Day Year
        dat.reg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- check.bad(d$DATE)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%m-%d-%y %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%m-%d-%Y %H:%M")
        d <- d[, -which(names(d) == "DATE")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT1")){
        if (do.date){
            stop(dup.date)
        }
        ## DAT1   day month year
        dat.reg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- check.bad(d$DAT1)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%d-%m-%y %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%d-%m-%Y %H:%M")
        d <- d[, -which(names(d) == "DAT1")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT2")){
        ## DAT2   year month day
        if (do.date){
            stop(dup.date)
        }
        dat.reg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- check.bad(d$DAT2)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%y-%m-%d %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%Y-%m-%d %H:%M")
        d <- d[, -which(names(d) == "DAT2")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT3")){
        ## DAT3   year day month
        if (do.date){
            stop(dup.date)
        }
        dat.reg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- check.bad(d$DAT3)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%y-%d-%m %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%Y-%d-%m %H:%M")
        d <- d[, -which(names(d) == "DAT3")];
        do.date <- TRUE;
    }
    if (do.date){
        if (any(is.na(d$DATE.TIME))){
            stop("The date time format was not correctly specified.")
        }
    }
    if (any(col.names == "SS") && any(d$SS == 2)){
        stop("Steady state = 2 not supported.");
    }

    if (any(col.names == "RATE") && any(d$RATE < 0)){
        stop("RxODE does not currently estimating infusion rates or durations. (Found RATE < 0)");
    }

    if (!any(col.names == "EVID") && !any(col.names == "MDV")){
        d$EVID <- 0;
        warning("Assumed all DV values are observations. (EVID=0)")
    }
    if (any(d$EVID == 2)){
        warning("EVID=2 is dropped from RxODE/nlmixr datasets.")
        d <- d[d$EVID != 2, ];
    }

    if (!any(col.names == "AMT")){
        warning("Assuming AMT=0.")
        d$AMT <- 0;
    }
    if (!any(col.names == "ID")){
        warning("ID=1 added to dataset");
        d$ID <- 1;
    }

    drop.cmt <- FALSE
    if (any(d$EVID > 0) && !any(col.names == "CMT")){
        drop.cmt <- TRUE
        d$CMT <- 1;
        d$CMT[d$EVID == 0] <- NA;
    }

    if (do.date){
        ## Sort by date/time (though this should have been done already...)
        d <- d[order(d$ID, d$DATE.TIME, -d$EVID), ];
        d$TIME <- as.vector(unlist(sapply(unique(d$ID), function(id){
            d0 <- d[d$ID == id, ];
            return(as.numeric(difftime(d0$DATE.TIME, d0$DATE.TIME[1], units="hours")))
        })))
        d <- d[, -which(names(d) == "DATE.TIME")];
    }

    ## TODO list:
    ## messages to user (ambiguous names; unique ID; etc.)

    ## handle II/ADDL
    if (any(col.names == "II") & any(col.names == "ADDL")) {
        ds.addl <- d[d$II>0 & d$ADDL>0, ]
        ds.exp <- ds.addl[rep(row.names(ds.addl), ds.addl$ADDL), 1:ncol(ds.addl)]
        ## expand dosing records
        ds.exp$TIME <- unlist(sapply(1:nrow(ds.addl),
                                     function(idx) {
            seq(ds.addl$ADDL[idx])*ds.addl$II[idx]+ds.addl$TIME[idx]
        }))
        ## attach expanded dosing records to the original dataset
        d <- rbind(d, ds.exp)

        ## sort data and remove II/ADDL to avoid confusion
        d <- d[order(d$ID, d$TIME, -d$EVID), -which(names(d) %in% c("II","ADDL"))]
    }

    ## handle EVID
    if (any(col.names == "RATE")) { ## have infusion records
        ## actually bolus
        d.bol <- d[is.na(d$RATE) | d$RATE == 0, ]
        d.bol$EVID <- ifelse(d.bol$EVID>0, d.bol$CMT*100+d.bol$EVID, d.bol$EVID)
        ## another case: EVID with no CMT => assume CMT=1

        ## infusion
        d.inf <- d[d$RATE>0, ]
        d.inf$EVID <- .Call(`_nlmixr_convertEvidRate`, as.integer(d.inf$EVID), as.integer(d.inf$CMT));
        ## add end of infusion records
        d.inf.end <- d.inf[d.inf$AMT>0, ]
        d.inf.end$TIME <- d.inf.end$TIME + d.inf.end$AMT/d.inf.end$RATE
        d.inf.end$RATE <- -d.inf.end$RATE
        ## attach end of infusion records
        d.inf.new <- rbind(d.inf, d.inf.end)
        d.inf.new$AMT <- d.inf.new$RATE

        ## bind all records together
        dat <- rbind(d.bol, d.inf.new)
        dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]

    } else { ## no infusion: EVID is the only thing to change
        ## last two digits: EVID=01
        ## mid two digits: label for CMT in original NONMEM dataset
        ## upper two digits: infusion indicator, set to 0
        d$EVID <- .Call(`_nlmixr_convertEvid`, as.integer(d$EVID), as.integer(d$CMT));
        dat <- d
    }
    dat
}

