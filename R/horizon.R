#' Create .\emph{skyscapeR.horizon} object from Az/Alt data
#'
#' This function creates a \emph{skyscapeR.horizon} object from measurements of
#' azimuth and altitude.
#' @param az Array of azimuth values
#' @param alt Array of altitude values.
#' @param loc Location, a vector containing the latitude and longitude of
#' the location, in this order.
#' @param name Name of site.
#' @seealso \code{\link{plotHor}}
#' @export
#' @examples
#' # Create a skyscapeR.horizon from 5 measurements:
#' az <- c(0,90,180,270,360)
#' alt <- c(0,5,5,0,0)
#' hor <- createHor(az, alt, c(40.1,-8), 'Test')
#' plotHor(hor)
createHor = function(az, alt, loc, name) {
  # return result
  hor <- c()
  hor$alt <- alt
  hor$az <- az
  hor$georef <- loc; names(hor$georef) <- c('Lat','Lon')
  hor$name <- name
  class(hor) <- "skyscapeR.horizon"
  return(hor)
}



#' Exports a \emph{skyscapeR.horizon} object into \emph{Stellarium} format
#'
#' This function exports any \emph{skyscapeR.horizon} object into the landscape
#' format of \emph{Stellarium}, ready to be imported.
#' @param hor Horizon data in \emph{skyscapeR.horizon} format.
#' @param name Horizon name to be displayed in \emph{Stellarium}, if different
#' from one in \emph{skyscapeR.horizon} object.
#' @param author (Optional) Author, to be included in \emph{landscape.ini} file.
#' @param description (Optional) Description, to be included in \emph{landscape.ini} file.
#' @param ground_col Colour of ground. Defaults to \emph{Stellarium}'s default.
#' @param hor_col Colour of horizon line. Defaults to \emph{Stellarium}'s default.
#' @seealso \code{\link{createHor}}, \code{\link{download.HWT}}, \code{\link{plotHor}}
#' @references \href{http://www.stellarium.org/}{Stellarium: a free open source planetarium}
#' @export
#' @import utils
#' @examples
#' # Downloads horizon data from HeyWhatsThat and exports it into Stellarium:
#' hor <- download.HWT('HIFVTBGK')
#' exportHor(hor, name='Test', description='Test horizon export to Stellarium')
exportHor = function(hor, name, author="skyscapeR", description, ground_col, hor_col) {
  if (class(hor) != 'skyscapeR.horizon') { stop('No skyscapeR.hor object found.') }

  if (missing(name)) { name = hor$name }
  if (missing(description)) { description <- paste0("Horizon created using skyscapeR ", packageVersion('skyscapeR'), ".") }
  if (substr(description,nchar(description),nchar(description)) != ".") { description <- paste0(description,". Horizon created using skyscapeR ", packageVersion('skyscapeR'), ".") }
  if (missing(ground_col)) {ground_col = ".15,.45,.15"}
  if (missing(hor_col)) {hor_col = ".75,.45,.45"}

  # Horizon data
  data <- data.frame(x = hor$az, y = hor$alt)
  write.table(data, file="horizon.txt", sep = " ", row.names=F, col.names=F)

  # Landscape.ini file
  fileConn<-file("landscape.ini")
  string.text <- c("[landscape]", paste0("name = ",name), "type = polygonal", paste0("author = ",author),
                   paste0("description = ", description), "polygonal_horizon_list = horizon.txt",
                   "polygonal_angle_rotatez = 0", paste0("ground_color = ",ground_col),
                   paste0("horizon_line_color =  ",hor_col), "",
                   "[location]", "planet = Earth", paste0("latitude = ",hor$Lat,"d"),
                   paste0("longitude = ",hor$Lon,"d"), paste0("altitude = ",round(hor$Elev*0.3,0)))
  writeLines(string.text, fileConn)
  close(fileConn)

  # Save as zip file
  zip(zipfile=paste0(name,'-horizon.zip'), files=c("landscape.ini", "horizon.txt"))
  file.remove(c('landscape.ini', 'horizon.txt'))
}



#' Download horizon data from \emph{HeyWhatsThat}
#'
#' This function downloads horizon data from \emph{HeyWhatsThat},
#' given its ID, and saves it as a \emph{skyscapeR.horizon} object.
#' @param HWTID This is the 8 character ID attributed by
#' \emph{HeyWhatsThat.com}
#' @export
#' @import utils
#' @references \href{http://heywhatsthat.com/}{HeyWhatsThat.com}
#' @examples
#' # Retrieve horizon data for \href{https://www.heywhatsthat.com/?view=HIFVTBGK}{Liverpool Cathedral}:
#' hor <- download.HWT('HIFVTBGK')
download.HWT = function(HWTID) {
  if (nchar(HWTID) != 8) { stop('Incorrect HeyWhatsThat ID.') }

  ## Horizon metadata
  test <- readLines(paste0("http://www.heywhatsthat.com/iphone/pan.cgi?id=",HWTID))
  hor <- c()

  # Lat/Lon/Elev
  mypattern = '<div class=\"details_data\">([^<]*)</div>'
  datalines = grep(mypattern,test,value=TRUE)
  getexpr = function(s,g)substring(s,g,g+attr(g,'match.length')-1)
  gg = gregexpr(mypattern,datalines)
  matches = mapply(getexpr,datalines,gg)
  result = gsub(mypattern,'\\1',matches)
  names(result) = NULL
  result[c(1,2,4)]

  hor$Lat <- as.numeric(strtrim(result[1],regexpr("&deg", result[1])[1]-1))
  if (substr(result[1],nchar(result[1]),nchar(result[2])) == "S") {hor$Lat <- -hor$Lat}
  hor$Lon <- as.numeric(strtrim(result[2],regexpr("&deg", result[2])[1]-1))
  if (substr(result[2],nchar(result[2]),nchar(result[2])) == "W") {hor$Lon <- -hor$Lon}
  aux <- strtrim(result[4],regexpr("&nbsp;", result[4])[1]-1)
  hor$Elev <- as.numeric(substr(aux,2,nchar(aux)))

  # Site Name
  grep("pan_top_title", test)
  mypattern = '<div id=\"pan_top_title\" class=\"ellipsis\" style=\"position: absolute; top: 46px; width: 296px\">([^<]*)</div>'
  datalines = grep(mypattern,test,value=TRUE)
  getexpr = function(s,g)substring(s,g,g+attr(g,'match.length')-1)
  gg = gregexpr(mypattern,datalines)
  matches = mapply(getexpr,datalines,gg)
  result = gsub(mypattern,'\\1',matches)
  names(result) = NULL
  hor$Name <- result

  ## Horizon data
  hor.ex <- unique(substr(list.files(tempdir()),1,8)) # check if already downloaded


  if (sum(HWTID == hor.ex) == 0) {
    if (NROW(hor.ex) > 500) {
      # delete oldest
      details <- file.info(file.path(tempdir(),list.files(file.path(tempdir()))))
      details <- details[with(details, order(as.POSIXct(ctime))), ]
      files <- unique(substr(rownames(details),10,17))[1]
      files <- paste0(file.path(tempdir()),"/",c(files, paste0(files,"-0-1")), ".png")
      file.remove(files)
    }

    # download new one
    curdir <- getwd()
    setwd(tempdir())
    download.file(paste0('http://www.heywhatsthat.com/results/',HWTID,'/image_north.png'), mode ='wb', destfile=paste0(HWTID,'.png'), quiet = T)
    download.file(paste0('http://www.heywhatsthat.com/bin/image-north-with-alts.png?alts=0%201&id=',HWTID), mode ='wb', destfile=paste0(HWTID,'-0-1.png'), quiet = T)
    setwd(curdir)
  }

  horizon <- png::readPNG(file.path(tempdir(), paste0(HWTID,'.png')))
  horizon.alt <- png::readPNG(file.path(tempdir(), paste0(HWTID,'-0-1.png')))

  # clean-up : removal of triangles and cardinals
  ind <- which(round(horizon[,,1]*256) == 199, arr.ind = T)
  if (NROW(ind)>0) {
    for (i in 1:NROW(ind)) {
      horizon[ind[i,1],ind[i,2],1:3] <- c(251,251,251)/256
    }
  }
  ind <- which(round(horizon[,,1]*256) == 0 & round(horizon[,,2]*256) == 0 & round(horizon[,,3]*256) == 0, arr.ind = T)
  if (NROW(ind)>0) {
    for (i in 1:NROW(ind)) {
      horizon[ind[i,1],ind[i,2],1:3] <- c(251,251,251)/256
    }
  }

  # get altitude values
  ind <- which(round(horizon.alt[,,1]*256)==140, arr.ind=T);
  lines <- ind[1:2,1]
  altitude <- array(NA, dim=dim(horizon)[1:2]);
  if (diff(lines)==0) { lines[2] <- dim(altitude)[1] }
  altitude[lines[2],] = 0; altitude[lines[1],] = 1;

  # linear extrapolation of altitude
  m <- 1./(lines[1]-lines[2]);
  b <- 1 - m*lines[1];
  for (k in 1:dim(altitude)[1]) {
    altitude[k,] <- m*k + b;
  }

  # identification of altitude
  alt <- array(NA,c(1,dim(altitude)[2]))
  for (l in 1:dim(altitude)[2]) {
    for (k in 1:dim(altitude)[1]) {
      if (round(horizon[k,l,1]*256) != 251) {
        alt[l] <- altitude[k,l]
        break
      }
    }
  }

  # return result
  hor$alt <- as.vector(alt)
  hor$az <- seq(0,359,length.out = NCOL(horizon))-180
  ind <- which(hor$az<0)
  hor$az <- c(hor$az,hor$az[ind]+360); hor$az <- hor$az[-ind]
  hor$alt <- c(hor$alt,hor$alt[ind]); hor$alt <- hor$alt[-ind]
  hor$georef <- c(hor$Lat,hor$Lon); names(hor$georef) <- c('Lat','Lon')
  hor$ID <- HWTID
  class(hor) <- "skyscapeR.horizon"
  return(hor)
}
