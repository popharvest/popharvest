#' Sustainability harvest index visualizations from PEG and PTL functions
#'
#' Create visualizations of sustainability harvest indices calculated from \code{PEG} or \code{PTL} functions
#'
#' @param data output dataframe from \code{PEG} and \code{PTL} functions. Other inputs will give an error.
#' @param Sp number ; the identification number assigned to each studied species
#' @param IDSp boolean ; default \code{TRUE}. If \code{TRUE}, print the identification number assigned to each studied species. If \code{FALSE}, doesn't print the identification number assigned to the studied species.
#' @param NameSp character string ; the name of the studies species. Add the names of the studies species to the main title of the plot. Example : c("Species A", "Species B")
#' @param title character string ; default "Distribution of sustainability harvest index" ; main title
#' @param title.x character string ; default "Sustainability harvest index" ; title x label
#' @param title.y character string ; default "Simulations number" ; title y label
#' @param color1 character string ; default green (#008000)
#' @param color2 character string ; default red (#B40000)
#'
#' @return visualizations
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}.
#' @import ggplot2 gridExtra
#' @examples
#' df <- PEG(Nsim = 20, NSp = 2, living.rate = c("long", "short"),
#' surv.fixed = c(0.8, 0.65),
#' alpha.unif = TRUE, min.alpha = c(2, 1), max.alpha = c(3, 2),
#' pop.fixed = c(3605244, 55805898),
#' harvest.fixed = c(107802, 8447950),
#' Fs = c(0.1, 0.3, 0.5))
#'
#' SHI.plot(df)
SHI.plot <- function(data = NULL,
                     Sp = NULL,
                     IDSp = TRUE,
                     NameSp = NULL,
                     title = "Distribution of sustainability harvest index",
                     title.x = "Sustainability harvest index",
                     title.y = "Number of simulations",
                     color1 = "#00BFC4",
                     color2 = "#F8766D"
){

  ## I. Error messages ----
  # Print error messages if the vector types of the datas are not correct
  if(is.null(data)){ # If 'data' argument is null then print error message
    stop("'data' must be specified")
  }
  if (!is.data.frame(data)) { # If 'data' argument is not a data.frame then print error message
    stop("'data' must be a dataframe")
  }
  if(!is.null(Sp) && !is.numeric(Sp)){ # If 'Sp' argument is not null and is not a numeric vector then print error message
    stop("'Sp' must be a numeric vector")
  }
  if(!is.null(NameSp) && !is.character(NameSp)){ # If 'NameSp' argument is not null and is not a character string then print error message
    stop("'NameSp' must be a character string")
  }

  ## II. Creation of graphs ----

  ## Step 1 : Filter of datas ----
  # This step allows to create a graph for the species wanted
  if(!is.null(Sp)){ # If 'Sp' is not null then the data frame named 'data' is filtered according to 'Sp'
    dataFilterLogical <- with( # Create a data filter logical compararing the values of the column 'SP' of the data frame 'data' and the value of the object 'Sp'
      data,
      SP %in% Sp
    )
    newdata <- data[dataFilterLogical, ] # Extract the data frame for which 'dataFilterLogical' is TRUE
  } else { newdata <- data } # If 'Sp' is null then the data frame 'data' is renamed 'newdata'

  newdata$SP <- factor(newdata$SP) # Convert columns 'SP' of the data frame named 'newdata' to a factor vector



  ## Step 2 :  Right legend and facet boxes
  # This step allows to write the right legend for the facet boxes according to the PEG or PTL function used previously.
  # It also simplifies the creation of the facet boxes by rename the column "Fs" or "Fobj" in a unique column named "param".
  col <- which(colnames(data) == "Fs") # Search for a column named "Fs"



  ## Condition "if" : If the object 'col' is > 0 ----
  # There is a column named "Fs" in the data frame 'data' and the PEG function was used previously
  if(length(col)> 0){
    legend <- "Safety factor" # If the function used previously is PEG then the legend of the facet boxes is "Safety factor"
    colnames(newdata)[which(names(newdata)=="Fs")] <- "param" # Rename the column "Fs" by "param"



    ## Condition "else" : If the object 'col' is not > 0 ----
    # There is a column named "Fobj" in the data frame 'data' and the PTL function was used previously
  } else {
    legend <- "Objective factor" # If the function used previously is PTL then the legend of the facet boxes is "Objective factor"
    colnames(newdata)[which(names(newdata)=="Fobj")] <- "param" # Rename the column Fobj" by "param"
  }


  ## Step 3 : Data visualizations ----
  # This step allows to create the plot from a loop according to the level 'SP'
  for(i in levels(newdata$SP)){

    newdataT <- subset(newdata, newdata$SP == i) # Extract the data frame for which the value of the 'SP' column is equal to the level i. It allows to create a plot by 'SP'

    index <- as.numeric(i) # Convert the index i as a numeric vector for the next line
    Species <- NameSp[index] # Extract the name of the species study if it is given
    if (is.null(NameSp)){
      addtitle <- NULL
    } else {addtitle <- paste0(" of ",Species)} # Additionnal title for the title of the graphic with the name of the study species

    if (isTRUE(IDSp)){
      main.title <- paste0("SP", i, "\n",title, addtitle, "\n")
    } else { main.title <- paste0(title, addtitle, "\n") }

    nb_ticks <- function(n) {function(limits) pretty(limits, n)} # Function which allows to adapt the scale of the horizontal axis  according to the datas
    newdataT$high_rev <- factor((newdataT[,"SHI"]>1)*1) # The new column named "high_rev" allows to adapt the colour during the creation of the plot according to if SHI is lower or upper than 1


    theme_set(theme_classic()) # The theme of the plot

    print(
      ggplot(newdataT, aes_(x=~SHI, fill=~high_rev)) + # Plot with values of the column 'SHI' in horizontal axis from the data frame 'newdataT' and colorization according to the values of the column 'high_rev'
        geom_histogram(position="identity", colour="grey15", bins=50) + # Create an histogram
        facet_grid(param ~.) + # Facet according to the values of the column "Fs"
        geom_vline(aes(xintercept = 1), color="black", linetype="dashed", size = 0.7) + # Dotted line at 1
        labs(title= main.title , x=paste0("\n", title.x), y= paste0(title.y,"\n"), subtitle=legend) + # Main title, legend title and axes title
        scale_fill_manual(values=c(color1, color2)) + # Change the graph colours according to the values of "high_rev"
        scale_x_continuous(breaks=nb_ticks(10)) + # Change the scale of the horizontal axis
        theme(plot.subtitle = element_text(hjust = 1.05), # Position of facets legend
              panel.border = element_blank(), # No borders for facet boxes
              panel.grid.minor.y = element_blank(), # No background graduation
              panel.grid.major.x = element_blank(), # No background graduation
              panel.grid.minor.x = element_blank(), # No background graduation
              panel.spacing = unit(1.2, "lines"), # Change the spacing between facets
              plot.title = element_text(colour = "grey15", hjust = 0.5, size = 16), # Changes the title plot
              axis.text = element_text(size=11, color = "grey15"), # Changing the axes text
              axis.title = element_text(size=13, color = "grey15"), # Changing the axes title
              legend.background = element_blank(), # No background for the legend
              legend.position="none" , # No legend
              strip.background = element_rect(color="grey80", fill="white", size=0.5, linetype="solid"), # Change the facets box
              strip.text = element_text(color = "grey15", size = 11)) # Change the facets text
    )

  } # End of loop i


} # End of function
