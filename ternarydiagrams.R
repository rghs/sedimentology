# Making sedimentary ternary diagrams

library(ggplot2)
library(ggtern)

#' buildtern
#' 
#' Constructs ternary diagrams, with presets based on the standard sedimentary
#' ternary diagrams after Dickinson et al. (1983)
#'
#' @param dataset Dataframe containing recalculated modes
#' @param x_col Character of column name containing data for bottom left ternary (e.g., feldspar)
#' @param y_col Character of column name containing data for top ternary (e.g., quartz)
#' @param z_col Character of column name containing data for bottom right ternary (e.g., lithics)
#' @param categories Character of column name containing classifying variable to shade points by (e.g. locality); optional
#' @param fields Character specifying lines on diagram, select 'qtfl' or 'qmflt'; optional
#' @param limits Vector containing limits for ternary i.e. c(1,0.2,0.2); optional
#'
#' @return tern: ggplot2 plot
#' @export
#'
#' @examples
#' buildtern(df, 'f', 'qt', 'l', 'locality', 'qtfl', c(1,0.2,0.2))
buildtern = function(dataset, x_col, y_col, z_col, categories = NULL, fields = NULL, limits = NULL){
  if(is.character(x_col) == FALSE | is.character(y_col) == FALSE | is.character(z_col) == FALSE){
    stop('Column names must be passed as characters.')
  }
  
  cols = c('diagram','x_name','y_name','z_name','t_arr','l_arr','r_arr')
  diags = c('qtfl','qmflt','qmpk','qplvls')
  ynames = c('Qt','Qm','Qm','Qp')
  xnames = c('F','F','P','Lv')
  znames = c('L','Lt','K','Ls')
  tarrs = c('Qt = Qm + Qp','Qm','Qm','Qp')
  larrs = c('F = P + K','F = P + K','P','Lv')
  rarrs = c('L = Lv + Ls','Lt = L + Qp','K','Ls')
  type = data.frame(diags, xnames, ynames, znames, tarrs, larrs, rarrs)
  colnames(type) = cols
  type = type[type$diagram == fields,]
  
  tern = ggtern(data = dataset, aes_string(x=x_col, y=y_col, z=z_col, color=categories)) +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(tern.panel.grid.minor = element_blank(),
          tern.panel.grid.major = element_blank(),
          tern.axis.text = element_text(size = 8))
  
  if(fields == 'qtfl'){
    suppressWarnings({tern = tern + geom_segment(aes(x=0, y=0.97, z= 0.03, xend=0.85, yend=0, zend=0.15), linetype = 1, color = "black") + # Continental block zone
      geom_segment(aes(x=0.4, y=0.512, z=0.087, xend=0, yend=0.25, zend=0.75), linetype = 1, color = "black") + # Recycled orogen zone
      geom_segment(aes(x=0.18, y=0.82, z=0., xend=0.158, yend=0.799, zend=0.054), linetype = 2, color = "black") + # Craton interior
      geom_segment(aes(x=0.45, y=0.55, z=0., xend=0.4, yend=0.512, zend=0.087), linetype = 2, color = "black") + # Transitional continental
      geom_segment(aes(x=0.705, y=0.175, z=0.134, xend=0.128, yend=0.333, zend=0.543), linetype = 2, color = "black") + # Dissected arc
      geom_segment(aes(x=0.5, y=0., z=0.5, xend=0., yend=0.25, zend=0.75), linetype = 2, color = "black") # Transitional arc
    })
  } else if(fields == 'qmflt'){
    suppressWarnings({tern = tern + geom_segment(aes(x=0., y=0.89, z=0.11, xend=0.77, yend=0., zend=0.23), linetype = 1, color = "black") + # Continental block
      geom_segment(aes(x=0.191, y=0.671, z=0.142, xend=0.13, yend=0., zend=0.87), linetype = 1, color = "black") + # Recycled orogen
      geom_segment(aes(x=0.349, y=0.491, z=0.165, xend=0.16, yend=0.294, zend=0.55), linetype = 1, color = "black") + # Mixed arc
      geom_segment(aes(x=0.2, y=0.8, z=0., xend=0.132, yend=0.742, zend=0.134), linetype = 2, color = "black") + # Craton interior
      geom_segment(aes(x=0.167, y=0.498, z=0.331, xend=0., yend=0.58, zend=0.42), linetype = 2, color = "black") + # Quartzose recycled
      geom_segment(aes(x=0.142, y=0.211, z=0.639, xend=0., yend=0.29, zend=0.71), linetype = 2, color = "black") + # Transitional recycled
      geom_segment(aes(x=0.43, y=0.57, z=0., xend=0.349, yend=0.491, zend=0.165), linetype = 2, color = "black") + # Transitional continental
      geom_segment(aes(x=0.584, y=0.219, z=0.204, xend=0.16, yend=0.294, zend=0.55), linetype = 2, color = "black") + # Dissected arc
      geom_segment(aes(x=0.47, y=0., z=0.53, xend=0.144, yend=0.127, zend=0.731), linetype = 2, color = "black") # Transitional arc
    })
  }
  tern = tern + geom_point()
  if(is.null(limits) == FALSE){
    tern = tern + tern_limits(limits[1], limits[2], limits[3]) # Order: top, left, right
  }
  tern = tern + labs(x = type[1,2],
                     y = type[1,3],
                     z = type[1,4],
                     Tarrow = type[1,5],
                     Larrow = type[1,6],
                     Rarrow = type[1,7])
  
  return(tern)
}