# define function to remove results with certain flags.
flag <- function(data) 
  {data[!grepl('CurveFitFail|MultipleMeltPeak|NoAmplification',
                                     data$flags),]}

# define function to strip a prefix from a string.
prefix_strip <- function(str, prefix) 
{substr(str, 1 + nchar(prefix), nchar(str))}