#bioconductor library is required for packages
library(scales)
library(spatial)
library(RColorBrewer)
library(viridis)
library(plotrix)
library(factoextra)
#############################################################################
##############################DEFINE LOCATIONS###############################
#Generate new plot
plot.new()
op <- par(mfrow = c(2, 2))
#############Load in data###########

# input file should be a CSV with cluster ID as column 1, 
#subsequent columns represent each marker and each row a new cluster
#e.g.
#clusterID   CD1   CD2   CD3
# 001        30     3     20
# 002        20     40    5

phe_file <- file.choose() 
all_data <- read.csv(phe_file)

#create marker correlation matrix

phenotype = all_data[,-1]
cluster_correlations = cor(phenotype, method = "spearman")

data_set <<- cluster_correlations
data_set = data_set[, 2:length(data_set[1, ])] 

#defining scale function
range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
} 
#apply scale function
data_set = apply(data_set, 2, range01) # applying scale function
data_set = as.data.frame(data_set)
k.max <- length(all_data[1,])-1
data <- data_set

k = 4 #manually assign K based on elbow plot
set.seed(1)

#re-run Kmeans with selected K
km = kmeans(data, k) 

#defining parameters
cluster_length = length(km$cluster[which(km$cluster == 1)])
num_parameters = length(data_set[1,])
# use km$cluster to define colours
broad_colours <<- km$cluster*100
broad_colours <<- as.data.frame(broad_colours) #defining clusters as a colour
cluster_types = unique(broad_colours)
number_clusters = length(cluster_types$broad_colours)
palette_2 = inferno((100*(number_clusters+1))) #defining palette


#applying colour range within a cluster
for(cluster_group in cluster_types$broad_colours) {
  iteration = 1
  name = cluster_group
  for (i in 1:length(broad_colours$broad_colours)) {
    if (broad_colours$broad_colours[i] == name) {
      temp_colour_label <- broad_colours$broad_colours[i]
      print(broad_colours$broad_colours[i])
      broad_colours$broad_colours[i] <- temp_colour_label + (as.integer(50 / (num_parameters/(k))) * iteration) 
      iteration = iteration + 1
    }
  }
} 

#performing principal component analysis
pc = (princomp(data)) 

#plotting PCA variation curve
plot(pc, type = "l")
x = pc$scores[, 1]
y = pc$scores[, 2]

#plotting initial locations
plot(
  x,
  y,
  col = alpha(palette_2[km$cluster*100], 0.8),
  pch = 19,
  xlab = "",
  ylab = "",
  axes = FALSE,
  cex = 2
) 

#plotting PCA variance contribution and marker location

colours_list = vector(mode = "list", length = k)
for(i in 1:k){
  colours_list[i] <-palette_2[100*i]
}

fviz_pca_ind(pc,col.ind = km$cluster,gradient.cols = colours_list, repel = TRUE)

print("Markers Clustered")


############################################################################
##############################CREATE BOXES##################################
#Generate new plot

plot.new()
op <- par(mfrow = c(2, 2))

print("cluster data assigned")

#scaling data
scaled_data <- all_data[,-1]
scaled_data = as.matrix(scaled_data)
markers <- colnames(scaled_data)

############################defining functions##############################

#make brick function - a function to define the size, colour and position of a given brick
make_brick = function(edge,
                      place_brick = FALSE,
                      i = 1, adjust = FALSE, adjusted_a, adjusted_b) {
  
  start_a = values[i, 3] * 50 #alter value to change box spread (default is 50)
  start_b = values[i, 4] * 50 # this one two and update line 159
  edge = edge
  a = start_a
  b = start_b
  width = a + edge
  height = b + edge
  list_positions = c(a, b, width, height)
  test_box <<- list_positions
  if (place_brick == TRUE) {
    rect(a, b, width, height, col = alpha(palette[cols], 0.5))
  }
  if (adjust == TRUE){
    a = adjusted_a
    b = adjusted_b
    width = a + edge
    height = b + edge
    list_positions = c(a, b, width, height)
    test_box <<- list_positions
  }
  return(test_box)
}

print("brick function defined")

#place brick function - a function to place the pre-defined brick on the plot

place_brick <- function(a, b, width, height, cols) {
 # print("drawing Brick")
  rect(a,
       b,
       width,
       height,
       col = alpha(palette_2[cols], 0.6),
       border = NA)
}

print("place brick function defined")

#check overlap function

brick_list <<- matrix(data = -100, length(scaled_data[1,])+1, 4)

check_overlap <- function(brick, current_bricks) {
  print("checking overlap")
  overlap_score = 0
  for (i in 1:length(current_bricks[,1 ])) {
    if (((brick[1] <= current_bricks[i, 3] && brick[3] >= current_bricks[i, 1])
         ||
         (brick[3] >= current_bricks[i, 1] && brick[1] <= current_bricks[i, 3]))
    &&
    ((brick[2] <= current_bricks[i, 4] && brick[4] >= current_bricks[i, 2])
     ||
     (brick[4] >= current_bricks[i, 2] && brick[2] <= current_bricks[i, 4])
    ))
    {
      print("overlap detected")
      overlap_score = overlap_score + 1
      overlapped_a <<- current_bricks[i,1]
      overlapped_height <<- current_bricks[i,4]
    }
  }
  print(overlap_score)
  if(overlap_score > 0){return(TRUE)} else{return(FALSE)} #change TRUE to FALSE in this line to stop overlap adjustment
}

####THIS LOGIC IS NOT QUITE WORKING CORRECTLY

##################Generating boxes for each cluster and labelling############################

for (cluster in 1:length(scaled_data[, 1])) {
  print(paste("plotting cluster", cluster, sep = " "))

  values = (scaled_data[cluster,])
  for (i in 1:length(values)) {
    if (values[i] <= 1.0) {  ############ Selecting minium cutoff value change the 1 to alter this value.
      values[i] <- 0
    }
  }
  values = as.data.frame(values)
  colours = seq(1, length(scaled_data[1,]), 1)
  colours = as.data.frame(colours)
  rownames(colours) = rownames(values)
  broad_colours = as.data.frame(broad_colours)
  rownames(broad_colours) = rownames(values)
  values = cbind(values, colours, x, y,broad_colours)
  
  test_box <<- matrix(data = NA, 1, 4)
  brick_list <<- matrix(data = -100, length(scaled_data[1,])+1, 4)
  plot(
    range(x)*100,#if changed spread, also adjust this value by the same fold (default is 100)
    range(y)*150,
    type = "n",
    xlab = paste("Cluster", all_data[cluster, 1]),
    ylab = "",
    axes = FALSE
  )
  for (i in 1:length(values[, 1])) {
    edge = sqrt(values[i, 1]) * (((range(x)[2] * 100) - (range(x)[1] * 100)) /
                                   20)
    make_brick(edge = edge, FALSE, i = i)
    brick_list[i, 1] = test_box[1]
    brick_list[i, 2] = test_box[2]
    brick_list[i, 3] = test_box[3]
    brick_list[i, 4] = test_box[4]
    
    if ((test_box[3] - test_box[1]) >= 0.6 * (((range(x)[2] * 100) - (range(x)[1] *
                                                                      100)) / 20)) {
 
      check_overlap(brick_list[i, ], brick_list[-i,])
      while (check_overlap(brick_list[i, ], brick_list[-i,]) == TRUE) {
        print("defining alternate brick")
        print(paste(test_box,rownames(values)[i], sep = " " ))
        brick_list[i, ] =0
        make_brick(edge = edge,
                   FALSE,
                   i = i,
                   adjust = TRUE, adjusted_a = overlapped_a+0.1, adjusted_b = overlapped_height+0.1) 
        brick_list[i, 1] <- test_box[1]
        brick_list[i, 2] <- test_box[2]
        brick_list[i, 3] <- test_box[3]
        brick_list[i, 4] <- test_box[4]
      }
      
      place_brick( brick_list[i, 1],  brick_list[i, 2], brick_list[i, 3],brick_list[i, 4], values[i, 5])
      center_x = (( brick_list[i, 1]) + (brick_list[i, 3])) / 2
      center_y = ( brick_list[i, 2] + brick_list[i, 4]) / 2
      text(center_x,
           center_y,
           labels = rownames(values)[i],
           cex = 0.7)
    }
  }
}