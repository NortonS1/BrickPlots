library(ComplexHeatmap)
library(dendextend)
library(scales)
library(spatial)
library(RColorBrewer)
library(viridis)
library(plotrix)
#############################################################################
##############################DEFINE LOCATIONS###############################
plot.new()
op <- par(mfrow = c(2, 2))

cor_mat_file <- file.choose()
data_set <<- read.csv(cor_mat_file)
data_set = data_set[, 2:length(data_set[1, ])] 
range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

data_set = apply(data_set, 2, range01)
data_set = as.data.frame(data_set)
k.max <- 30
data <- data_set
wss <- sapply(1:k.max, function(k) {
  kmeans(data, k)$tot.withinss
})

plot(
  1:k.max,
  wss,
  type = "b",
  pch = 19,
  frame = FALSE,
  xlab = "Number of clusters K",
  ylab = "Total within-clusters sum of squares"
)
k = 3
km = kmeans(data, k)
cluster_length = length(km$cluster[which(km$cluster == 1)])
# USE km$cluster to define colours
broad_colours = km$cluster*100
broad_colours = as.data.frame(broad_colours)

cluster_types = unique(broad_colours)
number_clusters = length(cluster_types$broad_colours)
palette_2 = magma((100*(number_clusters+1)))

for(cluster_group in cluster_types$broad_colours) {
  iteration = 1
  name = cluster_group
  for (i in 1:length(broad_colours$broad_colours)) {
    # cluster_length = length(broad_colours$broad_colours[which(broad_colours$broad_colours ==broad_colours$broad_colours[i])])
    if (broad_colours$broad_colours[i] == name) {
      broad_colours$broad_colours[i] <<- broad_colours$broad_colours[i] + (as.integer(60 / length(markers)) * iteration) #FIX THIS LINE
      iteration = iteration + 1
    }
  }
}
range(x)
pc = (princomp(data))
plot(pc, type = "l")
x = pc$scores[, 1]
y = pc$scores[, 2]
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

print("Markers Clustered")


############################################################################
##############################CREATE BOXES##################################


plot.new()
op <- par(mfrow = c(2, 2))

phe_file <- file.choose()

all_data <- read.csv(phe_file)

print("cluster data assigned")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

scaled_data <- apply(all_data[,-1], 2, scale)
scaled_data = scaled_data
markers <- colnames(scaled_data)
make_brick = function(edge,
                      place_brick = FALSE,
                      i = 1) {
  print("making brick")
  start_a = values[i, 3] * 50
  start_b = values[i, 4] * 50
  edge = edge
  a = start_a
  b = start_b
  width = a + edge
  height = b + edge
  list_positions = c(a, b, width, height)
  print("positions defined")
  test_box <<- list_positions
  if (place_brick == TRUE) {
    print("drawing Brick")
    rect(a, b, width, height, col = alpha(palette[cols], 0.5))
  }
  return(test_box)
}

print("brick function defined")

place_brick <- function(a, b, width, height, cols) {
  print("drawing Brick")
  rect(a,
       b,
       width,
       height,
       col = alpha(palette_2[cols], 0.8),
       border = NA)
}

print("place brick function defined")

for (cluster in 1:length(scaled_data[, 1])) {
  print(paste("plotting cluster", cluster, sep = " "))
  values = (scaled_data[cluster,])
  for (i in 1:length(values)) {
    if (values[i] < 0) {
      values[i] <- 0
    }
  }
  values = as.data.frame(values)
  print("1")
  colours = seq(1, length(scaled_data[1,]), 1)
  print("2")
  colours = as.data.frame(colours)
  print("3")
  rownames(colours) = rownames(values)
  print("4")
  broad_colours = as.data.frame(broad_colours)
  print("5")
  rownames(broad_colours) = rownames(values)
  print("6")
  values = cbind(values, colours, x, y,broad_colours)
  print("7")
  
  positions <<- matrix(data = NA, length(scaled_data[1,]), 4)
  positions[1,] = c(0, 0, 0, 0)
  test_box <<- matrix(data = NA, 1, 4)
  
  plot(
    range(x)*100,
    range(y)*100,
    type = "n",
    xlab = paste("Cluster", all_data[cluster, 1]),
    ylab = "",
    axes = FALSE
  )
  for (i in 1:length(values[, 1])) {
    edge = sqrt(values[i, 1]) * (((range(x)[2]*100)-(range(x)[1]*100))/20)
    make_brick(edge = edge, FALSE, i = i)
    Number_bricks <<- colSums(!is.na(positions))[1]

    if ((test_box[3] - test_box[1]) >= 0.6*(((range(x)[2]*100)-(range(x)[1]*100))/20)) {
      place_brick(test_box[1], test_box[2], test_box[3], test_box[4], values[i, 5])
      center_x = (test_box[1] + test_box[3]) / 2
      center_y = (test_box[2] + test_box[4]) / 2
      text(center_x,
           center_y,
           labels = rownames(values)[i],
           cex = 0.7)
    }
  }
}


