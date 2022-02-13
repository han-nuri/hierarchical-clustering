# Unsupervised ML Hierarchical clustering

#packages
pacman::p_load(data.table, dplyr, ggplot2, 
               readr, stringr, bannerCommenter,
               cluster, R6, tidyverse, gridExtra, stargazer,
               ggdendro, pdfCluster, mclust, epitools, ggpubr, stringi, stringr, grid)

############################################################################
############################################################################
###                                                                      ###
###                              SECTION 0:                              ###
###                          DATA PREPROCESSING                          ###
###                                                                      ###
############################################################################
############################################################################

#load party manifesto data
#data <- fread("https://manifesto-project.wzb.eu/down/data/2021a/datasets/MPDataset_MPDS2021a.csv", encoding = "UTF-8")
#write_csv(data, file = "party_manifesto.csv")
data <- fread(file = "party_manifesto.csv", encoding = "UTF-8")

#only take into account the parties that have at least 5% electoral votes
data = data[pervote > 5]

#create new id variable that combines party name and countryname
data$id <- paste(data$partyname, data$countryname, sep = " - ")

# for each party, take only one coding year
data <- data %>% group_by(party) %>% slice(., which.max(coderyear))

# now select unique variable "id"
data <- data %>% group_by(id) %>% slice(., which.max(pervote))

#only take relevant variable columns
data <- data[,c(1,2,7,8,9,10, 13, 26:81, 175)]


X <- data[,8:63]
rnames <- data[,"id"]
rownames(X)[] <- rnames$id

summary_tably <- X %>% map_df(.f = ~ broom::tidy(summary(.x)), .id = "variable")

row.names(summary_tably)[] <- summary_tably$variable


###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 1:                             ###
###                           DATA INSPECTION                           ###
###                                                                     ###
###########################################################################
###########################################################################


# pdf("data_gridExtra.pdf")       # Export PDF
# grid.table(summary_tably)
# dev.off()
# g <- tableGrob(summary_tably[c(seq(1,56,2)),])
# 
# grid.arrange(g)

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 2:                             ###
###                   CREATING DISSIMILARITY MATRICES                   ###
###                                                                     ###
###########################################################################
###########################################################################



##----------------------------------------------------------------
##                    Section 2.1: MINKOWSKI                    -- requires normalization
##----------------------------------------------------------------


minkowski_dist <- dist(X, method = "minkowski", p = 1.5)

##----------------------------------------------------------------
##                   Section 2.2: MALAHANOBIS                   -- does not require normalization
##----------------------------------------------------------------

# mahala <- function(X) {
#   n <- nrow(X)
#   C <- cov(X)
#   P <- t(X)
#   M <- matrix(NA, nrow = n, ncol = n)
#   for (i in 1:n) {
#     for (j in 1:n) {
#       M[i,j] <-
#         sqrt((X[i,]-X[j,])%*%solve(C)%*%(P[,i]-P[,j]))
#     }
#   }
#   return(M)
# }
# 
# X <- as.matrix(X)
# M_1 <- mahala(X)
# 
# save(M_1, file = "M_matrix_1percent.rds")

# load mahalanobis distance matrix
load("M_matrix.rds")
load("M_matrix_5percent.rds")
load("M_matrix_1percent.rds")
# if five percent
M = M_5
#if one percent
#M = M_1
row.names(M) <- rnames$id
colnames(M) <- rnames$id

# option 1
mahalanobis_dist <- as.dist(M)

##----------------------------------------------------------------
##               Section 2.3: EUCLIDIAN DISTANCES               --
##----------------------------------------------------------------

# Section 2.3.1: Create Euclidian Distance Matrix
euclid_dist <- dist(X, method = "euclidean")


##---------------------------------------------------------------
##                 Section 2.4: GOWER DISTANCE                 --
##---------------------------------------------------------------
gower_dist <- daisy(X, metric = "gower")
gower_dist <- as.dist(gower_dist)

# Section 2.5 Manhattan distance
manhattan_dist <- dist(X, method = "manhattan")

# Section 2.6 Canberra distance
canberra_dist = dist(X, method = "canberra")

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 3:                             ###
###                       HIERARCHICAL CLUSTERING                       ###
###                                                                     ###
###########################################################################
###########################################################################
##----------------------------------------------------------------
##              Section 3.2: AGGLOMERATIVE NESTING              --
##----------------------------------------------------------------

# for optimal par.method values, see website: https://www.davidzeleny.net/anadat-r/doku.php/en:hier-agglom_examples

agnes_models = R6Class(classname = "agnes_models",
                       private = list(
                         ..x = NULL
                       ),
                       public = list(
                              single = NULL,
                              complete = NULL,
                              average = NULL,
                              ward = NULL,
                              flexible = NULL,
                              gaverage = NULL,
                              distance_matrix = NULL,
                              initialize  = function(x) { # this function creates a new object with the function new()
                                
                                self$single <- agnes(x, diss = F, method = "single", #single
                                                     keep.diss = F, keep.data = T, trace.lev = 1)
                                
                                self$complete <- agnes(x, diss = F, method = "complete", #complete
                                                       keep.diss = F, keep.data = T, trace.lev = 1)
                                
                                self$average <- agnes(x, diss = inherits(x, "dist"), method = "average", par.method, #average
                                                       keep.diss = F, keep.data = T, trace.lev = 1)
                                
                                self$ward <- agnes(x, diss = F, method = "ward", #ward
                                                       keep.diss = F, keep.data = T, trace.lev = 1)

                                self$distance_matrix <- deparse(substitute(x)) #deparse substitute trick
                                
                                private$..x = x
                              },
                              
                              print = function() { # this function tells me what kind of distance matrix the object is
                                cat("Cluster Model: \n")
                                cat("   Distance Matrix:   ", self$distance_matrix, sep = "")
                              },
                              
                              plot = function(type, 
                                              xlim = NULL, 
                                              ylim = NULL,
                                              horiz = NULL) {
                                if (type == "single") {                         # single
                                  hcd = as.dendrogram(self$single)
                                } else if (type == "complete") {              # complete
                                  hcd = as.dendrogram(self$complete)
                                } else if (type == "average") {              # average
                                  hcd = as.dendrogram(self$complete)
                                } else if (type == "ward") {              # ward
                                  hcd = as.dendrogram(self$complete)
                                } else {
                                  print("no linkage method specified")
                                }
                              
                                plot(hcd, type = "rectangle", ylab = "Height", 
                                     xlim = xlim, 
                                     ylim = ylim, 
                                     horiz = horiz)
                              },
                              
                              class = T))



euclid_agnes = agnes_models$new(euclid_dist)
gower_agnes = agnes_models$new(gower_dist)
mahalanobis_agnes = agnes_models$new(mahalanobis_dist)
minkowski_agnes = agnes_models$new(minkowski_dist)
manhattan_agnes = agnes_models$new(manhattan_dist)
canberra_agnes = agnes_models$new(canberra_dist)

###########################################################################
###########################################################################
###                                                                     ###
###                    SECTION 4: CLUSTER EVALUATION                    ###
###                                                                     ###
###########################################################################
###########################################################################

##----------------------------------------------------------------
##                   Section 4.1: CH FUNCTION                   --
##----------------------------------------------------------------

# ETHEN LIU's CHCriterion
Distance <- function(cluster)
{
  # the center of the cluster, mean of all the points
  center <- colMeans(cluster)
  
  # calculate the summed squared error between every point and 
  # the center of that cluster 
  distance <- apply(cluster, 1, function(row)
  {
    sum( ( row - center )^2 )
  }) %>% sum()
  
  return(distance)
}

# calculate the within sum squared error manually for hierarchical clustering 
# [WSS] : pass in the dataset, and the resulting groups(cluster)

WSS <- function( data, groups )
{
  k <- max(groups)
  
  # loop through each groups (clusters) and obtain its 
  # within sum squared error 
  total <- lapply( 1:k, function(k)
  {
    # extract the data point within the cluster
    cluster <- subset( data, groups == k )
    
    distance <- Distance(cluster)
    return(distance)
  }) %>% unlist()
  
  return( sum(total) )
}


CHCriterion <- function( data, kmax, clustermethod, dist = dist, ...)
{
  if( !clustermethod %in% c( "kmeanspp", "hclust" ) )
    stop( "method must be one of 'kmeanspp' or 'hclust'" )
  
  # total sum squared error (independent with the number of cluster k)
  tss <- Distance( cluster = data )
  
  # initialize a numeric vector storing the score
  wss <- numeric(kmax)
  
  # k starts from 2, cluster 1 is meaningless
  if( clustermethod == "kmeanspp" )
  {
    for( k in 2:kmax )
    {
      results <- Kmeanspp( data, k, ... )
      wss[k]  <- results$tot.withinss 
    }		
  }else # "hclust"
  {
    d <- dist
    clustering <- hclust( d, ... )
    for( k in 2:kmax )
    {
      groups <- cutree( clustering, k )
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }		
  
  # between sum of square
  bss <- tss - wss[-1]
  
  # cluster count start from 2! 
  numerator <- bss / ( 1:(kmax-1) )
  denominator <- wss[-1] / ( nrow(data) - 2:kmax )
  
  criteria <- data.frame( k = 2:kmax,
                          CHIndex = numerator / denominator)
  
  # convert to long format for plotting 
  criteria_long <- gather( criteria, "index", "value", -1 )
  
  plot <- ggplot( criteria_long, aes( k, value, color = index ) ) + 
    geom_line() + geom_point( aes( shape = index ), size = 3 ) +
    facet_wrap( ~ index, scale = "free_y" ) + 
    guides( color = FALSE, shape = FALSE ) + 
    ggtitle("Calinski-Harabasz Criterion") +
    ylab("CH Criterion Score") +
    xlab("Number of Clusters k") +
    geom_vline(xintercept = criteria$k[which.max(criteria$CHIndex)], linetype = "dashed") +
    geom_text(aes(x = criteria$k[which.max(criteria$CHIndex)], 
                  y = max(criteria$CHIndex), 
                  label = paste("k = " , criteria$k[which.max(criteria$CHIndex)])), 
                  vjust = 1.8, hjust = -1, size = 4, color = "black")
  
  return( list( data = criteria, 
                plot = plot ) )
}

# FIRST CHECK: WHICH LINKAGE METHOD TO USE
plot_ac = function(agnes, distance) {
  ac_data = rbind(agnes$single$ac, agnes$average$ac, agnes$complete$ac, agnes$ward$ac)
  colnames(ac_data)[1] <- "score"
  ac_data = as.data.frame(ac_data, row.names = c("single", "average", "complete", "ward"))
  
  #plot
  ggplot(ac_data, aes(x = rownames(ac_data), y = score, group = 1)) + geom_point() +
    geom_line() +
    theme_bw() + ggtitle(paste("Agglomerative Coefficients", distance ,"Distance")) +
    xlab("Linkage Methods") +
    ylab("Agglomerative Coefficient") + geom_text(aes(label = round(score, 3)), vjust = 1.2, hjust = 1.2, size = 3.2)
  
}

ac1 = plot_ac(euclid_agnes, "Euclid")
ac2 = plot_ac(gower_agnes, "Gower")
ac3 = plot_ac(mahalanobis_agnes, "Mahalanobis")
ac4 = plot_ac(minkowski_agnes, "Minkowski")
ac5 = plot_ac(manhattan_agnes, "Manhattan")
ac6 = plot_ac(canberra_agnes, "Canberra")

ggarrange(ac1, ac2, ac3, ac4, ac5, ac6, ncol = 2, nrow = 3)

# ==> use Ward linkage method


# Second Check: Which distance to use
plot_CH = function(x, distance, dist) {
  # inject inner function into outer function
  environment(get_distance_name) <- environment()
  data = CHCriterion(x$ward$data, kmax = 30, clustermethod = "hclust", dist = dist)
  data$plot + ggtitle(paste(distance, "Distance"))
}

ch1 = plot_CH(euclid_agnes, "Euclid", dist = euclid_dist)
ch2 = plot_CH(gower_agnes, "Gower", dist = gower_dist)
ch3 = plot_CH(mahalanobis_agnes, "Mahalanobis", dist = mahalanobis_dist)
ch4 = plot_CH(minkowski_agnes, "Minkowski", dist = minkowski_dist)
ch5 = plot_CH(manhattan_agnes, "Manhattan", dist = manhattan_dist)
ch6 = plot_CH(canberra_agnes, "Canberra", dist = canberra_dist)


ggarrange(ch1, ch2, ch3, ch4, ch5, ch6, ncol = 2, nrow = 3)


# Discussion:
#

#==> Mahalanobis is the best when party threshold set to 10 percent OR 1 percent, 5 clusters THus: Use Mahalanobis Distance with Ward Linkage Algorithm
# ==> Euclidean distance is best when party threshold is set to 1 percent. use in this case euclidean.



############################################################################
############################################################################
###                                                                      ###
###                              SECTION 5:                              ###
###                               PLOTTING                               ###
###                                                                      ###
############################################################################
############################################################################

# plot function

colored_dendrogram = function(distance = NULL, n = NULL) {
  hc <- hclust(distance, "ward.D2")
  dendr    <- dendro_data(hc, type="rectangle") # convert for ggplot
  clust    <- cutree(hc,k=n)                    # find 2 clusters
  clust.df <- data.frame(label=names(clust), cluster=factor(clust))
  
  # dendr[["labels"]] has the labels, merge with clust.df based on label column
  dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
  # plot the dendrogram; note use of color=cluster in geom_text(...)
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=label(dendr), aes(x, y, label=label, hjust=0.1, color=cluster), 
              size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank()) +
    ggtitle("Cluster Groupings")
}


colored_dendrogram(distance = canberra_dist, n = 6)
asdf <- colored_dendrogram(distance = canberra_dist, n = 6)

lol = asdf$plot_env$dendr$labels

#merging cluster information to main dataset
data[lol, cluster := i.cluster, on = c("id" = "label")]

#merging percentage vote info
data2 <- fread(file = "party_manifesto.csv", encoding = "UTF-8")

#only take into account the parties that have at least 5% electoral votes
data2 = data2[pervote > 5]

#create new id variable that combines party name and countryname
data2$id <- paste(data2$partyname, data2$countryname, sep = " - ")

# for each party, take only one coding year
data2 <- data2 %>% group_by(party) %>% slice(., which.max(coderyear))

# now select unique variable "id"
data2 <- data2 %>% group_by(id) %>% slice(., which.max(pervote))

setDT(data2)

#merge here
data[data2, pervote := pervote, on = c("party")]



asdf + 
  ylim(100, 0) + geom_text(data = label(dendr),
                          aes(x, y, label = label, size = 2, hjust = -0.2))

clusterings = list()
test = list()


test_test = data %>% group_by(cluster) %>% arrange(., pervote, .by_group = T) %>% slice_max(., n = 15, order_by = pervote)

for (i in 1:6) {
  set.seed(333)
  clusterings[[i]] = test_test$id[test_test$cluster == i]
}

df <- data.frame(matrix(unlist(clusterings), nrow=length(clusterings), byrow= T))

lololol = as.data.frame(t(df))
colnames(lololol)[] <- c(paste("cluster", c(1:6), sep = "_"))
rownames(lololol)[] <- c(1:15)

# VISUALIZING THE PARTIES IN TEXT

lololol = sapply(lololol, function(x) stri_trans_general(x, "ascii"))

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1.2)),
  colhead = list(fg_params=list(cex = 1)),
  rowhead = list(fg_params=list(cex = 0.5)))



pdf("data_gridExtra.pdf")       # Export PDF
grid.table(lololol)
dev.off()
g1 <- tableGrob(lololol[,c(1:2)], theme = mytheme)
g2 <- tableGrob(lololol[,c(3:4)], theme = mytheme)
g3 <- tableGrob(lololol[,c(5:6)], theme = mytheme)

g1$widths <- unit(rep(1/ncol(g2), ncol(g2)), "npc")
g2$widths <- unit(rep(1/ncol(g2), ncol(g2)), "npc")
g3$widths <- unit(rep(1/ncol(g2), ncol(g2)), "npc")

grid.arrange(g3, ncol = 1)

###########################################################################
###########################################################################
###                                                                     ###
###               SECTION 6: COMPARISON WITH LABELED DATA               ###
###                                                                     ###
###########################################################################
###########################################################################

# Mahalanobis distance
ARI = rep(NA, 20)

for (i in 1:20) {
  # Mahalanobis
  mahala <- colored_dendrogram(canberra_dist, n = i)
  
  # first merge cluster labelling with original data: In the original data, there are 12 parfams with partythreshold set to 10%
  setDT(data)
  
  # merge_cluster
  merge_cluster = function(plot_object, data) {
    cluster_data = plot_object$plot_env$dendr$labels
    setDT(cluster_data)
    data[cluster_data, cluster := i.cluster, on = c("id" = "label")]
  }
  
  newly_labeled <- merge_cluster(mahala, data)
  
  #make labels to factors
  newly_labeled$parfam = as.factor(newly_labeled$parfam)
  newly_labeled$cluster = as.factor(newly_labeled$cluster)
  
  ARI[i] = adjustedRandIndex(newly_labeled$parfam, newly_labeled$cluster)
}


test = as.data.frame(ARI)
test$N = c(1:20)

ggplot(test, aes(x = N, y = ARI, group = 1)) + geom_point() +
  geom_line() +
  theme_bw() + ggtitle(paste("Rand Index")) +
  xlab("Number of Clusters in Canberra Model") +
  ylab("Score") + geom_text(aes(label = round(ARI, 3)), vjust = 1.8, hjust = 1.5, size = 3.2)



# visualize ARI
test = as.data.frame(ARI)
test$N = c(1:20)

ggplot(test, aes(x = N, y = ARI, group = 1)) + geom_point() +
  geom_line() +
  theme_bw() + ggtitle(paste("Rand Index")) +
  xlab("Number of Clusters") +
  ylab("Score") + geom_text(aes(label = round(ARI, 3)), vjust = 1.8, hjust = 1.5, size = 3.2)
