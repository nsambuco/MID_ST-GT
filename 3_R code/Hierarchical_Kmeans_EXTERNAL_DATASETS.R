# ANALYSIS 3: HIERARCHICAL K-MEANS CLUSTERING - REPLICATION COHORT 
# Import the dataset
df_ext <- as.data.frame(readxl::read_excel("External_dataset.xlsx"))
str(df_ext)

# Scaling -1 1 Replication
df <- df_ext[,4:25]
str(df)
df <- as.data.frame(lapply(df, rescale, to = c(-1, 1)))
df <- round(df, 3)
head(df)

# load the Discovery dataset
df_discovery <- as.data.frame(readxl::read_excel("dataset_clustering.xlsx"))
str(df_discovery)

# Scaling -1 1 Discovery
data <- as.data.frame(lapply(df_discovery[,5:26], rescale, to = c(-1, 1)))
data <- round(data, 3)
head(data)

# Regress out the effect of the site from the Discovery
data$site <- df_discovery$site
str(data)
data$site <- as.factor(data$site)
data <- datawizard::adjust(data, effect = "site", select = names(data)[1:22])
head(data)
data$site <- NULL
data <- round(data, 3)
str(data)

# Hierarchical k-means on the Discovery and centers on Replication

################### Discovery - K=3 ##################
str(data) # Discovery scaled and reg out

library(factoextra)
library(tidyr)
library(ggplot2)
library(ggpubr)

### K = 3 
set.seed(1234)
res.hk <- hkmeans(data, 3, hc.metric = "euclidean", hc.method = "ward.D")
data$cluster.hm3 <- as.factor(res.hk$cluster)
# visualize clusters
df_cluster <- data
str(df_cluster)
medie.km <- aggregate(. ~ cluster.hm3, data = df_cluster, FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -cluster.hm3, names_to = "BOLD_roi", values_to = "BOLD")
colori <- c(
  "yellow","gold","darkgoldenrod1", "darkorange", "coral", "chocolate2","chocolate","brown2","red","brown","darkred","darkorange4",
  
  "cyan2","cadetblue2","cornflowerblue","cyan4","blue4","chartreuse1","chartreuse4","darkgreen","magenta","blueviolet")
plot_barre.roi.cl <- ggplot(medie_long.km, aes(x = cluster.hm3, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Clusters", y = "Mean BOLD signal") +
  ylim(-0.8,0.8)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),  
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5) 
  )
plot_barre.roi.cl

# rename levels
data$cluster.hm3 <- factor(data$cluster.hm3, levels = c(1,2,3), labels = c("ST", "INT","GT"))

############### Replication - HKmeans #################
### projection on the Replication of the centroids by Discovery of only ST e GT
str(df) # Replication dataset
df$predicted_cluster <- apply(df, 1, function(row) {
  distances <- apply(res.hk$centers[c(1,3),], 1, function(centroid) sum((row - centroid)^2))
  return(which.min(distances))  
})
df$predicted_cluster <- as.factor(df$predicted_cluster)

# visualization 
df_cluster <- df
str(df_cluster)
medie.km <- aggregate(. ~ predicted_cluster, data = df_cluster, FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -predicted_cluster, names_to = "BOLD_roi", values_to = "BOLD")

plot_barre.roi.cl <- ggplot(medie_long.km, aes(x = predicted_cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Clusters", y = "Mean BOLD signal", title = "Distribuzione delle medie delle ROI per Cluster") +
  theme_pubclean() +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_barre.roi.cl

# rename levels
df$predicted_cluster <- factor(df$predicted_cluster, levels = c(1,2), labels = c("ST","GT"))

# add cluster variable to the dataset
df_ext_def <- df_ext
df_ext_def$pred.cluster.hm3.drop <- df$predicted_cluster

######## HIERARCHICAL K-MEANS ON THE REPLICATION WITHOUT PROJECTION ##########
# Scaled Replication dataset 
str(df)
# remove the predicted cluster
df$predicted_cluster <- NULL
### Compute hierarchical kmeans (K=2) ###
library(factoextra)
library(tidyr)
set.seed(1234)
res.hk <- hkmeans(df, 2, hc.metric = "euclidean", hc.method = "ward.D")

### VISUALIZE CLUSTERS ###
df.cluster <- df
df.cluster$cluster <- as.factor(res.hk$cluster)
str(df.cluster)
medie.km <- aggregate(. ~ cluster, data = df.cluster, FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -cluster, names_to = "BOLD_roi", values_to = "BOLD")
colori <- c(
  "yellow","gold","darkgoldenrod1", "darkorange", "coral", "chocolate2","chocolate","brown2","red","brown","darkred","darkorange4",
  
  "cyan2","cadetblue2","cornflowerblue","cyan4","blue4","chartreuse1","chartreuse4","darkgreen","magenta","blueviolet")

plot_barre.roi.cl <- ggplot(medie_long.km, aes(x = cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Clusters", y = "Mean BOLD signal", title = "Distribuzione delle medie delle ROI per Cluster") +
  theme_pubclean() +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_barre.roi.cl

# rename factor levels 
df.cluster$cluster <- factor(df.cluster$cluster, levels = c(1,2), labels = c("ST","GT"))
levels(df.cluster$cluster)
str(df.cluster)

# add cluster variable to the dataset
df_ext_def$true.cluster.hm2 <- df.cluster$cluster

# Contingency table
table(df_ext_def$pred.cluster.hm3.drop, df_ext_def$true.cluster.hm2)
library(caret)
cm <- confusionMatrix(df_ext_def$pred.cluster.hm3.drop, df_ext_def$true.cluster.hm2)
cm