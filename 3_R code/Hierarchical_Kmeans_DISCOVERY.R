# MANUSCRIPT: ST-GT 

# ANALYSIS 1: HIERARCHICAL K-MEANS CLUSTERING - DISCOVERY COHORT 

# Import dataset 
setwd()
data <- as.data.frame(readxl::read_excel("dataset_clustering.xlsx"))

# make "site" variable as factor
data$site <- as.factor(data$site) 
str(data)

# Hierarchical K-means - Full dataset

### STEP 1: SCALING -1 1 ###
library(scales)
df <- as.data.frame(lapply(data[,5:26], rescale, to = c(-1, 1))) # select only continuous variable, i.e., BOLD signal values (22 columns)
head(df)
df <- round(df, 3)

### STEP 2: ADJUST (regress out site) ###
df$site <- data$site # add site to the df
str(df)
df <- datawizard::adjust(df, effect = "site", select = names(data)[1:22])
head(df)
df$site <- NULL
df <- round(df, 3)

### STEP 3: compute hierarchical kmeans clustering ###
library(factoextra)
library(tidyr)
str(df) # make sure that are 22 columns, only BOLD signal in the 22 ROIs
set.seed(1234)
res.hk <- hkmeans(df, 3, hc.metric = "euclidean", hc.method = "ward.D")
names(res.hk)
res.hk
# Visualize the tree
fviz_dend(res.hk, cex = 0.6, palette = "jco", 
          rect = TRUE, rect_border = "jco", rect_fill = TRUE)
# Visualize the hkmeans final clusters
fviz_cluster(res.hk, palette = "uchicago", repel = FALSE,
             ggtheme = theme_classic())
# customize: visualize plot on the two PCs
fviz_cluster(res.hk, data = df, 
             ellipse = TRUE,
             ellipse.type = "norm",
             ellipse.alpha = 0,
             show.clust.cent = TRUE,
             stand = TRUE,
             geom = c("point"),
             palette = c("#1B9E77", "#D95F02", "#7570B3"),
             ggtheme = theme_minimal(),
             main = "Hierarchical K-means cluster plot")

### STEP 4: VISUALIZE CLUSTERS ###
library(ggpubr)
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

# rename factor levels based on the visualization
df.cluster$cluster <- factor(df.cluster$cluster, levels = c(1,2,3), labels = c("ST","INT","GT"))
levels(df.cluster$cluster)
str(df.cluster)
table(df.cluster$cluster) 
