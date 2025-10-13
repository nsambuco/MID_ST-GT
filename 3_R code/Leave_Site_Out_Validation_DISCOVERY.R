# ANALYSIS 2: CLUSTERING VALIDATION IN THE DISCOVERY COHORT
################## LEAVE-SITE-OUT VALIDATION ####################
############## SENSITIVY, SPECIFICITY, ACCURACY #################
#################################################################
# import raw dataset
setwd()
data <- as.data.frame(readxl::read_excel("dataset_clustering.xlsx"))
data$site <- as.factor(data$site)
str(data)
### validation on one site ###
# Split the dataset into training and testing set  
train_data <- data[-which(data$site=="PARIS"), c(2, 5:26)]
test_data <- data[which(data$site=="PARIS"), 5:26]
str(train_data)
str(test_data)
test_site <- "PARIS"  # change each time for each site 

# STEP 1: SCALING -1 1 
# training 
library(scales)
library(dplyr)
train_site <- as.factor(train_data$site)
# drop level site
table(train_data$site) 
train_data$site <- droplevels(train_data$site)
levels(train_data$site)
table(train_data$site)
train_site <- as.factor(train_data$site)
# scaling -1 1
train_data_scaled <- as.data.frame(lapply(train_data[,-1], rescale, to = c(-1, 1)))
str(train_data_scaled)
train_data_scaled <- train_data_scaled %>%
  bind_cols(train_site)
colnames(train_data_scaled)[23] <- "site"
train_data_scaled <- train_data_scaled %>%
  select(site, everything())
# testing set
test_data_scaled <- as.data.frame(lapply(test_data, rescale, to = c(-1,1)))
str(test_data_scaled)

# STEP 2: ADJUST (regress out site)
train_data_scaled_adjust <- datawizard::adjust(train_data_scaled, effect = "site", select = names(train_data_scaled)[2:23])
levels(train_data_scaled_adjust$site)

# STEP 3: Hierarchical K-means training set
set.seed(1234)
hkmeans_train <- factoextra::hkmeans(train_data_scaled_adjust[,2:23], 3, hc.metric = "euclidean", hc.method = "ward.D")
train_data_scaled_adjust$cluster <- as.factor(hkmeans_train$cluster)
levels(train_data_scaled_adjust$cluster)

# STEP 3.1: Visualize the clusters
df_cluster <- train_data_scaled_adjust[,2:24]
str(df_cluster)
# train
medie.km <- aggregate(. ~ cluster, data = df_cluster, FUN = mean)
library(tidyr)
library(ggplot2)
library(ggpubr)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -cluster, names_to = "BOLD_roi", values_to = "BOLD")
colori <- c(
  "yellow","gold","darkgoldenrod1", "darkorange", "coral", "chocolate2","chocolate","brown2","red","brown","darkred","darkorange4",
  
  "cyan2","cadetblue2","cornflowerblue","cyan4","blue4","chartreuse1","chartreuse4","darkgreen","magenta","blueviolet")

plot_barre.cl.train <- ggplot(medie_long.km, aes(x = cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Train Clusters", y = "Mean BOLD signal", title = paste0("Training clusters: ", test_site)) +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
plot_barre.cl.train

# STEP 4: Prediction on the testing set using the training centroids (ST and GT)
str(test_data_scaled)
test_data_scaled$predicted_cluster <- apply(test_data_scaled, 1, function(row) {
  distances <- apply(hkmeans_train$centers[c(1,3),], 1, function(centroid) sum((row - centroid)^2))
  return(which.min(distances))  
})
test_data_scaled$predicted_cluster <- as.factor(test_data_scaled$predicted_cluster)

# STEP 4.1: Visualize the predicted cluster on the test
df_cluster <- test_data_scaled
# predicted
medie.km <- aggregate(. ~ predicted_cluster, data = df_cluster, FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -predicted_cluster, names_to = "BOLD_roi", values_to = "BOLD")
plot_barre.cl.pred <- ggplot(medie_long.km, aes(x = predicted_cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Predicted Clusters", y = "Mean BOLD signal", title = paste0("Projection Test Clusters: ", test_site)) +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
plot_barre.cl.pred

# STEP 5: hierarchical K-means on the testing set (k=2)
set.seed(1234)
hkmeans_test_clean <- factoextra::hkmeans(test_data_scaled[,-23], 2, hc.metric = "euclidean", hc.method = "ward.D")
test_data_scaled$true_cluster <- hkmeans_test_clean$cluster
test_data_scaled$true_cluster <- as.factor(test_data_scaled$true_cluster)

# STEP 5.1: Visualize true cluster 
df_cluster <- test_data_scaled
# true
medie.km <- aggregate(. ~ true_cluster, data = df_cluster[,-23], FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -true_cluster, names_to = "BOLD_roi", values_to = "BOLD")
plot_barre.cl.true <- ggplot(medie_long.km, aes(x = true_cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "True Clusters", y = "Mean BOLD signal", title = paste0("Test True Clusters: ", test_site)) +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
plot_barre.cl.true

# STEP 6.1: Compute centroids of predicted clusters and true clusters
predicted_centroids <- aggregate(. ~ predicted_cluster, data = test_data_scaled[, -24], FUN = mean)
true_centroids <- aggregate(. ~ true_cluster, data = test_data_scaled[, -23], FUN = mean)
predicted_centroids <- predicted_centroids[,-1]
true_centroids <- true_centroids[,-1]

# STEP 6.2: Compute distance between predicted centroids and true ones
distance_matrix <- as.matrix(dist(rbind(predicted_centroids, true_centroids)))

distance_pred_true <- distance_matrix[1:nrow(predicted_centroids), (nrow(predicted_centroids)+1):nrow(distance_matrix)]

# STEP 6.3: Find optimal correspondance between predicted and true
matching <- clue::solve_LSAP(distance_pred_true, maximum = FALSE)

# STEP 6.4: Reassign levels of clusters based on the matching
test_data_scaled$predicted_cluster <- factor(test_data_scaled$predicted_cluster,
                                             levels = 1:length(matching),
                                             labels = as.character(matching))

# STEP 7: Contingency table updated with alligned clusters
contingency_table <- table(test_data_scaled$predicted_cluster, test_data_scaled$true_cluster)
print("Confusion Matrix:")
print(contingency_table)

# STEP 8: Metrics
test_data_scaled$predicted_cluster <- factor(test_data_scaled$predicted_cluster, levels = levels(test_data_scaled$true_cluster))
cm <- caret::confusionMatrix(test_data_scaled$predicted_cluster, test_data_scaled$true_cluster)
print("Confusion Matrix Metrics:")
print(cm)

# STEP 8.1: Save metrics
results <- data.frame(
  Site = character(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Balanced.Accuracy = numeric(),
  Overall.Accuracy = numeric(),
  Accuracy.pvalue = numeric(),
  Kappa = numeric(),
  stringsAsFactors = FALSE
)

sensitivity <- round(cm$byClass["Sensitivity"], digits = 3)
specificity <- round(cm$byClass["Specificity"], digits = 3)
balancedaccuracy <- round(cm$byClass["Balanced Accuracy"], digits = 3)
accuracy <- round(cm$overall["Accuracy"], digits = 3)
pvalue <- cm$overall["AccuracyPValue"]
kappa <- round(cm$overall["Kappa"], digits = 3)

results <- rbind(results, data.frame(
  Site = test_site,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Balanced.Accuracy = balancedaccuracy,
  Overall.Accuracy = accuracy,
  Accuracy.pvalue = pvalue,
  Kappa = kappa
))

### STEP 9: visualize clusters
df_cluster <- test_data_scaled
str(df_cluster)
# predicted
medie.km <- aggregate(. ~ predicted_cluster, data = df_cluster[,-24], FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -predicted_cluster, names_to = "BOLD_roi", values_to = "BOLD")
plot_barre.cl.test <- ggplot(medie_long.km, aes(x = predicted_cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "Predicted Clusters", y = "Mean BOLD signal", title = paste0("Realigned Projection Clusters: ", test_site)) +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
plot_barre.cl.test
# true
medie.km <- aggregate(. ~ true_cluster, data = df_cluster[,-23], FUN = mean)
medie_long.km <- medie.km %>%
  pivot_longer(cols = -true_cluster, names_to = "BOLD_roi", values_to = "BOLD")
plot_barre.cl.test.2 <- ggplot(medie_long.km, aes(x = true_cluster, y = BOLD, fill = BOLD_roi)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.5) +
  labs(x = "True Clusters", y = "Mean BOLD signal", title = paste0("Realigned True Clusters: ",test_site)) +
  ylim(-1,1)+
  scale_fill_manual(values = colori)+
  theme_bw() +
  theme_pubclean(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
plot_barre.cl.test.2

# repeat for each site starting from STEP 1