
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method = "spearman")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8
  text(0.5, 0.5, txt, cex = cex.cor)
}

# Function to perform cross-validation with glmnet and calculate log loss
cv_log_loss <- function(X, y, nfolds = 10, alpha = 1, lambda = NULL) {
  
  folds <- sample(rep(seq_len(nfolds), length.out = nrow(X)))
  log_loss_values <- numeric(nfolds)
  accuracy_values <- numeric(nfolds)
  roc_curves <- list()
  auroc_values <- numeric(nfolds)
  
  for (i in 1:nfolds) {
    
    test_idx <- which(folds == i)
    train_idx <- which(folds != i)
    fit_glmnet <- glmnet::glmnet(x = X[train_idx, ], y = y[train_idx], family = "binomial", alpha = alpha, lambda = lambda, intercept = FALSE, standardize = TRUE)
    predicted_probabilities <- as.vector(predict(fit_glmnet, newx = X[test_idx, ], type = "response"))
    
    true_labels <- y[test_idx]
    log_loss_values[i] <- log_loss(true_labels, predicted_probabilities)
    
    # Calculate accuracy
    predicted_labels <- ifelse(predicted_probabilities > 0.5, 1, 0)  # Thresholding at 0.5
    accuracy_values[i] <- mean(predicted_labels == true_labels)
    
    # Compute ROC curve
    roc_curves[[i]] <- pROC::roc(true_labels, predicted_probabilities)
    
    # Compute AUROC
    auroc_values[i] <- auc(roc_curves[[i]])
  }
  return(list(log_loss = log_loss_values, accuracy = accuracy_values, auroc = auroc_values, roc_curves = roc_curves))
}

log_loss <- function(y_true, y_pred) {
  return(-mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred)))
}

plot_roc_curves <- function(roc_curves) {
  roc_df <- data.frame(fpr = c(), tpr = c(), fold = c())
  
  for (i in seq_along(roc_curves)) {
    roc_df_add <- data.frame(fpr = c(roc_df$fpr, (1 - roc_curves[[i]]$specificities)),
                             tpr = c(roc_df$tpr, roc_curves[[i]]$sensitivities),
                             fold = c(roc_df$fold, rep(i, length(roc_curves[[i]]$sensitivities))))
    if(i == 1){
      roc_df <- roc_df_add
    }else{
      roc_df <- rbind(roc_df, roc_df_add)
    }
  }
  
  roc_plt <- roc_df %>% arrange(fold, tpr) %>%
    ggplot2::ggplot(ggplot2::aes(x = fpr, y = tpr, col = as.factor(fold))) +
    ggplot2::geom_line(alpha = 1, size = 1) +
    ggplot2::theme_bw() +
    ggsci::scale_color_d3() +
    ggplot2::geom_abline(slope = 1, linetype = "dashed", size = 1) +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.text = ggplot2::element_text(colour = "black", size = 14),
                   strip.background = ggplot2::element_rect(fill = "white"),
                   axis.title = ggplot2::element_text(face = "bold", size = 16),
                   strip.text = ggplot2::element_text(face = "bold", size = 14)) +
    ggplot2::geom_text(x = 0.7, y = 0.05, label = paste0("AUC: ", round(mean(fit$auroc), 2), " ± ", round(sd(fit$auroc), 2)), col = "black", size = 5, hjust = 0) +
    ggplot2::geom_text(x = 0.7, y = 0.0, label = paste0("Log-Loss: ", round(mean(fit$log_loss), 2), " ± ", round(sd(fit$log_loss), 2)), col = "black", size = 5, hjust = 0)
  
  return(list(roc_df = roc_df, roc_plot = roc_plt))
}

plot_density <- function(df_ana, feature, logfeature){
  
  igg_df <- df_ana  %>%
    mutate(!!logfeature := log10(!!sym(feature) + 1)) %>%
    dplyr::select(!!sym(feature), !!sym(logfeature)) %>%
    tidyr::pivot_longer(cols = matches(paste0(feature, "|", logfeature)), names_to = "transform", values_to = "values")
  theo_params <- igg_df %>% group_by(transform) %>% summarise(mu_est = mean(values),
                                                              sd_est = sd(values),
                                                              xmin = min(values),
                                                              xmax = max(values))
  dens_plts <- list()
  igg_names <- unique(igg_df$transform)
  
  mains <- c("Density raw data", "Density of logarithmized")
  xlabs <- c(feature, logfeature)
  
  m <- rbind(c(0, 0.5, 0, 1),
             c(0.5, 1, 0, 1))
  
  split.screen(m)
  
  for(i in seq_along(igg_names)){
    
    screen(i)
    
    if(i == 1){
      par(mar = c(4.1, 4.1, 2.1, .5))
    }
    
    if(i == 2){
      par(mar = c(4.1, 4.1, 2.1, .5))
    }
    
    values <- igg_df %>% filter(transform %in% igg_names[i]) %>% pull(values)
    x_vals <- seq(min(values), max(values), length.out = 1000)
    
    hist(values, xlab = xlabs[i], prob = T, col = "white", main = mains[i], breaks = 30)
    lines(x = x_vals, y = dnorm(x = x_vals, mean = mean(values), sd = sd(values)), lty = 1, lwd = 2)
    box()
  }
  
  close.screen(all.screens = TRUE)
  
}

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL,
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.02,0.96),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02),
                     bottomcenter = c(0.5525, 0.02),
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}
