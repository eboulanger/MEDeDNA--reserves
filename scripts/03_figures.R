#'##############################################################################
#'  
#'  codes to reproduce figures of the Boulanger et al. (under review) paper: "Environmental DNA metabarcoding reveals
#'  and unpacks a biodiversity conservation paradox in Mediterranean marine reserves." 
#'  
#' Script to run:
#' * figure 2: diversity indices (richness & dissimilarity)
#' * figure 3: partial Jaccard dbRDA
#' * figure 4: species scores x traits
#' * figure 5: richness GLM by trait category
#'   
#'  by Emilie Boulanger
#'  
#'##############################################################################

# set general plotting parameters

ggplot2::theme_set(theme_classic())

# purplescale for protection level
colres    <-  "#3f007d"
colout5   <-  "#6a51a3"
colout10  <-  "#9e9ac8"

# species categories
col.crypto   <- "#e31a1c" 
col.pelagic  <- "#1f78b4" 
col.benthic  <- "#33a12c" 
col.demersal <- "#b2e08a" 

#### figure 2. diversity indices ####

# violin plot raw richness
gtot <- ggplot(richness, aes(Protection, total)) +
  geom_violin(aes(fill = Protection), draw_quantiles = 0.5, alpha = 0.9) + # set boxplot colour by protection level
  #geom_boxplot(width = 0.1) +
  geom_jitter(aes(group = Protection), alpha = 0.5, width = 0.1, height = 0.1) +
  scale_fill_manual(values = c(colres, colout5, colout10)) + # set legend items order
  scale_x_discrete(labels = c("reserve", "5km outside", "10km outside")) +
  labs(y = "Overall species richness", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + # centralize title and hide legend
  coord_cartesian(ylim=c(0,50))

# conditional predicted richness from model
cond_total   <- cplot(mod_total   , "Protection", draw = F)
cond_total$xvals <- factor(cond_total$xvals, levels = c("reserve", "outside5", "outside10"))

cplot_tot <- ggplot(data=cond_total) +
  geom_point(aes(x= xvals, y=yvals, col = xvals), cex = 3) +
  geom_errorbar(aes(x=xvals, ymin=lower, ymax=upper, col=xvals), width = 0.5, size = 1) +
  scale_color_manual(values=c(colres, colout5,colout10)) +
  scale_x_discrete(labels = c("reserve", "5km outside", "10km outside")) +
  geom_text(aes(x=xvals, y=upper + 1), label = c("", "","**"), size = 6) +
  labs(x = "", y = "Predicted species richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color="none") +
  coord_cartesian(ylim=c(0,50))

# turnover vs nestedness
plot_beta <- function(df, comparison, title) {
    filter(df, Pairwise == comparison) %>% 
    ggplot(aes(x = Component, y = value)) +
    geom_violin(aes(fill = Component), alpha = 0.9, draw_quantiles = 0.5) +
    labs(title = title,
         x = "", y = "") +
    scale_fill_brewer(palette = "Dark2") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    coord_cartesian(ylim = c(0.00, 1))
}

g_res_out5   <- plot_beta(replicate_beta_all, "res_out5", "No-take reserve\nvs 5 km outside") + labs(y= "Beta diversity")
g_res_out10  <- plot_beta(replicate_beta_all, "res_out10", "No-take reserve\nvs 10 km outside")
g_out5_out10 <- plot_beta(replicate_beta_all, "out5_out10", "5 km outside\nvs 10 km outside")

g_res_out5  + g_res_out10 + g_out5_out10 + plot_layout(guides = "collect")

# combine all
fig2 <- (gtot + cplot_tot) / (g_res_out5  + g_res_out10 + g_out5_out10 + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 15))
fig2

# export figure
ggsave(plot = fig2, filename  = here::here("figures", "figure2.png"), width = 10,  height = 10, dpi = 600)

#### figure 3. partial dbRDA ####

# get scores
site_scores <- scores(part.dbrda)$sites     ## separating out the site scores
species_scores <- scores(part.dbrda)$species %>% data.frame()   ## separating out the species

# get most differentiated species along first axis
quant75 <- quantile(abs(species_scores$CAP1), probs = c(0.75))
species_scores_diff75 <- species_scores[which(abs(species_scores$CAP1) > quant75["75%"]),]
# add colour variable by vertical distribution
species_scores_diff75$col <- rep("other", nrow(species_scores_diff75))
species_scores_diff75$col[rownames(species_scores_diff75) %in% crypto_sp]   <- "crypto"
species_scores_diff75$col[rownames(species_scores_diff75) %in% pelagic_sp]  <- "pelagic"
species_scores_diff75$col[rownames(species_scores_diff75) %in% demersal_sp] <- "demersal"
species_scores_diff75$col[rownames(species_scores_diff75) %in% benthic_sp]  <- "benthic"
species_scores_diff75$col <- factor(species_scores_diff75$col, levels = c("crypto", "benthic", "demersal", "pelagic"))

# extract the percentage variability explained by axes
sumdbrda <- summary(part.dbrda)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(meta_env$Code), rownames(site_scores)) # verify that data in same order
site_scores_protection <- cbind(site_scores,select(meta_env, Protection))

# plot in ggplot
grda_sites <- ggplot(site_scores_protection, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = Protection,linetype = Protection,fill= Protection), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = species_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= species_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = Protection, fill = Protection), cex = 4, col = "black") +
  scale_fill_manual(values = c(colres, colout5, colout10),
                    name = "Protection", labels = c("reserve", "5km outside", "10km outside")) +
  scale_shape_manual(values = c(23:21),
                     name = "Protection", labels = c("reserve", "5km outside", "10km outside")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0, 1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_sites

grda_species <- ggplot() + 
  geom_segment(data= species_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc"))) + # all species
  geom_segment(data= species_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data= species_scores_diff75, 
                   aes(x= CAP1, y=CAP2, #hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2)),
                       col = col, fontface=3), #size = 3,
                   label = rownames(species_scores_diff75),
                   show.legend = F) +
  scale_color_manual(values = c(col.crypto, col.benthic, col.demersal, col.pelagic), name = "Vertical position", labels = c("Cryptobenthic", "Benthic", "Demersal", "Pelagic")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11)) +
  # re-add legend and change text legend key by making invisible points and overriding its key shape
  geom_point(data= species_scores_diff75, 
             aes(x=CAP1, y=CAP2, col=col),
             size = 0, stroke = 0) + 
  guides(colour = guide_legend(override.aes = list(size = 5, 
                                                   shape = c(utf8ToInt("C"), utf8ToInt("B"), utf8ToInt("D"), utf8ToInt("P")))))
grda_species

# combine

fig3 <-  (grda_sites + grda_species) + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face="bold", size = 15))
fig3

# export figure
ggsave(plot = fig3, filename = here::here("figures", "figure3.png"), width = 15,  height = 8, dpi = 600)

#### figure 4: species scores x traits ####

# prepare plotting data
sp_scores_traits <- cbind(sp_scores, sp_traits) %>% 
  mutate(Vertical_Distribution_adj = fct_relevel(Vertical_Distribution_adj, c("Cryptobenthic", "Benthic", "Demersal", "Pelagic")))
# prepare annotations
an_trophic    <- paste0("Kendall tau = ", round(cor.test(sp_scores$CAP1, sp_traits$Trophic_level,      method = "kendall")$estimate, 2),
                        "\nP = ",         round(cor.test(sp_scores$CAP1, sp_traits$Trophic_level,      method = "kendall")$p.value,  3))
an_length     <- paste0("Kendall tau = ", round(cor.test(sp_scores$CAP1, log(sp_traits$Common_length), method = "kendall")$estimate, 2),
                        "\nP = ",         round(cor.test(sp_scores$CAP1, log(sp_traits$Common_length), method = "kendall")$p.value,  3))
an_vulnerable <- paste0("Kendall tau = ", round(cor.test(sp_scores$CAP1, sp_traits$Vulnerability,      method = "kendall")$estimate, 2),
                        "\nP = ",         round(cor.test(sp_scores$CAP1, sp_traits$Vulnerability,      method = "kendall")$p.value,  3))
an_vertical   <- paste0("ANOVA P = ",     round(summary(aov(sp_scores$CAP1 ~ sp_traits$Vertical_Distribution_adj))[[1]]$`Pr(>F)`[1], 4))

grob_trophic    <- grobTree(textGrob(an_trophic   , x=0.03,  y=0.07, hjust=0, gp=gpar(fontsize=12)))
grob_length     <- grobTree(textGrob(an_length    , x=0.03,  y=0.07, hjust=0, gp=gpar(fontsize=12)))
grob_vulnerable <- grobTree(textGrob(an_vulnerable, x=0.03,  y=0.07, hjust=0, gp=gpar(fontsize=12)))
grob_vertical   <- grobTree(textGrob(an_vertical  , x=0.03,  y=0.07, hjust=0, gp=gpar(fontsize=12)))

# plots
gtrophic <- ggplot(sp_scores_traits) + 
  geom_point( aes(x= Trophic_level, y=CAP1), col = "black")  +
  geom_smooth(aes(x= Trophic_level, y=CAP1), col = "black", method = "lm") +
  labs(y = "dbRDA CAP1 species scores", x = "Trophic Level") +
  annotation_custom(grob_trophic)
  
glength <- ggplot(sp_scores_traits) +
  geom_point( aes(x = log(Common_length), y=CAP1), col = "black") +
  geom_smooth(aes(x = log(Common_length), y=CAP1), col = "black", method = "lm") +
  labs(y = "", x = "log(Common Length)")  +
  annotation_custom(grob_length)

gvulnerable <- ggplot(sp_scores_traits) +
  geom_point( aes(x = Vulnerability, y=CAP1), col = "black") +
  geom_smooth(aes(x = Vulnerability, y=CAP1), col = "black", method = "lm") +
  labs(y = "dbRDA CAP1 species scores")  +
  annotation_custom(grob_vulnerable)

gvertical <- ggplot(sp_scores_traits) + 
  geom_boxplot(aes(x=Vertical_Distribution_adj, y = CAP1, col = Vertical_Distribution_adj)) +
  geom_jitter(aes(x=Vertical_Distribution_adj, y = CAP1, col = Vertical_Distribution_adj), width = 0.1, height = 0, alpha = 01) +
  scale_color_manual(values = c(col.crypto, col.benthic, col.demersal, col.pelagic)) +
  theme(legend.position = "none") +
  labs(y = "", x = "Vertical distribution")  +
  annotation_custom(grob_vertical)
# add shapes representing post hoc results to vertical distribution boxplots
bp <- ggplot(sp_scores_traits, aes(x=Vertical_Distribution_adj, y=CAP1)) + geom_boxplot() 
shape_pos <- ggplot_build(bp)$data[[1]]$ymax
shape_df <- data.frame(x= levels(sp_scores_traits$Vertical_Distribution_adj), 
                        y = shape_pos + 0.05,
                        shape = c("up","up","up", "down")) # post hoc grouping

gvertical_shp <- gvertical + geom_point(data=shape_df, aes(x=x, y=y, shape = shape), cex = 2, fill = "black") +
  scale_shape_manual(values = c(24,25))

fig4 <- (gtrophic + glength) / (gvulnerable + gvertical_shp) + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face="bold", size = 15))
fig4

# export figure
ggsave(plot = fig4, filename = here::here("figures", "figure4.png"), width = 10,  height = 8, dpi = 600)

#### figure 5: category GLMs ####
# get conditional predicted richness from model
cond_crypto  <- cplot(mod_crypto    , "Protection", draw = F)  %>% mutate(index = rep("crypto",3))
cond_pelagic <- cplot(mod_pelagic   , "Protection", draw = F)  %>% mutate(index = rep("pelagic",3))
cond_rare    <- cplot(mod_rare      , "Protection", draw = F)  %>% mutate(index = rep("rare",3))
cond_vuln    <- cplot(mod_vulnerable, "Protection", draw = F)  %>% mutate(index = rep("vulnerable",3))
cond_res <- rbind(cond_crypto, cond_pelagic, cond_rare, cond_vuln) %>% 
  mutate(xvals = factor(xvals, levels = c("reserve", "outside5", "outside10")))

# create plotting function
ggcplot <- function(df, category) {
  df %>% filter(index == category) %>% 
  ggplot() +
    geom_point(aes(x= xvals, y=yvals, col = xvals), cex = 3) +
    geom_errorbar(aes(x=xvals, ymin=lower, ymax=upper, col=xvals), width = 0.5, size = 1) +
    scale_color_manual(values=c(colres, colout5,colout10)) +
    scale_x_discrete(labels = c("reserve", "5km outside", "10km outside")) +
    geom_text(aes(x=xvals, y=upper + 1), label = c("", "","**"), size = 6) +
    labs(x = "", y = "") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color="none")
}

# plot for all categories
cplot_crypto     <- ggcplot(cond_res, "crypto")     + labs(title = "Cryptobenthic", y="Predicted species richness")
cplot_pelagic    <- ggcplot(cond_res, "pelagic")    + labs(title = "Pelagic")
cplot_rare       <- ggcplot(cond_res, "rare")       + labs(title = "Rare")
cplot_vulnerable <- ggcplot(cond_res, "vulnerable") + labs(title = "Vulnerable")

fig5 <- (cplot_crypto| cplot_pelagic | cplot_rare | cplot_vulnerable) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face="bold", size = 15))
fig5

# export figure
ggsave(plot = fig5, filename = here::here("figures", "figure5.png"), width = 15,  height = 7, dpi = 600)

