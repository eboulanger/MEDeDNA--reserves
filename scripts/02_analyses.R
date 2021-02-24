#'##############################################################################
#'  
#'  codes to reproduce analyses and figures of the Boulanger et al. (under review) paper: "Environmental DNA metabarcoding reveals
#'  and unpacks a biodiversity conservation paradox in Mediterranean marine reserves." 
#'  
#' Script to run:
#' * richness analyses 
#' * dissimilarity analyses
#' * species trait analyses
#'   
#'  by Emilie Boulanger
#'  
#'##############################################################################

#### (a) Species richness paradox ####

# calculate total species richness and add metadata
richness <- species %>% 
  rowSums(.) %>%       # richness per samle
  as.data.frame() %>% 
  `colnames<-` ("total") %>% 
  mutate(Replicate = rownames(.)) %>% 
  left_join(select(meta_env,  Replicate, Site, Region, Protection, envPC1, envPC2, envPC3, envPC4), by = "Replicate") %>% 
  mutate(Protection = fct_relevel(Protection, c("reserve", "outside5", "outside10"))) %>%   # set the correct protection level order
  mutate(Region = fct_inorder(Region)) # reorder factor levels in same order as regions appear in dataset

# compare richness by protection
boxplot(richness$total ~ richness$Protection)

summary.stats <- function(x) {
  c(mean = mean(x), sd = sd(x))
}
richness %>% filter(Protection == "reserve")   %>% pull(total) %>% summary.stats()
richness %>% filter(Protection == "outside5")  %>% pull(total) %>% summary.stats()
richness %>% filter(Protection == "outside10") %>% pull(total) %>% summary.stats()

kruskal.test(richness$total ~ richness$Protection)

# model richness by protection, taking environmental variables (PCA axes) into account

mod_total  <- lm(total ~ Protection * (envPC1 + envPC2 + envPC3 + envPC4), data = richness)
rsq(mod_total)
# overall p-value
null_mod_total <- lm(total ~ 1, data = richness)
anova(mod_total, null_mod_total, test="Chisq")
# average marginal effects
margins(mod_total) %>% summary()
# test effect region on model residuals
kruskal.test(residuals(mod_total) ~ richness$Region)

#### (b) Species dissimilarity between assemblages ####

# calculate jaccard dissimilarity, turnover and nestedness components
totalBETA <- beta.pair(species, index.family = "jaccard")
turnover <- as.matrix(totalBETA$beta.jtu)
diag(turnover) <- NA
nestedness <- as.matrix(totalBETA$beta.jne)
diag(nestedness) <- NA
totaldis <- as.matrix(totalBETA$beta.jac)
diag(totaldis) <- NA

# compare between protection levels
diss.prot.levels <- function(dissimilarity.df, samples, level.a, level.b) {
  dissimilarity.df %>% melt() %>% 
    filter(X1 %in% vars_select(samples, contains(paste0(level.a,"_")))) %>% 
    filter(X2 %in% vars_select(samples, contains(paste0(level.b,"_")))) %>% 
    select(value)
} # function to extract pairwise dissimilarities 

res_out5_tot  <- diss.prot.levels(totaldis,   rownames(species), level.a = "1", level.b = "2") # reserves vs 5 km outside
res_out5_turn <- diss.prot.levels(turnover,   rownames(species), level.a = "1", level.b = "2")
res_out5_nest <- diss.prot.levels(nestedness, rownames(species), level.a = "1", level.b = "2")
res_out5_all  <- cbind(res_out5_tot, res_out5_turn, res_out5_nest) %>% 
  `colnames<-`(c("Dissimilarity", "Turnover", "Nestedness")) %>% 
  mutate(Pairwise = "res_out5")

res_out10_tot  <- diss.prot.levels(totaldis,   rownames(species), level.a = "1", level.b = "3") # reserves vs 10 km outside
res_out10_turn <- diss.prot.levels(turnover,   rownames(species), level.a = "1", level.b = "3")
res_out10_nest <- diss.prot.levels(nestedness, rownames(species), level.a = "1", level.b = "3")
res_out10_all  <- cbind(res_out10_tot, res_out10_turn, res_out10_nest) %>% 
  `colnames<-`(c("Dissimilarity", "Turnover", "Nestedness")) %>% 
  mutate(Pairwise = "res_out10")

out5_out10_tot  <- diss.prot.levels(totaldis,   rownames(species), level.a = "2", level.b = "3") # 5 km outside vs 10 km outside
out5_out10_turn <- diss.prot.levels(turnover,   rownames(species), level.a = "2", level.b = "3")
out5_out10_nest <- diss.prot.levels(nestedness, rownames(species), level.a = "2", level.b = "3")
out5_out10_all  <- cbind(out5_out10_tot, out5_out10_turn, out5_out10_nest) %>% 
  `colnames<-`(c("Dissimilarity", "Turnover", "Nestedness")) %>% 
  mutate(Pairwise = "out5_out10")

replicate_beta_all <- rbind(res_out5_all, res_out10_all, out5_out10_all) %>% 
  melt(id.vars = "Pairwise", variable_name = "Component") 

# extract summary statistics
# - mean pairwise dissimilarity between protection levels per comparison
replicate_beta_all %>% filter(Pairwise == "res_out5"   & Component == "Dissimilarity") %>% pull(value) %>% summary.stats()
replicate_beta_all %>% filter(Pairwise == "res_out10"  & Component == "Dissimilarity") %>% pull(value) %>% summary.stats()
replicate_beta_all %>% filter(Pairwise == "out5_out10" & Component == "Dissimilarity") %>% pull(value) %>% summary.stats()
# - mean pairwise turnover for all comparisons
replicate_beta_all %>% filter(Component == "turnover") %>% pull(value) %>% summary.stats() * 100
# - proportion of dissimilarity that is turnover
prop.turn <- c(  res_out5_all$Turnover /   res_out5_all$Dissimilarity * 100,
                res_out10_all$Turnover /  res_out10_all$Dissimilarity * 100,
               out5_out10_all$Turnover / out5_out10_all$Dissimilarity * 100)
summary.stats(prop.turn)
# - proportion of dissimilarity that is nestedness
prop.nest <- c(  res_out5_all$Nestedness /   res_out5_all$Dissimilarity * 100,
                res_out10_all$Nestedness /  res_out10_all$Dissimilarity * 100,
               out5_out10_all$Nestedness / out5_out10_all$Dissimilarity * 100)
summary.stats(prop.nest)

# distance-based redundancy analysis
# with jaccard dissimilarity
tot.dbrda  <- capscale(species ~ Protection + envPC1 + envPC2 + envPC3 + envPC4, data = meta_env, distance = "jaccard")
RsquareAdj(tot.dbrda)
anova(tot.dbrda, by = "axis",   permutations = 9999)
anova(tot.dbrda, by = "margin", permutations = 9999)

# on turnover between samples
turn.dbrda      <- capscale(turnover   ~ Protection + envPC1 + envPC2 + envPC3 + envPC4, data = meta_env)
anova(turn.dbrda, permutations = 9999)
anova(turn.dbrda, by = "margin", permutations = 9999)
# on nestedness between samples
nest.dbrda      <- capscale(nestedness ~ Protection + envPC1 + envPC2 + envPC3 + envPC4, data = meta_env)
anova(nest.dbrda, permutations = 9999)

# partial distance-based redundancy analysis
part.dbrda <- capscale(species ~ Protection + Condition(envPC1 + envPC2 + envPC3 + envPC4), data = meta_env, distance = "jaccard")
RsquareAdj(part.dbrda)
anova(part.dbrda, permutations = 9999)

#### (c) Unpacking the paradox by traits ####

# species scores x traits
# extract species scores from the partial dbrda
sp_scores <- scores(part.dbrda)$species %>% as.data.frame()
# check if same order as species trait data
identical(rownames(sp_scores), rownames(sp_traits))

# test relationships & difference in means for different traits
cor.test(sp_scores$CAP1, sp_traits$Trophic_level, method = "kendall")
cor.test(sp_scores$CAP1, log(sp_traits$Common_length), method = "kendall")
cor.test(sp_scores$CAP1, sp_traits$Vulnerability, method = "kendall")

anova_vertical <- aov(sp_scores$CAP1 ~ sp_traits$Vertical_Distribution_adj)
summary(anova_vertical)
TukeyHSD(anova_vertical)

# calculate richness by species categories 
# cryptobenthic species richness
rich_crypto <- species[, crypto_sp] %>% 
  rowSums(.) %>%       # richness per samle
  as.data.frame() %>% `colnames<-` ("crypto")

# pelagic species richness
rich_pelagic <- species[, pelagic_sp] %>% 
  rowSums(.) %>% 
  as.data.frame() %>% `colnames<-` ("pelagic")

# rare species richness
rich_rare <- species[, rare_sp] %>% 
  rowSums(.) %>% 
  as.data.frame() %>% `colnames<-` ("rare")

# vulnerable species richness
rich_vulnerable <- species[,vulnerable_sp] %>% 
  rowSums(.) %>% 
  as.data.frame() %>% `colnames<-` ("vulnerable")

# verify if all in same order and add to one df
stopifnot(identical(Reduce(union, list(richness$Replicate,row.names(rich_crypto),
                                       row.names(rich_pelagic), row.names(rich_rare))),
                           row.names(rich_vulnerable)) == TRUE)
richness2 <- cbind(richness, rich_crypto, rich_pelagic, rich_rare, rich_vulnerable) %>% 
  select(Replicate, total, crypto, pelagic, rare, vulnerable, everything())

# model different richness values by protection, taking environmental variables (PCA axes) into account

mod_crypto      <- glm(crypto     ~ Protection * (envPC1 + envPC2 + envPC3 + envPC4), data = richness2, family = "gaussian")
mod_pelagic     <- glm(pelagic    ~ Protection * (envPC1 + envPC2 + envPC3 + envPC4), data = richness2, family = "gaussian")
mod_rare        <- glm(rare       ~ Protection * (envPC1 + envPC2 + envPC3 + envPC4), data = richness2, family = "gaussian")
mod_vulnerable  <- glm(vulnerable ~ Protection * (envPC1 + envPC2 + envPC3 + envPC4), data = richness2, family = "poisson")

rsq(mod_crypto    )
rsq(mod_pelagic   )
rsq(mod_rare      )
rsq(mod_vulnerable)

# overall p-value
null_mod_crypto     <- glm(crypto     ~ 1, data = richness2, family = "gaussian")
null_mod_pelagic    <- glm(pelagic    ~ 1, data = richness2, family = "gaussian")
null_mod_rare       <- glm(rare       ~ 1, data = richness2, family = "gaussian")
null_mod_vulnerable <- glm(vulnerable ~ 1, data = richness2, family = "poisson")

anova(mod_crypto    , null_mod_crypto    , test="Chisq")
anova(mod_pelagic   , null_mod_pelagic   , test="Chisq")
anova(mod_rare      , null_mod_rare      , test="Chisq")
anova(mod_vulnerable, null_mod_vulnerable, test="Chisq")

# average marginal effects
margins(mod_crypto    ) %>% summary()
margins(mod_pelagic   ) %>% summary()
margins(mod_rare      ) %>% summary()
margins(mod_vulnerable) %>% summary()

# test effect region on model residuals
kruskal.test(residuals(mod_crypto    ) ~ richness$Region)
kruskal.test(residuals(mod_pelagic   ) ~ richness$Region)
kruskal.test(residuals(mod_rare      ) ~ richness$Region)
kruskal.test(residuals(mod_vulnerable) ~ richness$Region)
