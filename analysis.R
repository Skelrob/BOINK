library(tidyverse)

# load population and individual tables
p <- read_csv("pop_info_9-Feb-2019.csv") %>%
      rename(popnumber=pop) %>%
      rename_at(vars(contains("1951")), 
                funs(sub("1951_1980", "", .))) %>%
      filter(!is.na(newpopnumber),
             newpopnumber != "H") %>%
      select(popnumber, longitude, latitude, cwd:aet) %>%
      mutate(popnumber=as.integer(popnumber))
i <- read_tsv("Finalrevised2analysis_Qd2018.txt") %>%
      rename(embolism="percentembolism_%") %>%
      select(tag:embolism) %>%
      mutate(individual=factor(individual)) %>%
      distinct()

# view all response curves
ggplot(i, aes(lwp_MPa, embolism, color=individual, order=lwp_MPa)) +
      geom_line() +
      facet_grid(popnumber ~ site + tissue)

# investigate the weirdo
i %>%
      filter(individual==1,
             popnumber==18,
             site=="Wild",
             tissue=="stem") %>%
      distinct() %>%
      ggplot(aes(lwp_MPa, embolism, color=individual, order=lwp_MPa)) +
      geom_point() +
      facet_grid(popnumber ~ site + tissue)
# seems to have two curves. huh.

# isolate P50, ghetto-style
# and join to population data
p50 <- i %>%
      group_by(individual, popnumber,
               site, tissue) %>%
      #filter(embolism >= 50) %>%
      #arrange(embolism) %>%
      filter(embolism <= 50) %>%
      arrange(desc(embolism)) %>%
      slice(1) %>%
      ungroup() %>%
      rename(p50 = lwp_MPa) %>%
      left_join(p) %>%
      mutate(popnumber=factor(popnumber))

# some exploratory plots
p %>%
      ggplot(aes(longitude, latitude, color=cwd)) +
      geom_point() +
      scale_color_gradientn(colours=c("blue", "green", "yellow", "red"))

p50 %>%
      select(-embolism) %>%
      spread(tissue, p50) %>%
      ggplot(aes(leaf, stem, color=popnumber, 
                 shape=site, linetype=site)) +
      geom_abline(slope=1, intercept=0) +
      geom_line() +
      geom_point() +
      coord_fixed()

p50 %>%
      select(-embolism) %>%
      group_by(site, tissue, popnumber) %>% 
      summarize(p50=mean(p50), cwd=mean(cwd)) %>%
      spread(site, p50) %>%
      ggplot(aes(Hopland, Wild, color=cwd, group=popnumber)) +
      geom_abline(slope=1, intercept=0) +
      geom_line() +
      geom_point(aes(shape=tissue)) +
      coord_fixed() +
      scale_color_viridis_c()


p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), cwd=mean(cwd)) %>%
      ggplot(aes(site, p50, group=popnumber, color=cwd)) +
      geom_line() +
      facet_grid(.~tissue)

p50 %>%
      ggplot(aes(tissue, p50, color=site)) +
      geom_point() +
      facet_grid(.~popnumber)

p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50_se = sd(p50)/sqrt(n()),
                p50 = mean(p50), 
                cwd=mean(cwd)) %>%
      ggplot(aes(cwd, p50, ymin=p50-p50_se, ymax=p50+p50_se,
                 color=tissue, linetype=site)) +
      geom_line() +
      geom_errorbar(alpha=.5) +
      geom_point()

p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), 
                cwd=mean(cwd), djf=mean(djf)) %>%
      ggplot(aes(cwd, djf, color=p50)) +
      geom_point(size=5) +
      facet_grid(tissue ~ site) +
      scale_color_gradientn(colours=c("magenta", "gray20", "green"))

p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), 
                latitude=mean(latitude),
                longitude=mean(longitude)) %>%
      ggplot(aes(latitude, longitude, color=p50)) +
      geom_point(size=5) +
      facet_grid(tissue ~ site) +
      scale_color_gradientn(colours=c("magenta", "gray20", "green"))


# we'd expect stronger djf response in stems than leaves
# i.e. interaction between tissue and temperature
p50 %>% 
      lmer(p50 ~ djf + jja + tissue + site + (1 | popnumber/individual), data=.) %>%
      summary()

# we'd expect stronger response to tmin in stems
# and to tmax/cwd in leaves




p50 %>%
      filter(tissue=="stem") %>%
      lmer(p50 ~ cwd + djf + site + (1 | popnumber/individual), data=.) %>%
      summary()


p50 %>%
      filter(site=="Wild") %>%
      lmer(p50 ~ djf + cwd + tissue + (1 | popnumber/individual), data=.) %>%
      summary()

p50 %>%
      lmer(p50 ~ djf + cwd + tissue + (1 | popnumber/individual), data=.) %>%
      summary()

p50 %>%
      filter() %>%
      lmer(p50 ~ djf + jja + tissue + site + (1 | popnumber/individual), data=.) %>%
      summary()


# what if we look just at the most (or least) resilient individual per popn?
p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize_at(vars(p50:aet), min) %>%
      lm(p50 ~ jja + djf + tissue + site, data=.) %>%
      summary()


# attempt to control for distance using pairwise differences
pw_p50_diff ~ pw_cwd_diff + geo_dist

ps <- p50 %>% filter(site=="Wild")
diffs <- data.frame(p50 = outer(ps$p50, ps$p50, '-') %>% as.vector(),
             cwd = outer(ps$cwd, ps$cwd, '-') %>% as.vector(),
             djf = outer(ps$djf, ps$djf, '-') %>% as.vector(),
             dist = dist(select(ps, longitude, latitude)) %>% 
                   as.matrix() %>% as.vector())

lm(p50 ~ cwd + djf + dist, data=diffs %>% filter(cwd > 0)) %>%
      summary()

ggplot(diffs %>% filter(cwd>0), aes(cwd, p50)) + geom_point() + geom_smooth()




# smoothing neighborhood, using herbarium popns
