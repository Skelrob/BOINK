library(tidyverse)
library(lme4)
library(rgdal)
library(raster)
select <- dplyr::select


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
      xlim(-8, NA) +
      facet_grid(popnumber ~ site + tissue)
ggsave(filename = paste0("figures/1.png"), width=8, height=8, units="in")

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
ggsave(filename = paste0("figures/2.png"), width=8, height=8, units="in")
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
      mutate(popnumber=factor(popnumber),
             individual=paste(popnumber, individual))

# some exploratory plots
p %>%
      ggplot(aes(longitude, latitude, color=cwd)) +
      geom_point() +
      scale_color_gradientn(colours=c("blue", "green", "yellow", "red"))
ggsave(filename = paste0("figures/3.png"), width=8, height=8, units="in")


p50 %>%
      select(-embolism) %>%
      spread(tissue, p50) %>%
      ggplot(aes(leaf, stem, color=popnumber, 
                 shape=site, linetype=site)) +
      geom_abline(slope=1, intercept=0) +
      geom_line() +
      geom_point() +
      coord_fixed()
ggsave(filename = paste0("figures/4.png"), width=8, height=8, units="in")


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
ggsave(filename = paste0("figures/5.png"), width=8, height=8, units="in")


p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), cwd=mean(cwd)) %>%
      ggplot(aes(site, p50, group=popnumber, color=popnumber)) +
      geom_line() +
      facet_grid(.~tissue)
ggsave(filename = paste0("figures/6.png"), width=8, height=8, units="in")


p50 %>%
      ggplot(aes(tissue, p50, color=site)) +
      geom_point() +
      facet_grid(.~popnumber)
ggsave(filename = paste0("figures/7.png"), width=8, height=8, units="in")


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
ggsave(filename = paste0("figures/8.png"), width=8, height=8, units="in")


p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), 
                cwd=mean(cwd), djf=mean(djf)) %>%
      ggplot(aes(cwd, djf, color=p50)) +
      geom_point(size=5) +
      facet_grid(tissue ~ site) +
      scale_color_gradientn(colours=c("magenta", "gray20", "green"))
ggsave(filename = paste0("figures/9.png"), width=8, height=8, units="in")



p50 %>%
      group_by(popnumber, site, tissue) %>%
      summarize(p50 = mean(p50), 
                latitude=mean(latitude),
                longitude=mean(longitude)) %>%
      ggplot(aes(latitude, longitude, color=p50)) +
      geom_point(size=5) +
      facet_grid(tissue ~ site) +
      scale_color_gradientn(colours=c("magenta", "gray20", "green"))
ggsave(filename = paste0("figures/10.png"), width=8, height=8, units="in")


# we'd expect stronger djf response in stems than leaves
# i.e. interaction between tissue and temperature
p50 %>%
      group_by(popnumber, site, tissue) %>%
      #summarize(p50 = mean(p50), 
      #          jja=mean(jja), 
      #          djf=mean(djf)) %>%
      gather(var, clim, jja, djf) %>%
      ggplot(aes(clim, p50, color=tissue)) +
      geom_point() +
      geom_smooth(method=lm) +
      facet_grid(site ~ var, scales="free_x")

p50 %>% 
      lmer(p50 ~ djf * tissue + site + (1 | popnumber/individual), data=.) %>%
      summary()


p50 %>% 
      lmer(p50 ~ djf + jja + tissue + site + (1 | popnumber/individual), data=.) %>%
      summary()



p50 %>% 
      lm(p50 ~ djf * cwd * tissue, data=.) %>%
      summary()



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
      summarize_at(vars(p50:aet), max) %>%
      lm(p50 ~ jja + djf + tissue + site, data=.) %>%
      summary()



### quadratic term?

pca <- prcomp(select(p50, djf, jja) %>% scale())
md <- data.frame(clim=pca$x[,1],
                 p50=p50$p50) %>%
      mutate(clim2=clim^2)
m <- lm(p50 ~ clim + clim2, data=md)
md$pred <- predict(m, md)

ggplot(md, aes(clim, pred)) + geom_line()


###################################

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

ggplot(diffs %>% filter(cwd>0), aes(cwd, p50)) + 
      geom_point() + geom_smooth()
ggsave(filename = paste0("figures/10.png"), width=8, height=8, units="in")





#####################

d <- readOGR("E:/edges/range-edges/data/species_ranges/little_trees/Quercus douglasii",
             "querdoug")

climate <- stack("E:/edges/range-edges/data/species_ranges/climate/climate_data.tif") %>%
      crop(d)
names(climate) <- c("AET", "CWD", "PPT", "JJA", "DJF")
climate$JJA <- climate$JJA/10
climate$DJF <- climate$DJF/10

clim <- mask(climate, d)

climd <- clim %>%
      rasterToPoints() %>%
      as.data.frame()


climd %>%
      select(x, y, CWD, JJA, DJF) %>%
      gather(var, value, CWD:DJF) %>%
      group_by(var) %>%
      mutate(value=scales::rescale(value)) %>%
      ggplot() +
      geom_raster(aes(x, y, fill=value)) +
      geom_point(data=p, aes(longitude, latitude), size=3) +
      facet_wrap(~var) +
      theme_void() +
      scale_fill_gradientn(colours=c("darkgreen", "yellow", "darkred")) +
      coord_fixed(ratio=1.2) +
      theme(legend.position="top")
ggsave(filename = paste0("figures/13.png"), width=8, height=8, units="in")





# smoothing neighborhood

library(dismo)
radii <- c(2, 10, 20, 50, 100, 250)
circ <- radii %>%
      lapply(function(x){
            clm <- clim
            if(x==2) clm <- climate
            circles(select(p, longitude, latitude), 
                    d=x*1000, lonlat=T, dissolve=F) %>%
                  as("SpatialPolygons") %>%
                  raster::extract(clm, .) %>%
                  sapply(function(y) apply(y, 2, mean, na.rm=T)) %>%
                  t() %>%
                  as.data.frame() %>%
                  mutate(popnumber=p$popnumber) %>%
                  gather(var, value, -popnumber) %>%
                  mutate(radius=x)
            }) %>%
      do.call("rbind", .) %>%
      left_join(p)

circ %>%
      select(popnumber, var, value, radius) %>%
      distinct() %>%
      spread(var, value) %>%
      select(popnumber, radius, CWD, JJA, DJF) %>%
      ecoclim::pairsData(xy_vars=c("CWD", "JJA", "DJF"),
                         z_vars=c("popnumber", "radius"),
                         mirror=T) %>%
      mutate(x_var = factor(x_var, levels=c("CWD", "JJA", "DJF")),
             y_var = factor(y_var, levels=c("CWD", "JJA", "DJF"))) %>%
      ggplot(aes(x_value, y_value, color=radius, group=radius)) +
      geom_point() +
      geom_smooth(method=lm, se=F) +
      facet_grid(y_var ~ x_var, scales="free") +
      scale_color_gradientn(colours=c("black", "blue", "red", "goldenrod1"),
                            trans="log10")
ggsave(filename = paste0("figures/14.png"), width=8, height=8, units="in")


p50 <- i %>%
      group_by(individual, popnumber,
               site, tissue) %>%
      filter(embolism <= 50) %>%
      arrange(desc(embolism)) %>%
      slice(1) %>%
      ungroup() %>%
      rename(p50 = lwp_MPa) %>%
      left_join(circ) %>%
      mutate(popnumber=factor(popnumber),
             individual=paste(popnumber, individual))

p50 %>%
      group_by(popnumber, site, tissue, var, radius) %>%
      summarize(p50_se = sd(p50)/sqrt(n()),
                p50 = mean(p50), 
                value=mean(value)) %>%
      ggplot(aes(value, p50, ymin=p50-p50_se, ymax=p50+p50_se,
                 color=tissue, linetype=site, shape=site)) +
      geom_point() +
      geom_smooth(method=lm, se=F) +
      facet_grid(radius~var, scales="free") +
      labs(x="climate")
ggsave(filename = paste0("figures/15.png"), width=8, height=8, units="in")



# super basic LM with JJA and DJF, for stems in the wild
p50 %>%
      filter(radius==20, tissue=="stem", site=="Wild") %>%
      spread(var, value) %>%
      lm(p50 ~ JJA + DJF, data=.) %>%
      summary()

# super basic LM with JJA and DJF, across entire dataset
p50 %>%
      filter(radius==20) %>%
      spread(var, value) %>%
      lm(p50 ~ JJA + DJF, data=.) %>%
      summary()

# LMM with JJA and DJF
p50 %>%
      filter(radius==20) %>%
      spread(var, value) %>%
      lmer(p50 ~ DJF + JJA + site + tissue + 
                 (1|popnumber/individual), data=.) %>%
      summary()

# testing for interaction between climate and tissue
p50 %>%
      filter(radius==20) %>%
      spread(var, value) %>%
      lmer(p50 ~ tissue * DJF + tissue * JJA + site + 
                 (1|popnumber/individual), data=.) %>%
      summary()

md <- p50 %>%
      filter(radius==20) %>%
      spread(var, value) %>%
      select(popnumber, individual, site, tissue, p50, JJA, DJF)
m <- md %>% lm(p50 ~ tissue * DJF * JJA, data=.)

pd <- expand.grid(DJF=0:5, JJA=28:35, tissue=c("leaf", "stem")) %>%
      mutate(pred = predict(m, .))

ggplot() +
      geom_contour(data=pd, 
                   aes(DJF, JJA, z=pred, linetype=tissue, color=..level..),
                   size=1) +
      geom_point(data=md, aes(DJF, JJA), color="black", size=3) +
      scale_color_gradientn(colors=c("yellow", "orange", "red", "purple", "blue", "darkblue"))

ggplot(pd, aes(DJF, pred, group=paste(JJA, tissue), color=tissue)) +
      geom_line()
ggplot(pd, aes(JJA, pred, group=paste(DJF, tissue), color=tissue)) +
      geom_line()






### quadratic term?

x <- p50 %>% filter(var %in% c("JJA", "DJF")) %>%
      select(popnumber, individual, site, tissue, p50, var, value, radius) %>%
      spread(var, value) %>%
      split(.$radius)
x <- lapply(x, function(y){
      
      pca <- prcomp(select(y, JJA, DJF) %>% scale())
      md <- y %>% mutate(clim=pca$x[,1]) %>%
            mutate(clim2=clim^2)
      m <- lm(p50 ~ clim + clim2 + site + tissue, data=md)
      md$pred <- predict(m, md)
      md$radius <- y$radius
      return(md)
})
x <- do.call("rbind", x)


ggplot(x, aes(clim, pred, color=factor(radius))) + 
      geom_line() +
      facet_grid(tissue ~ site, scales="free")



# single continuum?

x <- p50 %>% filter(var %in% c("JJA", "DJF")) %>%
      select(popnumber, individual, site, tissue, p50, var, value, radius) %>%
      #spread(var, value) %>%
      split(.$radius)
x <- lapply(x, function(y){
      
      md <- y %>% mutate(clim=scale(value)) %>%
            mutate(clim2=clim^2) %>%
            mutate(clim=as.numeric(clim), clim2=as.numeric(clim2))
      m <- lm(p50 ~ clim + clim2 + site + tissue, data=md)
      md2 <- expand.grid(site=unique(md$site), 
                         tissue=unique(md$tissue),
                        clim=seq(min(md$clim), max(md$clim), length.out=100)) %>%
            mutate(clim2=clim^2)
      md2$pred <- predict(m, md2)
      md$pred <- predict(m, md)
      md2$radius <- y$radius[1]
      md$radius <- y$radius[1]
      return(list(md=md, md2=md2))
})
x2 <- map(x, "md2") %>% do.call("rbind", .)
x <- map(x, "md") %>% do.call("rbind", .)


ggplot() + 
      #geom_line(data=x2, aes(clim, pred, color=factor(radius))) +
      geom_line(data=x, aes(clim, pred, color=factor(radius))) +
      geom_point(data=x, aes(clim, pred, color=factor(radius))) +
      #geom_point(data=x, aes(clim, p50, color=factor(radius))) +
      facet_grid(tissue ~ site, scales="free") +
      labs(x="continuum from winter min to summer max\n(scaled for each radius distance)",
           y="predicted p50",
           color="radius (km)",
           title="quadratic regression, with jja and djf on a single continuum. probably terrible method.\n(weird but maybe cool? suggests both ends of the spectrum may drive p50)")
ggsave(filename = paste0("figures/16.png"), width=8, height=8, units="in")

