# ======================================
# Load Package
# ======================================

library(mtf)
library(tidyverse)



# # Test midge pulse:
# plot(test_midges(1000, a=1e9, b=10, r=400, w=40, d=1), type = 'l', ylab = "iM")

midge_pulse <- function(t, b, s, w) {
    f = ifelse(t > s & t <= s + w, b, 0)
    return(f)
}

curve(midge_pulse(x, b = 10, s = 10, w = 15), 0, 100)

web_output <- food_web(tmax = 100, b = 20, s = 10, w = 15)


# Plot absolute N
web_output %>%
    group_by(pool) %>%
    # Define 'minb' to set the minimum value for the y-axis.
    # This allows different y-scales for different facets, with the ymin set to 0
    mutate(minb = 0) %>%
    ggplot(aes(time, N)) +
    facet_wrap(~pool, scales="free_y") +
    # The horizontal lines show the initial states
    geom_hline(data = web_output %>% filter(time == min(time)),
               aes(yintercept=N), color="firebrick") +
    geom_line(size = 1) +
    geom_point(aes(time, minb), shape="") +
    theme_classic()



# ======================================
# Load Package
# ======================================


d = crossing(area = c(10, 100, 1000),
           w = c(10, 15, 20, 25, 30)) %>%
    mutate(b = area/w)


yy = lapply(1:nrow(d), function(x){
    xx = d[x,]
    y = food_web(tmax = 100,  s = 10, b = xx$b, w = xx$w, .lM = 0.2)
    y = y %>%
        mutate(area = xx$area, b = xx$b, w = xx$w, N = round(N, 10))
    return(y)
}) %>%
    bind_rows()

# fig 1 (through time, area = 10)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 10) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+
    geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()

# fig 2 (through time, area = 100)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 100) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()

# fig 3 (through time, area = 1000)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 1000) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_line(size = 0.8)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()




# fig 2 (through time, area = 100)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    group_by(pool) %>%
    mutate(N = (N - mean(N))/sd(N)) %>%
    ggplot(aes(time, N, color = pool, linetype = pool))+
    facet_grid(area~w, scales= "free_y")+
    geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values = c("gray", "burlywood4", "forestgreen", "burlywood4", "forestgreen", "black", "dodgerblue"))+
    scale_linetype_manual(values = c(1, 2, 2, 1,1, 1, 2))+
    theme_classic()


yy2 = lapply(c(1e-6, 1e-3, 1, 100), function(x){
    y = food_web(tmax = 100,  s = 10, b = 5, w = 20, .lM = 0.2, .aM = x)
    y = y %>%
        mutate(aM = x, N = round(N, 10))
    return(y)
}) %>%
    bind_rows()

yy2 %>%
    group_by(pool) %>%
    mutate(N = (N - mean(N))/sd(N)) %>%
    ggplot(aes(time, N, color = pool, linetype = pool))+
    facet_wrap(~aM, scales= "free_y", nrow=3)+
    geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values = c("gray", "burlywood4", "forestgreen", "burlywood4", "forestgreen", "black", "dodgerblue"))+
    scale_linetype_manual(values = c(1, 2, 2, 1,1, 1, 2))+
    theme_classic()






yy %>%
    filter(time==75, pool != "midge") %>%
    ggplot(aes(w, N, color = factor(area)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_line(size = 0.8)+
    scale_color_manual(values=c("gray50","dodgerblue","firebrick"))+
    scale_x_continuous("Duration")+
    theme_classic()

yy %>%
    filter(pool != "midge", area == 10) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    theme_classic()




# ======================================
# Load Package
# ======================================


d = crossing(area = c(10, 50, 250),
             w = c(10, 15, 20, 25, 30)) %>%
    mutate(b = area/w)


yy = lapply(1:nrow(d), function(x){
    xx = d[x,]
    y = food_web(tmax = 100,  s = 10, b = xx$b, w = xx$w, .lM = 0.2)
    y = y %>%
        mutate(area = xx$area, b = xx$b, w = xx$w, N = round(N, 10))
    return(y)
}) %>%
    bind_rows()

# fig 1 (through time, area = 10)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 10) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+
    geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()

# fig 2 (through time, area = 50)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 50) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+geom_line(size = 0.8)+
    scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()

# fig 3 (through time, area = 250)
yy %>%
    filter(w %in% c(10,20,30)) %>%
    filter(pool != "midge", area == 250) %>%
    ggplot(aes(time, N, color = factor(w)))+
    facet_wrap(~pool, scales= "free_y")+
    geom_hline(data = yy %>% filter(time==0,pool != "midge"),
               aes(yintercept = N), size = 0.3)+
    geom_line(size = 0.8)+
    geom_hline(data = yy %>%
                   group_by(pool) %>%
                   filter(N == max(N),pool != "midge"),
               aes(yintercept = N), size = 0)+scale_x_continuous("Time")+
    scale_color_manual(values=c("black","burlywood4","forestgreen"))+
    theme_classic()





