# Copyright 2022 Twitter, Inc.
# SPDX-License-Identifier: Apache-2.0

require(gridExtra)
require(dplyr)
require(ggplot2)
require(reshape2)

set.seed(2022)
K.language <- 40
n <- 3000
K.age <- 15

#make two categories, age and language with language having many more values
age.probs <- c(15, 18, 20, 20, 50, 75, 50, 20, 13, 13, 5, 5, 5, 5, 4)
age.probs <- age.probs/sum(age.probs)
age <- sample(1:K.age, n, replace = TRUE, prob = age.probs)

#hist(table(age))
#table(age)

language.probs <- c(rep(20, 20), rep(100, 10), rep(200, 5), rep(400, 5))
language.probs <- language.probs/sum(language.probs)

language <- sample(1:K.language, n, replace = TRUE, prob = language.probs)

#hist(table(language))
#sort(table(language))


dat <- data.frame(Y1 = rbinom(n, 1, .8), Y2 = rbinom(n, 1, .8), Y3 = rbinom(n, 1, .8), age, language)   

age.charts <- dat %>% group_by(age) %>%
  summarize(`Model 1` = mean(Y1), `Model 2` = mean(Y2), `Model 3` = mean(Y3))

age.charts <- melt(age.charts, id.vars = 'age', value.name = "accuracy")
ageplot <- ggplot(age.charts, aes(x=age, y = accuracy)) + geom_bar(stat = 'identity') + 
         facet_grid(vars(variable)) + 
  theme_bw()  +
  theme(text = element_text(size = 20)) 

language.charts <- dat %>% group_by(language) %>%
  summarize(`Model 1` = mean(Y1), `Model 2` = mean(Y2), `Model 3` = mean(Y3))
language.charts <- melt(language.charts, id.vars = 'language', value.name = "accuracy")


languageplot <- ggplot(language.charts, aes(x=language, y = accuracy)) + geom_bar(stat = 'identity') + 
  facet_grid(vars(variable))+
  theme_bw() + 
  theme(text = element_text(size = 20))   



png(file = '~/Desktop/toy-example.png', height = 300, width = 1000)
plt <- grid.arrange(ageplot, languageplot, nrow = 1)
dev.off()

grid.arrange(ageplot, languageplot, nrow = 1)

age.charts %>% group_by(variable) %>% mutate(mean = mean(value)) %>% summarize(var = var(value), mad = mean(abs(value - mean)))
160