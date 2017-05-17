#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Following allow from http://juliasilge.com/blog/Life-Changing-Magic/

## ---- libraries ----
library("tidyr")
library("ggplot2")
library("tidytext")
library("janeaustenr")
library("dplyr")
library("stringr")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_rect(size = 0.3, fill = "transparent"),
  panel.grid = element_blank(),
  text = element_text(family = "Arial", color = "#22211d"),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

## ---- start-analysis ----
data("stop_words")
data("sentiments")

original_books <- austen_books() %>%
  group_by(book) %>%
  mutate(
    linenumber = row_number(),
    chapter = cumsum(
      str_detect(text, regex("^chapter [\\divxlc]", ignore_case = TRUE))
    )
  ) %>%
  ungroup()

original_books <- original_books %>%
  mutate(index = linenumber %/% 80) %>%
  unnest_tokens(word, text)

books <- books %>%
  anti_join(stop_words)

bing <- sentiments %>%
  filter(lexicon == "bing") %>%
  select(-score)

sentiment_labels <- books %>%
  inner_join(bing) %>%
  arrange(book, linenumber)

sentiment_count <- sentiment_labels %>%
  count(book, index, sentiment) %>%
  spread(sentiment, n, fill = 0) %>%
  mutate(sentiment = positive - negative)

## order books
book_order <- sentiment_count %>%
  group_by(book) %>%
  summarise(mean = mean(sentiment)) %>%
  arrange(desc(mean)) %>%
  .[["book"]]
sentiment_count$book <- factor(
  sentiment_count$book,
  levels = book_order
)

## plot results
ggplot(sentiment_count) +
  geom_bar(
    aes(x = index, y = sentiment, fill = sentiment),
    stat = "identity"
  ) +
  scale_fill_gradient2(low = "#6f426f", high = "#426f6f") +
  facet_wrap(~book, ncol = 2, scales = "free_x")

## well damn
sentiment_labels %>%
  filter(
    book == "Mansfield Park",
    index == 179
  )

## the most depressing passage
book_index <- sentiment_labels %>%
  select(linenumber, index, word, sentiment) %>%
  right_join(original_books)

book_index %>%
  filter(book == "Mansfield Park", index == 179) %>%
  .[["word"]] %>%
  paste(collapse = " ")
