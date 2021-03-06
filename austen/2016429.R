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
library("jsonlite")

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

books <- original_books %>%
  anti_join(stop_words)

bing <- sentiments %>%
  filter(lexicon == "bing", word != "miss") %>%
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
  right_join(original_books) %>%
  group_by(book, index, linenumber) %>%
  mutate(word_ix = seq_len(n()) - 1, nchar = nchar(word) + 1, nchar_sum = cumsum(nchar),
         nchar_offset = c(0, nchar_sum[-length(nchar_sum)])) %>%
  group_by(book, index, sentiment) %>%
  mutate(sentiment_ix = seq_len(n()))

glimpse(book_index, width = 500)

## ---- write-data ----
p <- ggplot(book_index %>% filter(book == "Mansfield Park", index == 179) %>% replace_na(list("sentiment" = "unknown"))) +
  geom_text(aes(y = -linenumber, x = nchar_offset, label = word, col = sentiment), size = 3.8,  hjust = 0, family = "Courier")

p <- ggplot(book_index %>% filter(book == "Pride & Prejudice", linenumber < 30) %>% replace_na(list("sentiment" = "unknown"))) +
  geom_text(aes(y = -linenumber, x = nchar_offset, label = word, col = sentiment), size = 3.8,  hjust = 0, family = "Courier")

book_small <- rbind(
  book_index %>%
  filter(book == "Mansfield Park", index == 179),
  book_index %>%
  filter(book == "Pride & Prejudice", index < 3)
) %>%
  replace_na(list("sentiment" = "none"))

cat(sprintf("var sentiment = %s;", toJSON(sentiment_count, auto_unbox = FALSE)), file = "~/Desktop/lab_meetings/20170519/sentiment.js")
cat(sprintf("var book = %s;", toJSON(book_small, auto_unbox = FALSE)), file = "~/Desktop/lab_meetings/20170519/book.js")

book_index %>%
  filter(linenumber == 10, index == 0) %>%
  glimpse()

book_index %>%
  filter(index == 1, book == "Pride & Prejudice", sentiment == "positive") %>%
  glimpse()

## Write top of document term matrix
passage_term <- book_index %>%
  filter(!(word %in% stop_words$word)) %>%
  ungroup() %>%
  select(-sentiment) %>%
  count(book, index, word) %>%
  arrange(desc(n)) %>%
  spread(word, n, fill = 0)

passage_term <- passage_term[passage_term$book == "Pride & Prejudice", ]
col_sums <- apply(passage_term[, -c(1, 2)], 2, sum)
passage_sub <- passage_term[1:20, c(2, 1, 2 + order(col_sums, decreasing = TRUE)[1:8])]
cat(sprintf("var document_term = %s", toJSON(passage_sub, auto_unbox = FALSE)), file = "~/Desktop/lab_meetings/20170519/document_term.js")

rownames(passage_sub) <- NULL
passage_sub <- passage_sub %>%
  mutate_at(.vars = vars(-book), as.integer)
cat(print(xtable(passage_sub), include.rownames = FALSE), file = "~/Desktop/lab_meetings/20170911/dtm.tex")

## Write top of microbiome counts matrix
library("treelapse")
library("phyloseq")
data(abt)
sample_data(abt)
taxa_order <- order(taxa_sums(abt), decreasing = TRUE)
taxa_sub <- cbind(
  sample_data(abt)[1:20, c("ind", "time")],
  t(get_taxa(abt))[1:20, taxa_order[1:7]]
) %>%
  rename(subject = ind)

taxa_sub <- taxa_sub[, c(2, 1, 3:9)]
rownames(taxa_sub) <- NULL
taxa_sub$time <- as.integer(taxa_sub$time)
cat(print(xtable(taxa_sub[, 1:7]), include.rownames = FALSE), file = "~/Desktop/lab_meetings/20170911/taxa.tex")
cat(sprintf("var sample_rsv = %s", toJSON(taxa_sub, auto_unbox = FALSE)), file = "~/Desktop/lab_meetings/20170519/sample_rsv.js")

## topic evolution data
n <- 250
time <- seq(0, 1, length.out = n)
y <- sin(time * 2 * pi * 2.2) #+ runif(100, -0.2, 0.2)
x <- cos(time * 2 * pi * 3) #+ runif(100, -0.2, 0.2)
z <- (sin(time * 2 * pi * 3)) ^ 2 #+ runif(100, -0.2, 0.2)

p_data <- apply(
  cbind(x, y, z), 1,
  function(a) exp(a) / sum(exp(a))
) %>% t()

A <- matrix(c(0, 0.5, 1, 0, 1, 0, 0, 0, 0), byrow = F, ncol = 3)
p_data <- p_data %*% A
plot(A[, 1:2], asp = 1)
points(p_data[, 1:2], asp = 1)

cat(sprintf("var sketch_data = %s", toJSON(data.frame(i = 0:(n - 1), x = p_data[, 1], y = p_data[, 2]))), file = "~/Desktop/lab_meetings/20170519/sketch_data.js")
