# ======================================= FUNCTIONS FOR SEAHORSE PIPELINE
# Author: Matej Stevuliak

dm <- dm$rates
sam <- dm %>%
  group_by(sample_id) %>%
  summarise(N = n())


samf <- d_OCR %>%
  group_by(sample_id) %>%
  summarise

write_csv(sam, "final-sample_ids.csv")

