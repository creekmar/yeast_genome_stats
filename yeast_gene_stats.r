# Ming Creekmore
# Dr. Babbit
# Bioinformatics Languages
# R Assignment
# Takes in data from a gcContent.txt file and a log2Transcription.txt file
# Creates a bar chart of the mean and SE of this data for the 16 chromosomes
# Runs an ANOVA for both gcContent and log2Transcription
# Creates a linear regression between the mean of gcContent vs. mean of log2Transcription
########################################################################

library(ggplot2)

# Extracting gcContent info from file
gcData <- read.table(file = "gcContent.txt", sep="\t", fill=TRUE)

# Making lists of calculated data for gcContent
# gcDataFrame is the raw Data
means = list()
s_devs = list()
s_errs = list()
rownames <- gcData[,1]
rowID <- seq(1, ncol(gcData)-1)
gcDataFrame = data.frame(rowID="",gc=0)

# populating lists and dataframe
for(x in 1:16) {
  myline <- as.numeric(gcData[x, -1])
  gcDataFrame <- rbind(gcDataFrame, data.frame(rowID=as.character(rownames[x]), gc=na.omit(myline)))
  my_mean <- mean(myline, na.rm = TRUE)
  means <- append(means, my_mean)
  my_sd <- sd(myline, na.rm = TRUE)
  s_devs <- append(s_devs, my_sd)
  my_se <- my_sd / sqrt(length(na.omit(myline)))
  s_errs <- append(s_errs, my_se)
}
# gcstats is the calculated data
gcstats <- data.frame(unlist(means), unlist(s_devs), unlist(s_errs), unlist(rownames))
names(gcstats) <- c("Mean", "SD", "SE", "Yeast_Type")
gcDataFrame = gcDataFrame[-1,]


# Bar chart of GcContent data
ggplot(gcstats, aes(x=Yeast_Type, y=Mean)) +
  ggtitle("GC Content") +
  geom_bar(position=position_dodge(), stat="identity",
          colour='blue', fill="yellow") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2)




# Extracting log2Transcription info from file
logData <- read.table(file = "log2Transcription.txt", sep="\t", fill=TRUE)

# Making lists of calculated data for log2Transcription
# logDataFrame is the raw data for log2 values
means2 = list()
s_devs2 = list()
s_errs2 = list()
rownames <- logData[,1]
rowID <- seq(1, ncol(logData)-1)
logDataFrame = data.frame(rowID="",log=0)

for(x in 1:16) {
  myline <- as.numeric(logData[x, -1])
  logDataFrame <- rbind(logDataFrame, data.frame(rowID=as.character(rownames[x]), log=na.omit(myline)))
  my_mean <- mean(myline, na.rm = TRUE)
  means2 <- append(means2, my_mean)
  my_sd <- sd(myline, na.rm = TRUE)
  s_devs2 <- append(s_devs2, my_sd)
  my_se <- my_sd / sqrt(length(na.omit(myline)))
  s_errs2 <- append(s_errs2, my_se)
}

# logstats is the dataframe of log2 statistics
logstats <- data.frame(unlist(means2), unlist(s_devs2), unlist(s_errs2), unlist(rownames))
names(logstats) <- c("Mean", "SD", "SE", "Yeast_Type")
logDataFrame = logDataFrame[-1,]

# bar chart for log2 values 
ggplot(logstats, aes(x=Yeast_Type, y=Mean)) +
  ggtitle("Log2 Transcription") +
  geom_bar(position=position_dodge(), stat="identity",
          colour='blue', fill="yellow") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2)




# # Perform ANOVA test
print("ANOVA tests for gc Content")
gc_anova <- oneway.test(gc~rowID, gcDataFrame)
print(gc_anova)
print("ANOVA tests for log2 Transcription Content")
log_anova <- oneway.test(log~rowID, logDataFrame)
print(log_anova)




# dataframe of gc means and log2 means
df2 <- data.frame(unlist(means), unlist(means2))
names(df2) <- c("gc", "log")

# # Create a scatterplot of log2 using ggplot2 vs GC
scatterplot <- ggplot(df2, aes(x = gc, y = log)) +
  geom_point() +
  labs(x = "Average GC Content", y = "Average Log2 Transcription Values") +
  ggtitle("Scatterplot of Log2 Transcription vs. GC Content")

# # Add a linear regression line to the scatterplot
scatterplot_with_regression <- scatterplot +
  geom_smooth(method = "lm", se = FALSE)

scatterplot_with_regression
