library(glmm)
data(salamander)


#### organize original data ####

orgF <- cbind.data.frame(Findex = 1:length(unique(salamander$Female)), idF = unique(salamander$Female))

orgM <- cbind.data.frame(Mindex = 1:length(unique(salamander$Male)), idM = unique(salamander$Male))

salamander2 <- merge(salamander, orgF, by.x = "Female", by.y = "idF")
salamander3 <- merge(salamander2, orgM, by.x = "Male", by.y = "idM")


#### n = 180 ####

filterF <- salamander3 %>% filter(Findex %in% 1:30)

filterM <- salamander3 %>% filter(Mindex %in% 1:30) # %>% group_by(Mindex) %>% summarise(count = n()) %>% as.data.frame()

setwd("~/Desktop/MCMCmaxlik-dev/examples/salamander_new")
write.csv(filterF, "first30F_n180.csv", row.names = F)
write.csv(filterM, "first30M_n180.csv", row.names = F)

#### n = 720 ####

newF <- cbind.data.frame(Findex_og = unique(salamander3$Findex), Findex_new = 61:120)
newM <- cbind.data.frame(Mindex_og = unique(salamander3$Mindex), Mindex_new = 61:120)

salamander4 <- merge(salamander3, newF, by.x = "Findex", by.y = "Findex_og")

salamander5 <- merge(salamander4, newM, by.x = "Mindex", by.y = "Mindex_og")

top <- salamander5[, c("Mate", "Cross", "Mindex", "Findex")]
bottom <- salamander5[, c("Mate", "Cross", "Mindex_new", "Findex_new")]

names(bottom) <- names(top)

big_data <- rbind.data.frame(top, bottom)

write.csv(big_data, "n720.csv", row.names = F)

#### n = 1440 ####

newF <- cbind.data.frame(Findex_og = unique(big_data$Findex), Findex_new = 121:240)
newM <- cbind.data.frame(Mindex_og = unique(big_data$Mindex), Mindex_new = 121:240)

salamander4b <- merge(big_data, newF, by.x = "Findex", by.y = "Findex_og")

salamander5b <- merge(salamander4b, newM, by.x = "Mindex", by.y = "Mindex_og")

top <- salamander5b[, c("Mate", "Cross", "Mindex", "Findex")]
bottom <- salamander5b[, c("Mate", "Cross", "Mindex_new", "Findex_new")]

names(bottom) <- names(top)

big_data2 <- rbind.data.frame(top, bottom)

write.csv(big_data2, "n1440.csv", row.names = F)

