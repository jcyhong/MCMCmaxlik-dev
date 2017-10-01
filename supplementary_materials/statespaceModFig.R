require(ggplot2)
require(reshape)
iteratesFixed <- data.frame(
  iter=1:nrow(resultsSSModFixed$param),
  resultsSSModFixed$param,
  method=rep("Fixed (0.05)", nrow(resultsSSModFixed$param))
)
iteratesSmallFixed <- data.frame(
  iter=1:nrow(resultsSSModSmallFixed$param),
  resultsSSModSmallFixed$param,
  method=rep("Fixed (0.005)", nrow(resultsSSModSmallFixed$param))
)
# iteratesLogregNR <- data.frame(
#   iter=1:nrow(resultsSSModNR$param),
#   resultsSSModNR$param,
#   method=rep("Newton-Raphson", nrow(resultsSSModNR$param))
# )
iterates1D <- data.frame(
  iter=1:nrow(resultsSSMod1D$param),
  resultsSSMod1D$param,
  method=rep("1D sampling", nrow(resultsSSMod1D$param))
)
iteratesAdadelta <- data.frame(
  iter=1:nrow(resultsSSModAdadelta$param),
  resultsSSModAdadelta$param,
  method=rep("Adadelta", nrow(resultsSSModAdadelta$param))
)
iteratesAdam <- data.frame(
  iter=1:nrow(resultsSSModAdam$param),
  resultsSSModAdam$param, 
  method=rep("Adam", nrow(resultsSSModAdam$param))
)

iterates <- rbind(#iteratesFixed,
  iteratesSmallFixed,
                 # iteratesLogregNR,
                  iteratesAdadelta, iteratesAdam,
                  iterates1D) 
names(iterates) <- c("iter", paramNodesSS, "method")

write.csv(iterates, "ssModResultsAll.csv", row.names=F)

iteratesMelted <- melt(iterates, id.vars=c("iter", "method"))
levels(iteratesMelted$variable) <- c(expression(rho),
                                     expression(mu),
                                     expression(sigma[PN]),
                                     expression(sigma[OE]))
library(RColorBrewer)
# cbPalette <- brewer.pal(8, "Dark2")
# 
# 
# trajectoryPlot <- ggplot(iteratesMelted, aes(x=iter, y=value, colour=method)) +xlim(0,500)
# trajectoryPlot <- trajectoryPlot + 
#   geom_line(alpha=0.8) + 
#   geom_point(data=iteratesMelted, mapping=aes(x=iter, y=value, shape=method))
# trajectoryPlot<- trajectoryPlot + 
#   facet_grid(variable ~ ., labeller=label_parsed, scales="free_y")
# trajectoryPlot <- trajectoryPlot + theme(legend.position="none")
# trajectoryPlot <- trajectoryPlot + 
#   scale_color_manual(values=cbPalette)
# trajectoryPlot <- trajectoryPlot + xlab("number of iterations")
# trajectoryPlot<- trajectoryPlot + 
#   theme(axis.title=element_text(face="bold", size=20),
#         axis.text=element_text(size=15),
#         strip.text.y=element_text(size=15))



cbPalette <- brewer.pal(8, "Dark2")
trajectoryPlot <- ggplot(iteratesMelted,
                         aes(x=iter, y=value, colour=method)) +xlim(0,300)
trajectoryPlot <- trajectoryPlot + geom_line(alpha=0.8) + 
  geom_point(data=iteratesMelted, mapping=aes(x=iter, y=value, shape=method))
trajectoryPlot <- trajectoryPlot + 
  facet_grid(variable ~ .,
             labeller=label_parsed,
             scales="free_y")
trajectoryPlot<- trajectoryPlot
trajectoryPlot <- trajectoryPlot + 
  scale_color_manual(values=cbPalette,labels=unique(iterates$method)) + 
  scale_shape_manual(values=c(0, 5, 6, 15,8),labels=unique(iterates$method))
trajectoryPlot<- trajectoryPlot + xlab("number of iterations")
trajectoryPlot <- trajectoryPlot + 
  theme(axis.title=element_text(face="bold", size=20),
        axis.text=element_text(size=15),
        strip.text.y=element_text(size=15)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank(), legend.text=element_text(size=15))

pdf(paste0("statespace_", 300, ".pdf"))
trajectoryPlot
dev.off()

