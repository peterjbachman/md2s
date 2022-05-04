## ---------------------------------------------------------------------------------------------------
##
## 		Replication Files: Plot Figure 3 (Top Panel)
##		"Scaling Data from Multiple Sources"
##		Authors: Ted Enamorado, Gabriel Lopez-Moctezuma, Marc Ratkovic
## 
## ---------------------------------------------------------------------------------------------------
source('./00_MD2S.R')

## --------------------------------------------------------------
## Panel A
## --------------------------------------------------------------

color.plot<-TRUE

load('./results_Panel_A/outall_Panel_A.RData')
pdf('Figure3_Panel_A.pdf',h=6,w=18)

x.plot<-c(1:5, (7:11), 13:17, 19:21)
par(mar=c(3.5,3.5,2,3.5))
plot(c(0,0),ylim=c(0,1),xlim=c(1,21),xlab="",ylab="",xaxt="n",yaxt="n",type="n")
segments(x0=c(6,12,18),x1=c(6,12,18),y0=0, y1=1.1,col=gray(.8),lwd=6)
segments(x0=-3,x1=18,y0=.1,y1=.1,col=gray(.25),lty=2,lwd=2)

if(color.plot) {
  col1<-rgb(1,0,0); col2<-rgb(0,0,1)
} else{
  col1<-gray(.25); col2<-gray(.75)
}
col.isreal<-c(rep(col1,2),rep(col2,3), rep(col1,2),rep(col2,3), rep(col1,3),rep(col2,2))

for(i in 1:15){
  boxplot(out.all[,c(4:8, 14:18, 9:13)][,i], at=x.plot[i], add=TRUE, yaxt="n", col=col.isreal[i],
          pch=19, cex=.7, outcol = col.isreal[i], boxcol = col.isreal[i],
          border = col.isreal[i], outcex = 0.45)
}

axis(side=2,at=c(seq(0,1,by=0.1)),cex.axis=1.5)
axis(side=1,at=x.plot[1:15],label=c(rep(1:5,3)),cex.axis=1.5)

mtext(side=3,at=3,"Shared Subspace",cex=1.5,font=2,line=.1)
mtext(at=9, side=3, expression(bold(paste("Idiosyncratic, ",Y[1]))),cex=1.5,font=2)
mtext(at=15, side=3, expression(bold(paste("Idiosyncratic, ",Y[2]))),cex=1.5,font=2)

mtext(side=2,at=.5,line=2.35, "Estimated p-value",cex=1.5,font=2)

if(!color.plot) mtext(side=1,at=c(1:2,7:8,13:15,19),"--",line=1.65,cex=1.5)

mtext(side=1,at = c(3,9,15),"Dimension",line=2.5,cex=1.5)
legend(x=18.75,y=1.025, c("Systematic","Noise"),col=c(col1,col2),pch=15,border="white",bty="n",cex=1.3,title="Dimension",pt.cex=1.8)

out.all2 <- out.all
out.all2[out.all < 0.10] <- 1
out.all2[out.all >= 0.10] <- 0

colMeans(out.all2)
dev.off()
