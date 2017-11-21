
args<-commandArgs();
if(length(args)!=7) {
   cat("USAGE: R --slave --args raw_data out_fig fig_width < Guidance_Plot.R\n");
   q();
}
raw_data <- args[4];
out_fig <- args[5];
x_width <- as.integer(args[6]);
Program <-args[7];
#Note that the args vector contains everything you wrote in the commandline.  In this case:
#args[1]:  /usr/lib64/R/bin/exec/R
#args[2]:  --slave
#args[3]:  --args
#args[4]:  raw_data
#args[5]:  out_fig
#args[5]:  fig_width

fig_width=15*x_width;
x_thick_num=x_width/10;
print (raw_data);
print (out_fig);
print (fig_width);

guidance_col <- read.csv(file=raw_data,head=TRUE,sep=",");
png(out_fig, width=fig_width,height=500);
plot ( guidance_col[,1], guidance_col[,2],type="b",col = "blue",ylab=paste(Program," Score"),xlab="Column",ylim=c(0,1),xlim=c(0,x_width),xaxp=c(0,x_width,x_thick_num),yaxp=c(0,1,10));
dev.off()
