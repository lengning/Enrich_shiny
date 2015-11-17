library(shiny)
library(gdata)
library(shinyFiles)

#####################
# allez function
#####################


enrich_allez <- function(File, FileSep, Local, LocalSep, outprefix, lib.v, side,  Lowersetsize, Uppersetsize, pcut ){

# csv or txt
FileType=FileSep

if(FileType=="csv"){
	cat("\n Read in csv file \n")
	prefix=strsplit(File,split="\\.csv")[[1]][1]
	In=read.csv(File,stringsAsFactors=F,row.names=1)
}
if(FileType!="csv"){
	cat("\n Read in tab delimited file \n")
	prefix=strsplit(File,split=paste0("\\.",FileType))[[1]][1]
	In=read.table(File,stringsAsFactors=F,row.names=1)
	if(!is.numeric(In[[1]]))In=read.table(File,stringsAsFactors=F,row.names=1,header=T)
}

if(!is.null(Local)){
	FileType2=LocalSep

	if(FileType2=="csv"){
		cat("\n Read in csv file (gene list)\n")
		ListIn=read.csv(Local,stringsAsFactors=F,row.names=1)
	}
	if(FileType2!="csv"){
		cat("\n Read in tab delimited file (gene list)\n")
		ListIn=read.table(Local,stringsAsFactors=F,row.names=1)
	}
	if(nrow(ListIn)>1){
	List=sapply(1:nrow(ListIn),function(i)setdiff(as.vector(ListIn[i,]),c(""," ")))
	}
	if(nrow(ListIn)==1) {
		List=vector("list",1)
		List[[1]]=setdiff(unlist(ListIn),c(""," "))}

	names(List)=rownames(ListIn)

} else List=NULL

library(allez)
Score=In[[1]]
names(Score)=rownames(In)

#browser()
Out=allez(score=Score,lib=lib.v,idtype="SYMBOL",locallist=List)

Mat=Out$setscores[,c("Term","set.mean","set.sd","set.size","z.score")]
message(c("one tailed p value?", side))
if(side=="F")Mat$p.value <- pnorm(-abs(Mat$z.score))# two tailed
if(side=="T"){
	prb <- pnorm(Mat$z.score)# one tailed
	Mat$p.value <- ifelse(1-prb>prb, prb, 1-prb)*2
}
Mat$p.adj <- p.adjust(Mat$p.value, method="BH")
Mat <- Mat[which(Mat$set.size>Lowersetsize),]
Mat <- Mat[which(Mat$set.size<Uppersetsize),]
MatOut <- Mat[order(Mat$p.value),c("Term","p.value","p.adj","z.score","set.size","set.mean","set.sd")]


message("sets with size < ",Lowersetsize, " or > ", Uppersetsize, " are not considered" )
LocalOut <- MatOut[which(is.na(MatOut[,"Term"])),]
MatOut2 <-  cbind(rownames(MatOut), MatOut)
LocalOut2 <- cbind(rownames(LocalOut), LocalOut)
colnames(MatOut2)[1] = colnames(LocalOut2)[1] = "GO_ID"

write.table(MatOut2,file=paste0(outprefix,"_enrichment_allsets.txt"),sep="\t", row.names=F)
write.table(LocalOut2,file=paste0(outprefix,"_enrichment_localsets.txt"), sep="\t", row.names=F)

Mat.p <- MatOut2[which(MatOut2$p.adj<=pcut),]
Local.p <- LocalOut2[which(LocalOut2$p.adj<=pcut),]
write.table(Mat.p,file=paste0(outprefix,"_enrichment_allsets_sig.txt"),sep="\t", row.names=F)
write.table(Local.p,file=paste0(outprefix,"_enrichment_localsets_sig.txt"), sep="\t", row.names=F)

Out=list(Allres=MatOut2, Localres=LocalOut2, SigAllres=Mat.p, SigLocalres=Local.p)
}







# Define server logic for slider examples
shinyServer(function(input, output, session) {
		volumes <- c('home'="~")
		shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
						
		In <- reactive({
	
		print(input$Outdir)
		outdir <- paste0("~/",input$Outdir[[1]][[2]],"/")
		print(outdir)

		the.file <- input$filename$name #input
		file.toread <- input$filename$datapath	
		file.sep0 <- strsplit(the.file,split="\\.")[[1]]
		file.sep <- file.sep0[length(file.sep0)]

		Group.file <- input$markerlist$name # markerlist
  	GroupB <- ifelse(is.null(Group.file),FALSE,TRUE)
		mklist.toread <- NULL
		mklist.sep <- NULL
		if(GroupB==TRUE){
			mklist.toread <- input$markerlist$datapath
			mklist.sep0 <- strsplit(Group.file,split="\\.")[[1]]
			mklist.sep <- mklist.sep0[length(mklist.sep0)]
		}


		method.v <- c("allez","EACI","EASE")
		List <- list(
		Input = file.toread, Inputsep=file.sep,
		mkInput = mklist.toread, mkSep=mklist.sep,
		Dir=outdir, out_pre=input$exFileName,
		lower=as.numeric(input$llsize),
		upper=as.numeric(input$ulsize),
		species=ifelse(input$species_buttons=="1","human","mouse"),
		method_use = method.v[as.numeric(input$method_buttons)],
		one_tail=ifelse(input$tail_buttons=="1","T","F"),
		pval_cutoff=input$pcut
)
		# main function
  	str(List)
		db.v <- List$species
		if(db.v=="human")lib.v <- "org.Hs.eg"
		if(db.v=="mouse")lib.v <- "org.Mm.eg"
		
		if(List$method_use=="allez"){
		Res <- enrich_allez(File=List$Input, FileSep=List$Inputsep,
												outprefix=paste(List$Dir,List$out_pre,sep="/"),
							Local=List$mkInput, LocalSep=List$mkSep,
							lib.v=lib.v, side=List$one_tail,  
							Lowersetsize=List$lower, Uppersetsize=List$upper,
							pcut=List$pval_cutoff)
		
		}

	 	List=c(List, Res)	
}) 

  Act <- eventReactive(input$Submit,{
		      In()})
	# Show the values using an HTML table
  output$print0 <- renderText({
				tmp <- Act()
				paste("# significant sets:",
							nrow(tmp$SigAllres),
							"; output directory:", tmp$Dir)
  })

	output$tab <- renderDataTable({
		tmp <- Act()$SigAllres
		t1 <- tmp
		print("done")
		t1
		},options = list(lengthManu = c(4,4), pageLength = 20))

#	output$done <- renderText({"Done"})
})
