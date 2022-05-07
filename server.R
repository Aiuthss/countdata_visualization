#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(pathview)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  
  output$textOutput = renderText(input$textInput)
  output$MA_title <- renderText("MAPlot")
  output$vol_title <- renderText("Volcano Plot")
  output$plot_title = renderText(paste(input$group_name1, "vs", input$group_name2))
  
  observeEvent(input$confirm_button, {
    count_table <<- reactive(read.csv(input$fileInput$datapath, header = T, row.names = 1) %>% mutate_all(as.integer))
    output$grouping_input = renderUI({
          output = tagList()
          output[[1]] = textInput(inputId="group_name1", label="input group1 name", value = "Group A")
          output[[2]] = multiInput(inputId="group1", label=NULL, choiceNames=colnames(count_table()), choiceValues=1:(length(count_table())))
          output[[3]] = textInput(inputId="group_name2", label="input group2 name", value = "Group B")
          output[[4]] = multiInput(inputId="group2", label=NULL, choiceNames=colnames(count_table()), choiceValues=1:(length(count_table())))
          output[[5]] = actionButton(inputId = "drawing_button", label="draw")
          output[[6]] <- textInput(inputId="gene", label="select gene", value="gene name")
          output
    })
    
    dataset = switch(isolate(input$strain),
      "MMusculus" = "mmusculus_gene_ensembl",
      "MAuratus" = "mauratus_gene_ensembl",
      "HSapiens" = "hsapiens_gene_ensembl")
    
    hd <<- readRDS(paste0("./biomaRt/", dataset, ".obj"))
    print(hd)
    id_symbol <<- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                       filters = "ensembl_gene_id", values = rownames(count_table()), 
                       mart = hd, useCache = FALSE)
    output$table = renderDataTable(right_join(id_symbol, count_table() %>% mutate(ensembl_gene_id=rownames(count_table())), by="ensembl_gene_id"))
  })
  
  observeEvent(input$drawing_button, {
    
    output$plot_parameter <- renderUI({
      output <- tagList()
      output[[1]] <- textOutput(outputId = "MA_title")
      output[[2]] <- sliderInput(inputId="MA_xmin", label="X-min", min=-20, max= 20, value=-5)
      output[[3]] <- sliderInput(inputId="MA_xmax", label="X-max", min=-20, max= 20, value=15)
      output[[4]] <- sliderInput(inputId="MA_ymin", label="Y-min", min=-15, max= 15, value=-8)
      output[[5]] <- sliderInput(inputId="MA_ymax", label="Y-max", min=-15, max= 15, value=8)
      output[[6]] <- textOutput(outputId = "vol_title")
      output[[7]] <- sliderInput(inputId="vol_xmin", label="X-min", min=-20, max= 20, value=-10)
      output[[8]] <- sliderInput(inputId="vol_xmax", label="X-max", min=-20, max= 20, value=10)
      output[[9]] <- sliderInput(inputId="vol_ymin", label="Y-min", min=0, max= 15, value=0)
      output[[10]] <- sliderInput(inputId="vol_ymax", label="Y-max", min=0, max= 300, value=100)
      output[[11]] <- actionButton(inputId = "parameter_button", label="refresh")
      output
    })
    
    con = factor(c(rep(input$group_name1, length(input$group1)), rep(input$group_name2, length(input$group2))))
    group = data.frame(con = con)
    
    table2 = count_table()[,c(as.integer(input$group1), as.integer(input$group2))]
    avOut <- apply(table2, 1, mean)
    table2 = table2[avOut>0,]
    print(head(table2))
# overexpression analysis by DESeq2    
    x = as.matrix(table2)
    dds = DESeqDataSetFromMatrix(countData = x, colData = group, design = ~con)
    dds = DESeq(dds)
    res_id <<- as.data.frame(results(dds))
    res <<- right_join(id_symbol, res_id %>% mutate(ensembl_gene_id=rownames(res_id)), by= "ensembl_gene_id") %>% drop_na(ensembl_gene_id)
    output$res <- renderDataTable(res)
    
    plot1 <- ggplot(res, aes(x=log2(baseMean), y=log2FoldChange, colour=factor(padj<0.01))) +
      geom_point(size=2, alpha=0.1) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(-5, 15), ylim=c(-8, 8))
    output$plot1 <<- renderPlot(plot1)
    plot2 = ggplot(res, aes(x=log2(baseMean), y=log2FoldChange, colour=factor(external_gene_name==input$gene), alpha=factor(external_gene_name==input$gene))) +
      geom_point(size=2) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      scale_alpha_manual(values=c(0.05, 1)) +
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(-5, 15), ylim=c(-8, 8))
    output$plot2 = renderPlot(plot2)
      
    vol_plot1 <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), colour=factor(padj<0.01))) +
      geom_point(size=2, alpha=0.1) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(-8, 8), ylim=c(0, 30))
    output$vol_plot1 <- renderPlot(vol_plot1)
    vol_plot2 <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj),colour=factor(external_gene_name==input$gene), alpha=factor(external_gene_name==input$gene))) +
      geom_point(size=2) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      scale_alpha_manual(values=c(0.05, 1)) +
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(-8, 8), ylim=c(0, 30))
    output$vol_plot2 <- renderPlot(vol_plot2)
    
    # KEGG, REACTOME
    OrgDb <<- switch(isolate(input$strain),
                     "MMusculus" = org.Mm.eg.db,
                     "MAuratus" = org.Mm.eg.db,
                     "HSapiens" = org.Hs.eg.db)
    org <<- switch(isolate(input$strain),
                   "MMusculus" = "mouse",
                   "MAuratus" = "mouse",
                   "HSapiens" = "human")
    
    DEgene <<- res[res$padj<0.01,]
    Upgene <<- DEgene[DEgene$log2FoldChange>0,]
    Downgene <<- DEgene[DEgene$log2FoldChange<0,]
    
    #allgene.ensemmbl <- res$ensembl_gene_id
    DEgene.ensembl <- DEgene$ensembl_gene_id
    Upgene.ensembl <- Upgene$ensembl_gene_id
    Downgene.ensembl <- Downgene$ensembl_gene_id
    
    #allgene.entrez <- bitr(allgene.ensemmbl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = OrgDb)
    DEgene.entrez <<- bitr(DEgene.ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = OrgDb)
    
    DEgene_log2expression <<- DEgene$log2FoldChange
    names(DEgene_log2expression) <<- DEgene.entrez$ENTREZID
    
    kegg_res <- enrichKEGG(gene=DEgene.entrez$ENTREZID, organism=org, pAdjustMethod="BH", pvalueCutoff=0.01)
    output$kegg_res <- renderDataTable(as.data.frame(kegg_res))
    
    reactome_res <- ReactomePA::enrichPathway(gene=DEgene.entrez$ENTREZID, organism=org, pAdjustMethod="BH", pvalueCutoff=0.01)
    output$reactome_res <- renderDataTable(as.data.frame(reactome_res))
    
    output$reactome_dotplot <- renderPlot(dotplot(reactome_res, showCategory=15))
    output$reactome_emapplot <- renderPlot(emapplot(pairwise_termsim(reactome_res)))
    output$reactome_cnetplot <- renderPlot(cnetplot(reactome_res, categorySize="pvalue", foldChange=DEgene_log2expression, showCategory=8))
  })
  
  observeEvent(input$parameter_button, {
    plot1 <<- ggplot(res, aes(x=log2(baseMean), y=log2FoldChange, colour=factor(padj<0.01))) +
      geom_point(size=2, alpha=0.1) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(input$MA_xmin,input$MA_xmax), ylim=c(input$MA_ymin, input$MA_ymax))
    output$plot1 <- renderPlot(plot1)
    plot2 = ggplot(res, aes(x=log2(baseMean), y=log2FoldChange, colour=factor(external_gene_name==input$gene), alpha=factor(external_gene_name==input$gene))) +
      geom_point(size=2) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      scale_alpha_manual(values=c(0.05, 1)) +
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(input$MA_xmin,input$MA_xmax), ylim=c(input$MA_ymin, input$MA_ymax))
    output$plot2 = renderPlot(plot2)
    
    vol_plot1 <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), colour=factor(padj<0.01))) +
      geom_point(size=2, alpha=0.1) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(input$vol_xmin, input$vol_xmax), ylim=c(input$vol_ymin, input$vol_ymax))
    output$vol_plot1 <- renderPlot(vol_plot1)
    vol_plot2 <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj),colour=factor(external_gene_name==input$gene), alpha=factor(external_gene_name==input$gene))) +
      geom_point(size=2) +
      scale_color_manual(values=c("black", "red")) +
      scale_fill_manual(values = c("black", "red"), breaks=c("0", "1"), labels=c("non DEG", "DEG")) + 
      scale_alpha_manual(values=c(0.05, 1)) +
      theme(legend.title=element_blank()) + 
      coord_cartesian(xlim=c(input$vol_xmin, input$vol_xmax), ylim=c(input$vol_ymin, input$vol_ymax))
    output$vol_plot2 <- renderPlot(vol_plot2)
  })
  observeEvent(input$GSEA_button, {
    # GO enrichment analysis by clusteProfiler
    
    Upgene.symbol <- Upgene$external_gene_name
    Downgene.symbol <- Downgene$external_gene_name
    ont = switch(isolate(input$ontology),
                 "Biological Process" = "BP",
                 "Cellular Component" = "CC",
                 "Molecular Function" = "MF")
    
    go_up <- enrichGO(gene=Upgene.symbol, OrgDb=OrgDb, keyType="SYMBOL",ont=ont, pAdjustMethod="BH", pvalueCutoff=0.01)
    go_down <- enrichGO(gene=Downgene.symbol, OrgDb=OrgDb, keyType="SYMBOL",ont=ont, pAdjustMethod="BH", pvalueCutoff=0.01)
    output$go_up <- renderDataTable(as.data.frame(go_up))
    output$go_down <- renderDataTable(as.data.frame(go_down))
  })
  observeEvent(input$KEGGsearch_button, {
    pathview(gene.data = DEgene_log2expression, pathway.id = input$KEGG_ID, limit = list(gene=2, cpd=1), species=org)
    file.copy(from=paste0(input$KEGG_ID, ".pathview.png"), to=paste0("www/", input$KEGG_ID, ".pathview.png"))
    file.remove(paste0(input$KEGG_ID, ".pathview.png"))
    file.remove(paste0(input$KEGG_ID, ".png"))
    file.remove(paste0(input$KEGG_ID, ".xml"))
    output$KEGG_image <- renderUI({
      img(src=paste0(input$KEGG_ID, ".pathview.png"))
    })
  })
  session$onSessionEnded(function(){
    # stopApp()
    q("no")
  })
})