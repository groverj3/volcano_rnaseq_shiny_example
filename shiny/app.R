library(shiny)
library(DT)
library(tidyverse)
library(ggrepel)
library(DescTools)


# Volcano plot code based from https://github.com/groverj3/genomics_visualizations/blob/master/volcano_plotteR.r
volcplot <- function(data, padj_threshold = 0.05, fc = 1, plot_title = 'Volcano Plot', plot_subtitle = NULL) {

  # Set the fold-change thresholds
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)

  # Make a dataset for plotting, add the status as a new column
  plot_ready_data <- data %>%
    mutate_at('padj', ~replace(.x, is.na(.x), 1)) %>%
    mutate_at("padj", ~replace(.x, .x == 0, .Machine$double.xmin)) %>%  # When p values are zero, they're actually below the lowest value R can display
    mutate_at('log2FoldChange', ~replace(.x, is.na(.x), 0)) %>%
    mutate(
      log2fc_threshold = ifelse(log2FoldChange >= pos_log2fc & padj <= padj_threshold, 'up',
                         ifelse(log2FoldChange <= neg_log2fc & padj <= padj_threshold, 'down', 'ns')
        )
    )

  # Get the number of up, down, and unchanged genes
  up_genes <- plot_ready_data %>% filter(log2fc_threshold == 'up') %>% nrow()
  down_genes <- plot_ready_data %>% filter(log2fc_threshold == 'down') %>% nrow()
  unchanged_genes <- plot_ready_data %>% filter(log2fc_threshold == 'ns') %>% nrow()

  # Make the labels for the legend
  legend_labels <- c(
      str_c('Up: ', up_genes),
      str_c('NS: ', unchanged_genes),
      str_c('Down: ', down_genes)
  )

  # Set the x axis limits, rounded to the next even number
  x_axis_limits <- DescTools::RoundTo(
    max(abs(plot_ready_data$log2FoldChange)),
    2,
    ceiling
  )

  # Set the plot colors
  plot_colors <- c(
      'up' = 'firebrick1',
      'ns' = 'gray',
      'down' = 'dodgerblue1'
  )


  # Make the plot, these options are a reasonable strting point
  plot <- ggplot(plot_ready_data) +
    geom_point(
      alpha = 0.25,
      size = 1.5
    ) +
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = log2fc_threshold,
      label = gene_id
    ) +
    geom_vline(
      xintercept = c(neg_log2fc, pos_log2fc),
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = 'dashed'
    ) +
    scale_x_continuous(
      'log2(FC)',
      limits = c(-x_axis_limits, x_axis_limits)
    ) +
    scale_color_manual(
      values = plot_colors,
      labels = legend_labels
      ) +
    labs(
      color = str_c(fc, '-fold, padj ≤', padj_threshold),
      title = plot_title,
      subtitle = plot_subtitle
    ) +
    theme_bw(base_size = 24) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(color = 'black'),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),  # Reduces dead area around legend
      legend.spacing.x = unit(0.2, 'cm')
    )
  plot
}


# Define UI
ui <- page_sidebar(
  title = "Volcano PlotteR",
  sidebar = sidebar(
    fileInput("deseq2_results", "DESeq2 Results Table"),
    numericInput("foldchange_threshold", "Fold Change Threshold", value = 1),
    numericInput("padj_threshold", "Adjusted p-value Threshold", value = 0.1),
    textInput("plot_title", "Plot Title", value = "Volcano Plot"),
    textInput("plot_subtitle", "Plot Subtitle", value = NULL),
  ),
  card(
    plotOutput("volcano_plot"),
    min_height = 580  # Ensures you don't have to scroll within this card
  ),
  DTOutput("deseq2_table"),
)


# Server function
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)

  deseq2_results <- reactive({
    req(input$deseq2_results)
    read_csv(input$deseq2_results$datapath)
  })

  deseq2_results_filtered <- reactive({
    req(deseq2_results)
    deseq2_results() %>%
      filter(2^abs(log2FoldChange) >= input$foldchange_threshold & padj <= input$padj_threshold)
  })

  output$deseq2_table <- renderDT(
    deseq2_results_filtered()
  )

  output$volcano_plot <- renderPlot({
    req(deseq2_results)
    deseq2_results() %>%
      volcplot(
        padj_threshold = input$padj_threshold,
        fc = input$foldchange_threshold,
        plot_title = input$plot_title,
        plot_subtitle = input$plot_subtitle
      )
    },
    height = 550
  )
}

# Run the application
shinyApp(ui = ui, server = server)