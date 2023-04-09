## UI -> aparência visual do app
## Controles
library(shiny)
library(plotly)
ui <- fluidPage(title = NULL,
                fluidRow(
                  column(width = 4,
                         column(width = 4,
                                ## Controle do slider do parâmetro A
                                sliderInput("A", label = "A:", ticks = F,
                                            min = 0, max = 0.01, value = 0.005, step = 0.0002)),
                         column(width = 4,
                                sliderInput("B", label = "B:", ticks = F,
                                            min = 0, max = 1, value = 0.5, step = 0.02)),
                         column(width = 4,
                                sliderInput("C", label = "C:", ticks = F,
                                            min = 0, max = 0.2, value = 0.1, step = 0.005))),
                  
                  column(width = 4,
                         column(width = 4,
                                sliderInput("D", label = "D:", ticks = F,
                                            min = 0, max = 0.01, value = 0.005, step = 0.0002)),
                         column(width = 4,
                                sliderInput("E", label = "E:", ticks = F,
                                            min = 0, max = 10, value = 5, step = 0.5)),
                         column(width = 4,
                                sliderInput("F_", label = "F:", ticks = F,
                                            min = 15, max = 110, value = 20, step = 1))),
                  
                  column(width = 4,
                         column(width = 4,
                                sliderInput("G", label = "G:", ticks = F,
                                            min = 0, max = 0.001, value = 0.0005, step = 0.00002)),
                         column(width = 4,
                                sliderInput("H", label = "H:", ticks = F,
                                            min = 1.05, max = 1.10, value = 1.075, step = 0.005)),
                         column(width = 4,
                                sliderInput("K", label = "K:", ticks = F,
                                            min = -1, max = 1, value = 0, step = 0.05)))
                ),
                fluidRow(plotlyOutput("hp_curve"))
)

## Função que gera a curva HP e cada componente
hp_curve <- function(x, p){
  a = p[1];b = p[2];c = p[3];d = p[4];e = p[5];f = p[6];g = p[7];h = p[8]
  k = p[9]
  c1 = a^((x+b)^c)
  c2 = d*exp(-e*log(x/f)^2)
  # c3 = g*h^x
  c3 = g*h^x/(1+k*g*h^x)
  curva = c1 + c2 + c3
  return(list(c1 = c1, c2 = c2, c3 = c3, curva = curva))
}

## Servidor -> onde são criados os objetos que aparecem na UI
server <- function(input, output){
  p = ggplot() +
    scale_y_continuous(trans = "log10", breaks = 10^-seq(0,4), limits = 10^-c(4,0), labels = scales::comma) +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = "Idade", y = "", title = NULL) + theme_bw()
  x = 0:100
  
  ## Objeto que vai ser plotado na UI
  output$hp_curve <- renderPlotly({
    curva = hp_curve(x = x,
                     c(input$A, input$B, input$C, input$D, input$E, input$F_, input$G, input$H, input$K
                     ))
    data = data.frame(Idade = x, m = curva$curva, c1 = curva$c1, c2 = curva$c2, c3 = curva$c3)
    p <- p + geom_line(data = data, aes(x = Idade, y = m, col = "Curva HP")) + 
      geom_line(data = data, aes(x = Idade, y = c1, col = "Bloco 1"), lty = 2) +
      geom_line(data = data, aes(x = Idade, y = c2, col = "Bloco 2"), lty = 2) +
      geom_line(data = data, aes(x = Idade, y = c3, col = "Bloco 3"), lty = 2) +
      guides(col = guide_legend(""))
    ggplotly(p)
  })
}

shinyApp(ui = ui, server = server)
