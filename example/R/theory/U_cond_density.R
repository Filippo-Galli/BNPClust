library(shiny)

# Define the conditional density function
log_density_U <- function(u, K, n, a, sigma, tau = 1) {

  # Term 1: u^(n-1) -> (n-1)*log(u)
  term1 <- log(u) * (n - 1)
  
  # Term 2: 1/(u + tau)^(n - sigma*K) -> (n - sigma*K)*log(u + tau)
  term2 <- (n - sigma * K) * log(u + tau)
  
  # Term 3: exp(-(a/sigma) * ((u + tau)^sigma - tau^sigma))
  a_over_sigma <- a / sigma
  tau_power_sigma <- tau^sigma
  term3 <- -a_over_sigma * ((u + tau)^sigma - tau_power_sigma)
  
  # Return density
  term1 - term2 + term3
}

ui <- fluidPage(
  titlePanel("Conditional log-Density U(u | ·)"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("K", 
                  "K (number of categories):", 
                  min = 1, 
                  max = 40, 
                  value = 25, 
                  step = 1),
      
      sliderInput("n", 
                  "n (sample size):", 
                  min = 10, 
                  max = 200, 
                  value = 100, 
                  step = 10),
      
      sliderInput("a", 
                  "a (parameter):", 
                  min = 0.01, 
                  max = 5, 
                  value = 1, 
                  step = 0.1),
      
      sliderInput("sigma", 
                  "σ (sigma):", 
                  min = 0.01, 
                  max = 0.3, 
                  value = 0.1, 
                  step = 0.01),
      
      sliderInput("tau", 
                  "τ (tau):", 
                  min = 0, 
                  max = 10, 
                  value = 1, 
                  step = 0.1),
      
      hr(),
      
      sliderInput("v_range", 
                  "v range:", 
                  min = 0.01, 
                  max = 2000, 
                  value = c(0.01, 2000), 
                  step = 0.5)
    ),
    
    mainPanel(
      plotOutput("densityPlot", height = "500px"),
      hr(),
      verbatimTextOutput("paramInfo")
    )
  )
)

server <- function(input, output) {
  
  output$densityPlot <- renderPlot({
    # Generate v values
    u <- seq(input$v_range[1], input$v_range[2], length.out = 500)
    
    # Compute density for each v
    density_values <- sapply(u, function(u_i) {
      log_density_U(u_i, 
                           K = input$K, 
                           n = input$n, 
                           a = input$a, 
                           sigma = input$sigma, 
                           tau = input$tau)
    })
    
    # Handle numerical issues (Inf, NA)
    valid_idx <- is.finite(density_values) & density_values > 0
    
    if (sum(valid_idx) > 0) {
      # Plot
      plot(u[valid_idx], density_values[valid_idx], 
           type = "l", 
           lwd = 2, 
           col = "steelblue",
           xlab = "v",
           ylab = "Density",
           las = 1)
      
      grid()
      
      # Add vertical line at mode (approximate)
      mode_idx <- which.max(density_values[valid_idx])
      mode_u <- u[valid_idx][mode_idx]
      abline(v = mode_u, col = "red", lty = 2, lwd = 1.5)
      legend("topright", 
             legend = paste("Mode ≈", round(mode_u, 3)),
             col = "red", lty = 2, lwd = 1.5, bty = "n")
    } else {
      plot(0, 0, type = "n", 
           xlab = "v", ylab = "Density",
           main = "No valid density values (numerical overflow/underflow)")
      text(0, 0, "Try different parameter values", col = "red", cex = 1.5)
    }
  })
  
  output$paramInfo <- renderText({
    paste(
      "Current Parameters:",
      paste("K =", input$K),
      paste("n =", input$n),
      paste("a =", input$a),
      paste("σ =", input$sigma),
      paste("τ =", input$tau),
      paste("v range = [", input$v_range[1], ",", input$v_range[2], "]"),
      sep = "\n"
    )
  })
}

shinyApp(ui = ui, server = server)