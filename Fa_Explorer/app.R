#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



# Load data ---------------------------------------------------------------
library(shiny)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,tidyverse,neuralnet,DescTools,reshape2)



# Define Functions and inputs: ---------------------------------------------------------
# needs as input a named vector with the code FA:CC:DB

fa2cl.predictor <- function(profile){
    #profile=fa.profile.result
    
    fa.cutoff <- 0.95
    fraction <- which(abs(cumsum(rev(profile[order(profile)]))-fa.cutoff)==min(abs(cumsum(rev(profile[order(profile)]))-fa.cutoff))) 
    fraction <- mean(fraction)
    fa <- names(tail(profile[order(profile)],n = fraction)) # Top 10 frequent cardiolipin fragments
    fa.profile <- tail(profile[order(profile)],n = fraction)
    fa.profile <- fa.profile/sum(fa.profile)
    
    #fa.comb <- arrange(expand.grid(a=fa,b=fa,c=fa,d=fa),a)
    fa.comb <- CombSet(fa, 4, repl = T, ord = T)
    fa.comb <- t(apply(X = fa.comb, 1, FUN = function(x) x <- x[order(x)]))
    
    cl.profile <- as.matrix(fa.comb) # Initiate cardiolipin profile object
    for(i in 1:length(fa)){ # Call fatty acid relative abuncance data and write into cardiolipin profile object
        cl.profile[cl.profile == fa[i]] <- fa.profile[i]
    }
    
    #x=fa.comb[1,]
    cl.cumprod <- apply(X = cl.profile, MARGIN = 1, FUN = function(x) prod(as.numeric(x)))
    cl.name <- apply(X = fa.comb, MARGIN = 1, FUN = function(x) paste(c("CL", paste("(", rownames(as.matrix(table(x))),")", as.matrix(table(x)),sep = "")),collapse = ""))
    cl.species <- apply(X = fa.comb, MARGIN = 1, FUN = function(x) paste(c("CL", colSums(matrix(as.numeric(unlist(strsplit(x, split = ":"))), ncol = 3, byrow = T)[,2:3])), collapse = ":"))
    
    result <- aggregate(cl.cumprod, by = list(cl.species), FUN = function(x) sum(x, na.rm = T))
    
    result.vector <- result[,2]
    names(result.vector) <- result[,1]
    return(result.vector)
    
}

input_FA_profiles <-  tibble::tribble(
    ~X16.0,     ~X16.1,    ~X18.0,    ~X18.1,    ~X18.2,     ~X18.3,      ~X20.1,     ~X20.2,     ~X20.3,     ~X20.4,      ~X20.5,      ~X22.1,     ~X22.2,     ~X22.4,      ~X22.5,     ~X22.6,
    0.2082069, 0.08864113, 0.1510188, 0.2323291, 0.1474864, 0.03866782, 0.005784702, 0.00447171, 0.02214417, 0.06654678, 0.001536865, 0.000835116, 0.00032633, 0.01044775, 0.002828304, 0.01872817
)

# actual app: -------------------------------------------------------------


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("FA_Explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel("User Input",
            sliderInput("FA18.1",
                        "FA 18:1 content",
                        min = 0,
                        max = 1,
                        value = 0.2323291),
            sliderInput("FA18.2",
                        "FA 18:2 content",
                        min = 0,
                        max = 1,
                        value = 0.1474864),
            sliderInput("FA20.4",
                        "FA 20:4 content",
                        min = 0,
                        max = 1,
                        value = 0.06654678),
            sliderInput("FA22.6",
                        "FA 22:6 content",
                        min = 0,
                        max = 1,
                        value = 0.01872817)
        ),
            

        # Show a plot of the generated distribution
        mainPanel(
            h3("Input Profile monomeric PL acyls available"),
           plotOutput("FA_Input",height = "200px"),
           h4(icon(name = "arrow-alt-circle-down"),align = "center"),
            h3("Output Profile: cardiolipin acyls incorporated"),
           plotOutput("FA_Output",height = "200px"),
           h4(icon(name = "arrow-alt-circle-down"),align = "center"),
            h3("Output Profile: cardiolipin profile"),
           plotOutput("CL_Profile",height = "200px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    load(file = "best_5_ANNs.RData")
    input_FA_profiles <-  tibble::tribble(
        ~X16.0,     ~X16.1,    ~X18.0,    ~X18.1,    ~X18.2,     ~X18.3,      ~X20.1,     ~X20.2,     ~X20.3,     ~X20.4,      ~X20.5,      ~X22.1,     ~X22.2,     ~X22.4,      ~X22.5,     ~X22.6,
        0.2082069, 0.08864113, 0.1510188, 0.2323291, 0.1474864, 0.03866782, 0.005784702, 0.00447171, 0.02214417, 0.06654678, 0.001536865, 0.000835116, 0.00032633, 0.01044775, 0.002828304, 0.01872817
    )
    dataInput <- reactive({
        input_FA_profiles[,4] <- input$FA18.1
        input_FA_profiles[,5] <- input$FA18.2
        input_FA_profiles[,10] <- input$FA20.4
        input_FA_profiles[,16] <- input$FA22.6
        input_FA_profiles %>% pivot_longer(cols = everything()) %>% mutate(value=value/sum(value)) %>% pivot_wider()
    })
    
    output$FA_Output <- renderPlot({
        # generate bins based on input$bins from ui.R
        data <- tibble(Names=c("FA 16:0",
                       "FA 16:1",
                       "FA 18:0",
                       "FA 18:1",
                       "FA 18:2",
                       "FA 18:3",
                       "FA 20:1",
                       "FA 20:2",
                       "FA 20:3",
                       "FA 20:4",
                       "FA 20:5",
                       "FA 22:1",
                       "FA 22:2",
                       "FA 22:4",
                       "FA 22:5",
                       "FA 22:6"
        ),Values=predict(best_5_ANNs[[4]],dataInput() , rep = 18)[1,])
        ggplot(data = data,aes(x  = Names,y = Values))+geom_col()+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust=.5))
    })
    output$FA_Input <- renderPlot({
        data <- dataInput() %>% pivot_longer(cols = everything())
        data$name <- c(
            "FA 16:0",
            "FA 16:1",
            "FA 18:0",
            "FA 18:1",
            "FA 18:2",
            "FA 18:3",
            "FA 20:1",
            "FA 20:2",
            "FA 20:3",
            "FA 20:4",
            "FA 20:5",
            "FA 22:1",
            "FA 22:2",
            "FA 22:4",
            "FA 22:5",
            "FA 22:6"
        )
        ggplot(data = data,aes(x  = name,y = value))+geom_col()+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust=.5))
    })
    output$CL_Profile <- renderPlot({
        data2 <- predict(best_5_ANNs[[4]],dataInput() , rep = 18)[1,]
        # data2 <- predict(best_5_ANNs[[4]],input_FA_profiles , rep = 18)[1,]
        names(data2) <- c(
            "FA:16:0",
            "FA:16:1",
            "FA:18:0",
            "FA:18:1",
            "FA:18:2",
            "FA:18:3",
            "FA:20:1",
            "FA:20:2",
            "FA:20:3",
            "FA:20:4",
            "FA:20:5",
            "FA:22:1",
            "FA:22:2",
            "FA:22:4",
            "FA:22:5",
            "FA:22:6"
        )
        data2 <- fa2cl.predictor(profile = data2) %>% enframe()
        ggplot(data = data2 %>% filter(value>0.005),aes(x  = name,y = value))+geom_col()+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust=.5))
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
