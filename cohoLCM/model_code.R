
##########################################################################
################################ Packages ################################
##########################################################################

library(shiny)
library(ggplot2)
library(plyr)
library(scales)

##########################################################################
############################### Functions ################################
##########################################################################

## load functions (Rep, Sum, Win, Cap, mmsR, mmsW, msm) ##
beta.mom<-function(mean,sd){
  v<-sd**2
  x<-mean
  a<-x*(x*(1-x)/v-1)
  b<-(1-x)*(x*(1-x)/v-1)
  c(a,b)
}

Rep <- function(spawner, fecund, sex_ratio, interceptRep, Kegg, determ){
  females <- qbinom(0.5,spawner,sex_ratio)
  eggs <- females*fecund
  eggS <- 1/((1/interceptRep)+(eggs/Kegg)) # Beverton-Holt
  if(determ==0){
    fry <- rbinom(1, eggs, eggS) # fish surviving to fry given number of eggs and survival
  } else {
    # if deterministic grab the 50th percentile (mean) value
    fry <- qbinom(0.5, eggs, eggS) 
  }
  return(fry)
}

Sum <- function(fry, rescue_level, facility_cap, encounter_prob, interceptSum, Ksum, determ){
  if(rescue_level<=1){ # if rescue level is a proportion...
    fryWild <- fry-min(c(qbinom(0.5, fry, rescue_level), facility_cap))
  } else { # if rescue level is a specific quantity...
    if(fry >= min(rescue_level, facility_cap)/encounter_prob){ # if there are enough fry to rescue...
      fryWild <- fry-min(c(rescue_level, facility_cap))
    } else { # when there are too few fish to rescue
      fryWild <- qbinom(0.5, fry, (1-encounter_prob))
    }
  }
  
  if(determ==0){ # if stochastic...
    fryWildS <- 1/((1/interceptSum)+(fryWild/Ksum))
    parr <- rbinom(1, fryWild, fryWildS) # wild fish surviving to parr given number of wild fry and survival
  } else {  # if deterministic...
    fryWildS <- 1/((1/interceptSum)+(fryWild/Ksum))
    # if deterministic grab the 50th percentile value
    parr <- qbinom(0.5, fryWild, fryWildS) 
  }
  return(c(fryWild, parr))
}

WinST <- function(parrW, parrR, interceptWin, Kwin, determ){
  parrS <- 1/((1/interceptWin)+((parrW+parrR)/Kwin))
  if(determ==0){ # if stochastic...
    smoltW <- rbinom(1, parrW, parrS) # wild fish surviving to smolt given number of parr and survival
    smoltR <- rbinom(1, parrR, parrS) # rescued fish surviving to smolt given number of parr and survival
  } else {  # if deterministic
    # if deterministic grab the 50th percentile (mean) value
    smoltW<- qbinom(0.5, parrW, parrS)
    smoltR<- qbinom(0.5, parrR, parrS) 
  }
  return(c(smoltW, smoltR))
}

Win <- function(parr, interceptWin, Kwin, determ){
  parrWildS <- 1/((1/interceptWin)+(parr/Kwin))
  if(determ==0){ # if stochastic...
    smoltWild <- rbinom(1, parr, parrWildS) # wild fish surviving to smolt given number of parr and survival
  } else {  # if deterministic
    # if deterministic grab the 50th percentile (mean) value
    smoltWild <- qbinom(0.5, parr, parrWildS) 
  }
  return(smoltWild)
}

CapST <- function(fry, rescue_level, facility_cap, encounter_prob, interceptCap, determ){
  if(rescue_level<=1){  # rescue level is a proportion (0-1)
    fryRes <- min(c(qbinom(0.5, fry, rescue_level), facility_cap))
  } else {              # rescue level is a specific number (>1 fish rescued)
    if(fry >= min(rescue_level, facility_cap)/encounter_prob){ # if fry abundance is high enough...
      fryRes <- min(rescue_level, facility_cap)  # rescue amount specified up to the facility capacity
    } else {  # if fry abundance is too low...
      fryRes <- qbinom(0.5, fry, encounter_prob)  # rescue a proportion of available fry (attempt to avoid rescuing all fry since not realistic; cannot locate all fry when abundance is low)
    }
  }
  
  if(determ==0){ # if stochastic
    parrR <- rbinom(1, fryRes, interceptCap) # rescued fish surviving to smolt given number of rescued fry and survival
  } else {
    parrR <- qbinom(0.5, fryRes, interceptCap)
  }
  return(c(fryRes, parrR))
}

Cap <- function(fry, rescue_level, facility_cap, encounter_prob, interceptCap, determ){
  if(rescue_level<=1){  # rescue level is a proportion (0-1)
    fryRes <- min(c(qbinom(0.5, fry, rescue_level), facility_cap))
  } else {              # rescue level is a specific number (>1 fish rescued)
    if(fry >= min(rescue_level, facility_cap)/encounter_prob){ # if fry abundance is high enough...
      fryRes <- min(c(rescue_level, facility_cap))  # rescue  amount specified
    } else {  # if fry abundance is lower than the lowest value of rescue level and facility capacity...
      fryRes <- qbinom(0.5, fry, encounter_prob)  # rescue a proportion of available fry (attempt to avoid rescuing all fry since not realistic; cannot locate all fry when abundance is low; ** needs to be reworked **)
    }
  }
  
  if(determ==0){ # if stochastic
    smoltR <- rbinom(1, fryRes, interceptCap) # rescued fish surviving to smolt given number of rescued fry and survival
  } else {
    smoltR <- qbinom(0.5, fryRes, interceptCap)
  }
  return(c(fryRes, smoltR))
}

mmsW <- function(smoltW, alpha_mmsW, beta_mmsW, determ){
  if(determ==0){
    adultWS <- rbeta(1, alpha_mmsW, beta_mmsW) # wild fish survival
    adultW <- rbinom(1, smoltW, adultWS)       # wild fish surviving to adult given number of wild smolt and survival probability
  } else {
    # if deterministic grab the 50th percentile value
    adultWS <- qbeta(0.5, alpha_mmsW, beta_mmsW)
    adultW <- qbinom(0.5, smoltW, adultWS) 
  }
  return(adultW)
}

mmsR <- function(smoltR, alpha_mmsR, beta_mmsR, determ){
  if(determ==0){
    adultRS <- rbeta(1, alpha_mmsR, beta_mmsR) # rescued fish survival
    adultR <- rbinom(1, smoltR, adultRS) # rescued fish surviving to adult given number of rescued smolt and survival probability
  } else {
    # if deterministic grab the 50th percentile value
    adultRS <- qbeta(0.5, alpha_mmsR, beta_mmsR) 
    adultR <- qbinom(0.5, smoltR, adultRS) 
  }
  return(adultR)
}

msm <- function(adultW, adultR, alpha_msm, beta_msm, determ){
  adults <- adultW + adultR
  if(determ==0){
    msmS <- rbeta(1, alpha_msm, beta_msm) # late marine survival
    spawner <- rbinom(1, adults, msmS) # total fish surviving to spawner given number of total adults and survival probability
  } else {
    # if deterministic grab the 50th percentile (mean) value
    msmS <- qbeta(0.5, alpha_msm, beta_msm)
    spawner <- qbinom(0.5, adults, msmS) 
  }
  return(spawner)
}

safeSample <- function(x, reps) if(length(x) == 1) rep(x, reps) else sample(x, reps, replace = T)  # fixes sampling error associated with length = 1

extinctRisk <- function(data,threshold,gens){
  risk <- sum(data<threshold)/gens
  return(risk)
}

##########################################################################
############################# Static Inputs ##############################
##########################################################################

## model option for stocastic = 0 or deterministic = 1
determ <- 0

##########################################################################
##########################################################################
############################ User Interface ##############################
##########################################################################
##########################################################################

ui <- fluidPage(
  titlePanel("Effects of Fish Rescue on Coho Salmon Population Dynamics"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("submit", label = "Run Simulation", style="float: right;"),
      h4(strong("Model Settings")),
      sliderInput("reps",
                  label = "Number of simulations:",
                  min = 1000, max = 10000, value = 3000, step = 500),
      sliderInput("gens",
                  label = "Simulation length: (in generations)",
                  min = 10, max = 50, value = 33, step = 1),
      numericInput("initial",
                   label = "Initial spawner abundance:",
                   min = 4, value = 300),
      numericInput("extn_thresh",
                  label = "Extinction threshold:",
                  min = 0, value = 50),
      br(),
      
      ## Rescue Parameters
      h4(strong("Rescue Parameters")),
      selectInput("Rlength",
                  label = "Holding time",
                  choices = list("Short-term" = 0, "Long-term" = 1), selected = 1),
      numericInput("rescue_level",
                  label = "Number of fry rescued:",
                  min = 0,value = 10000),
      sliderInput("encounter_prob",
                  label = "Encounter Probability:",
                  min = 0.2, max = 1, value = 0.4, step = 0.05), 
      sliderInput("interceptCap",
                  label = "Captivity survival rate:",
                  min = 0.2, max = 1, value = 0.98),
      sliderInput("mmsPenalty",
                  label = "Marine survival penalty:",
                  min = 0.2, max =1, value = .6, step = .1),
      br(),
      
      ## Freshwater Rearing Capacities
      h4(strong("Freshwater Rearing Capacities")),
      numericInput("Kegg",
                   label = "Spawning capacity (# eggs):",
                   min = 5000, value = 902500),
      numericInput("Ksum",
                   label = "Summer capacity (# fry):",
                   min = 5000, value = 10000),
      numericInput("Kwin",
                   label = "Winter capacity (# fry):",
                   min = 5000, value = 20000),
      br(),
      
      ## Reproduction
      h4(strong("Reproduction")),
      sliderInput("fecund",
                  label = "Fecundity: (eggs per female)",
                  min = 2000, max = 3000, value = 2500),
      sliderInput("sex_ratio",
                  label = "Spawner sex ratio: (percent female)",
                  min = 0.25, max = 0.75, value = 0.5),
      br(),
      
      ## Density Independent Survival
      h4(strong("Density Independent Survival (Max Productivity)")),
      sliderInput("interceptRep",
                  label = "Egg to fry survival:",
                  min = 0.05, max = 0.75, value = 0.429),
      sliderInput("interceptSum",
                  label = "Oversummer survival:",
                  min = 0.05, max = 0.85, value = 0.26),
      sliderInput("interceptWin",
                  label = "Overwinter survival:",
                  min = 0.300, max = 1, value = 0.9),
      br(),
      
      ## Early Marine Survival for Wild Fish
      h4(strong("Early Marine Survival for Wild Fish")),
      sliderInput("meanS_mmsW",
                  label = "Mean:",
                  min = 0, max = 0.4, value = 0.098),
      sliderInput("sdS_mmsW",
                  label = "Standard deviation:",
                  min = 0, max = 0.1, value = 0.01),
      br(),
      
      # Late Marine Survival for All Fish
      h4(strong("Late Marine Survival for All Fish")),
      sliderInput("meanS_msm",
                  label = "Mean:",
                  min = 0, max = 0.4, value = 0.15),
      sliderInput("sdS_msm",
                  label = "Standard deviation:",
                  min = 0, max = 0.1, value = 0.01)
      ),
    mainPanel(
      tabsetPanel(
        tabPanel("Introduction",
                 h3("Welcome to our R Shiny Application!"),
                 br(),
                 p("In the Mediterranean climate of the Pacific Northwest, dry summers result in prolonged periods of low flows in streams. Climate change, habitat alterations, and increasing water demands are leaving less water available for streams. As water levels drop, some small streams become fragmented, transforming from a ribbon of continuous habitat into a series of isolated pools. Fragmented streams pose a serious threat to salmon; juveniles that become stranded in small pools are at increased risk to overheat, starve, or be consumed by predators."),
                 p(div(img(src="stream2.png", width = 408, height = 306), img(src="stream.png", width = 408, height = 306), style="text-align: center;")),
                 p("Healthy salmon populations can cope with fragmentation and recover from a bad drought-year.  However, for depressed populations (e.g., ESA-listed), there is a strong need for adaptation strategies that reduce the risks of drought-induced fragmentation. While a loud and high profile debate has ensued over the permanent translocation of species (i.e., assisted migration), there has been much less attention given to the use of temporary translocation to buffer species from seasonal risks.  One form of temporary translocation is known as fish rescue, which seeks to reduce drought related mortality in wild fish by manually moving individuals from fragmented areas to artificial rearing facilities that provide refuge during periods of low flow.  Fish rescue differs from the more typical fish salvage operations where individuals are instead moved to free-flowing habitat."),
                 br(),
                 p("This Shiny application explores how seasonal fish rescue affects the population dynamics of coho salmon. This application considers fish rescue that involves seasonal or annual translocation to captive rearing. Before using this application, we strongly recommend that you first read Brittany Beebe's 2019 Master's Thesis from Oregon State University entitled", em("Evaluating Fish Rescue as a Drought Adaptation Strategy for Imperiled Coho Salmon: A Life-Cycle Modeling Approach"), tags$a(href = "https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/73666b76b?locale=en", "(Link)."), "This document provides a full description of the life-cycle model and additional background on fish rescue and ecological drought. The thesis document also provides discussion of the potential costs and benefits associated with fish rescue.")
        ),
        
        tabPanel("Model Methods",
                 h3("Understanding Model Mechanics"),
                 p("The model framework reflects the basic life cycle of coho salmon and simulates abundance across five serial life stages: spawner, fry, parr, smolt, and adult. Wild fish move sequentially through each life stage and are subjected to corresponding survival processes. Fish rescue is simulated as an alternative pathway between the fry and smolt life stages. The model explores two captive rearing durations: long- and short-term rescue. For long-term rescue, fish are held for roughly a year through the fry-to-smolt stages. This transition encompasses both summer and winter survival. For the short-term rescue, fish are held in captivity only during the season of wetted habitat contraction, which was considered as the summer fry-to-parr survival stage."),
                 p(div(img(src="conceptual_model_Page_1.jpg", width = 605, height = 467.5), img(src="conceptual_model_Page_2.jpg", width = 605, height = 467.5), style="text-align: center;")),
                 h4("Survival"),
                 p("Survival is represented by the transition of individuals from one life stage to the next. The transition between the spawner and fry life stages includes fecundity in addition to survival. Three general methods were used to calculate survival. The Beverton-Holt function is used for freshwater survival (Moussalli & Hilborn 1986), a beta distribution is used for marine survival (Hill et al. 2003), and a constant parameter value is used for survival in captivity (NWFR, personal communication). The survival estimates calculated for each transition are used in a binomial distribution to incorporate demographic stochasticity to the simulated number of individuals at the subsequent life stage (Nickelson & Lawson 1998)."),
                 br(),
                 h4("Rescue"),
                 p("Although not bred in captivity, rescued fish spend much of their freshwater life stages in artificial conditions, and we assumed this could negatively affect smolt-to-adult return rates. Reported estimates of marine survival of hatchery fish compared to wild fish range from ~37-100% (Jonsson et al. 2003; Kallio-Nyberg et al. 2004; Jokikokko et al. 2006; Hyvarinen & Rodewald 2013). In the model, early marine survival of rescued fish is calculated by multiplying early marine survival of wild fish by a marine survival penalty value."),
                 p("The attempted rescue level remains constant for every generation of a given simulation. Since fry abundance may fall below the attempted rescue level, a capped linear functional response (i.e., hockey-stick) is used in which a percentage of the total fry are captured up to the attempted rescue level. In other words, when fry abundance is less than the captive rearing capacity, an encounter probability is used. This avoids rescuing an unrealistic percentage of available fry."),
                 br(),
                 h4("Drought"),
                 p("Drought conditions are explored by adjusting the summer rearing capacity. To explore more severe drought, set the summer capacity to a lower value.")
        ),
        
        tabPanel("Scenario Plots",
                 h4("Mean Smolt Abundance Over Time"),
                 plotOutput("smolt_timeseries_plot"),
                 br(),
                 h4("Mean Spawner Abundance Over Time"),
                 plotOutput("spawner_timeseries_plot"),
                 br(),
                 h4("Extinction Risk Over Time"),
                 plotOutput("extinction_timeseries_plot")#,
                 #br(),
                 #h4("Head of Time Data"),
                 #tableOutput("timeTable")
                 ),
        
        tabPanel("Comparative Plots",
                 h4("Mean Spawner Abundance Across Drought Conditions"),
                 plotOutput("point_abund_by_drought"),
                 br(),
                 h4("Extinction Risk Across Drought Conditions"),
                 plotOutput("point_extinct_by_drought")
                 ),
        
        tabPanel("Tables",
                 h4("Parameter Values"),
                 tableOutput("selected_parameters"),
                 br(),
                 h4("Simulation Results"),
                 downloadButton("downloadResults", "Download All Results"),
                 tableOutput("table")#,
                 #br(),
                 #h4("Data"),
                 #tableOutput("table2")
                 ),
        
        tabPanel("Metadata",
                 h3("Model Settings"),  ## Model Settings
                 br(),
                 p(strong("- Number of simulations:"), "The number of simulations to run. Default = 3,000 simulations"),
                 p(strong("- Simulation length:"), "The length of time each simulation is run. Default = 33 generations (i.e., 99 years)"),
                 p(strong("- Initial spawner abundance:"), "The initial abundance of total spawners (male plus female). Default = 300 spawners"),
                 p(strong("- Extinction threshold:"), "A value of spawner abundance (male plus female) below which extinction is likely. Default = 50 spawners"),
                 br(),
                 h3("Rescue Parameters"),  ## Rescue Parameters
                 p(strong("- Holding time:"), "The length of time that rescued fish are held in captivity. Long-term implies rescued fish are held over both summer and winter. Short-term implies rescued fish are held only over summer.  Default = Long-term"),
                 p(strong("- Number of fry rescued:"), "The number of fry rescued for each rescue period. Default = 10,000"),
                 p(strong("- Encounter probability:"), "The probability of encountering fry when abundance is lower than number of fry to be rescued. This avoids rescuing all fish when fry abundance is low. Default = 0.4"),
                 p(strong("- Captivity survival rate:"), "Rescued fish survival spanning from the time of stream removal in early summer until the time of release in fall for short-term rescue or the following spring for long-term rescue. Default = 0.98 (Source: personal communication with Northwest Wild Fish Rescue)"),
                 p(strong("- Marine survival penalty:"), "Penalty applied to early marine survival for rescued fish. The value selected represents survival of rescued fish compared to survival of wild fish. A lower value indicates lower survival and thus a greater penalty for being rescued. For example, a value of 40% indicates that rescued fish survival is 40% of the wild fish survival (i.e., rescued fish survival = 40% of wild fish survival). Default = 60%"),
                 br(),
                 h3("Freshwater Rearing Capacities"),  ## Freshwater Rearing Capacities
                 p(strong("- Spawning capacity:"), "The maximum number of eggs sustained by the stream. Calculated from multiplying stream length (km) by maximum spawner density (19 females/km). Default = 902,500 (Source: Bradford et al. 2000)" ),
                 p(strong("- Summer capacity:"), "The maximum number of fry sustained by the stream when in its most contracted state during the summer.", span("This parameter is used as a proxy for drought condition.", style = "color:blue"), "A low number indicates more severe drought, while a larger number indicates less severe drought. Default = 10,000"),
                 p(strong("- Winter capacity:"), "The maximum number of fry sustained by the stream during highest winter flows. Default = 20,000"),
                 br(),
                 h3("Reproduction"),  ## Reproduction
                 p(strong("- Fecundity:"), "The number of eggs produced by each female spawner. Default = 2500 eggs/female"),
                 p(strong("- Spawner sex ratio:"), "The percent of total spawners that are female. Default = 50% female (i.e., 1:1)"),
                 br(),
                 h3("Density Independent Survival (Max Productivity)"),  ## Density Independent Survival (Max Productivity)
                 p(strong("- Egg to fry survival:"), "Survival from egg to fry as density approaches zero. Default = 0.429 (Source: Nickelson 1998)"),
                 p(strong("- Oversummer survival:"), "Oversummer survival of wild fish from fry to parr as density approaches zero. Default = 0.26 (Source: Nickelson 1998)"),
                 p(strong("- Overwinter survival:"), "Overwinter survival of wild fish from parr to smolt as density approaches zero. Default = 0.9 (Source: Nickelson 1998)"),
                 br(),
                 h3("Early Marine Survival for Wild Fish"),  ## Early Marine Survival for Wild Fish
                 p(strong("- Mean:"), "Mean marine survival for wild fish from outmigration through first year in marine environment. Default = 0.098 (Source: Bradford 1995)"),
                 p(strong("- Standard deviation:"), "Standard deviation for early marine survival for wild fish. Default = 0.01"),
                 br(),
                 h3("Late Marine Survival for All Fish"),  ## Late Marine Survival for All Fish
                 p(strong("- Mean:"), "Mean late marine survival for all fish (wild and rescued) Default = 0.15"),
                 p(strong("- Standard deviation:"), "Standard deviation for late marine survival for all fish. Default = 0.01")
                 ),
        
        tabPanel("References",
                 p(strong("Barrowman"), "NJ, RA Myers, R Hilborn, DG Kehler, and CA Field. 2003. The variability among populatinos of coho salmon in the maximum reproductve rate and depensation. Ecological Applications 13:784-793."),
                 p(strong("Beebe."), "2019. Evaluating Fish Rescue as a Drought Adaptation Strategy for Imperiled Coho Salmon: A Life-Cycle Modeling Approach. Master's Thesis. Oregon State University, Corvallis, Oregon.", tags$a(href = "https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/73666b76b?locale=en", "Link")),
                 p(strong("Bradford."), "1995. Comparative review of Pacific salmon survival rates. Canadian Journal of Fisheries and Aquatic Sciences 52:1327-1338."),
                 p(strong("Bradford"), "MJ, RA Myers, and JR Irvine. 2000. Reference points for coho salmon (Oncorhynchus kisutch) harvest rates and escapement goals based on freshwater production. Canadian Journal of Fisheries and Aquatic Sciences 57:677-686."),
                 p(strong("Hill"), "MF, LW Botsford, and A Hastings. 2003. The effects of spawning age distribution on salmon persistence in fluctuating environments. Journal of Animal Ecology 72:736-744."),
                 p(strong("Hyvarinen"), "P and P Rodewald. 2013. Enriched rearing improves survival of hatchery-reared Atlantic salmon smolts during migration in the River Tornionjoki. Canadian Journal of Fisheries and Aquatic Sciences 70:1386-1395."),
                 p(strong("Jokikokko"), "E, I Kallio-Nyberg, I Saloniemi, and E Jutila. 2006. The survival of semi-wild, wild and hatchery-reared Atlantic salmon smolts of the Simojoki River in the Baltic Sea. Journal of Fish Biology 68:430-442."),
                 p(strong("Jonsson"), "N, B Jonsson, and LP Hansen. 2003. The marine survival and growth of wild and hatchery-reared Atlantic salmon: growth and survival of salmon. Journal of Applied Ecology 40:900-911."),
                 p(strong("Kallio-Nyberg"), "I, E Jutila, I Saloniemi, and E Jokikokko. 2004. Association between environmental factors, smolt size and the survival of wild and reared Atlantic salmon from the Simojoki River in the Baltic Sea. Journal of Fish Biology 65:122-134."),
                 p(strong("Moussalli"), "E & R Hilborn. 1986. Optimal Stock Size and Harvest Rate in Multistage Life History Models. Canadian Journal of Fisheries and Aquatic Sciences 43:135-141."),
                 p(strong("Nickelson"), "TE. 1998. A habitat-based assessment of coho salmon production potential and spawner escapement needs for Oregon coastal streams. Oregon Department of Fish Wildlife Fish Division Information Report 98-4."),
                 p(strong("Nickelson"), "TE and PW Lawson. 1998. Population viability of coho salmon, Oncorhynchus kisutch, in Oregon coastal basins: application of a habitat-based life cycle model. Canadian Journal of Fisheries and Aquatic Sciences 55:2383-2392."),
                 p(strong("Steel"), "EA, DW Jensen, KM Burnett, K Christiansen, JC Firman, BE Feist, KJ Alauf, and DP Larson. 2012. Landscape characteristics and coho salmon (Oncorhynchus kisutch) distributions: explaining abundance versus occupancy. Canadian Journal of Fisheries and Aquatic Sciences 69:457-468.")
        )
        )
      )
    )
  )


##########################################################################
##########################################################################
################################# SERVER #################################
##########################################################################
##########################################################################

server <- function(input,output){
  
  # Function that runs model with updated values
  modelData <- reactive({
  
  ##########################################################################
  ########################### Initial conditions ###########################
  ##########################################################################
  
  ## simulations/replicates
  reps <- input$reps           # number of simulations/replicates to run
  
  ## projection timeframe
  gens <- input$gens           # length of time to simulate (in generations)
  
  ## rescue length
  Rlength <- input$Rlength     # model option: short-term rescue = 0 or long-term rescue = 1
  
  ## starting spawner abundance
  initial <- input$initial     # includes males and females
  
  ## extinction threshold (total spawners)
  extn_thresh <- input$extn_thresh

  
  ##########################################################################
  ############################### Parameters ###############################
  ##########################################################################
  
  ### Rep (spawner to total fry) ###
  
  fecund <- input$fecund                  # fecundity
  Kegg <- input$Kegg                      # Egg capacity, calculated
  sex_ratio <- input$sex_ratio            # Percent of spawners that are female; assumes a 1:1 sex ratio
  interceptRep <- input$interceptRep      # Survival as density approaches 0 (adjustable parm calculated from Nickelson 1998)
  
  ##########################################################################
  ##########################################################################
  
  ### Sum (wild fry to wild parr) ###
  
  rescue_level <- input$rescue_level      # number of fry rescued
  mmsPenalty <- input$mmsPenalty          # marine survival penalty for rescued fish
  
  interceptSum <- input$interceptSum      # survival as density approaches 0 (adjustable parm calculated from Nickelson 1998)
  Ksum <- input$Ksum                      # late summer habitat capacity
  encounter_prob <- input$encounter_prob  # probability of encountering fish when low fry abundance; avoids catching all fish when abundance is low and rescue level is high
  
  ##########################################################################
  ##########################################################################
  
  ### Win (wild parr to wild smolt) ###
  
  interceptWin <- input$interceptWin  # survival as density approaches 0 (Nickelson 1998)
  Kwin <- input$Kwin                  # winter habitat capacity
  
  ##########################################################################
  ##########################################################################
  
  ### Cap (rescued fry to rescued smolt) ###
  
  interceptCap <- input$interceptCap     # survival as density approaches 0
  facility_cap <- rescue_level*1.3       # how many fish the facility can hold with near 100% survival; set at max rescue level (unresticted capacity)
  
  ##########################################################################
  ##########################################################################
  
  ### mmsR (rescued smolt to rescued adult) ###
  # *commented out because calculated within simulations
  
  #meanS_mmsR <- varies w/penalty  # mean early marine survival for rescued fish 
  sdS_mmsR <- input$sdS_mmsW       # sd for early marine survival for rescued fish; set equal to that of wild fish
  
  #alpha_mmsR       # alpha parameter calculated from meanS_mmsR and sdS_mmsR of rescued early marine survival using beta.mom function
  #beta_mmsR        # beta parameter calculated from meanS_mmsR and sdS_mmsR of rescued early marine survival using beta.mom function
  
  #alpha_mmsR <- beta.mom(meanS_mmsR, sdS_mmsR)[1]
  #beta_mmsR <- beta.mom(meanS_mmsR, sdS_mmsR)[2]
  
  ##########################################################################
  ##########################################################################
  
  ### mmsW (wild smolt to wild adult) ###
  # *commented out because calculated within simulations
  
  meanS_mmsW <- input$meanS_mmsW   # mean early marine survival for wild fish
  sdS_mmsW <- input$sdS_mmsW       # sd for early marine survival for wild fish
  
  #alpha_mmsW        # alpha parameter calculated from meanS_mmsW and sdS_mmsW of wild early marine survival using beta.mom function
  #beta_mmsW         # beta parameter calculated from meanS_mmsW and sdS_mmsW of wild early marine survival using beta.mom function
  
  #alpha_mmsW <- beta.mom(meanS_mmsW, sdS_mmsW)[1]
  #beta_mmsW <- beta.mom(meanS_mmsW, sdS_mmsW)[2]
  
  ##########################################################################
  ##########################################################################
  
  ### msm (adult to spawner) ###
  # *commented out because calculated within simulations
  
  meanS_msm <- input$meanS_msm    # mean late marine survival for all fish
  sdS_msm <- input$sdS_msm        # sd for late marine survival for all fish
  
  #alpha_msm        # alpha parameter calculated from meanS_msm and sdS_msm of late marine survival using beta.mom function
  #beta_msm         # beta parameter calculated from meanS_msm and sdS_msm of late marine survival using beta.mom function
  
  
  ##########################################################################
  ############################## Data Storage ##############################
  ##########################################################################  

  ## make array to hold simulated data
  data <- array(NA,c(gens,10,reps))
  colnames(data) <- c("T_fry", "W_fry", "W_parr", "W_smolt", "W_adult", "R_fry", "R_parr", "R_smolt", "R_adult", "T_spawner")
  
  ## make matrix to hold parameters and results
  parms <- matrix(nrow = reps, ncol = 18)
  colnames(parms) <- c("Kegg", "Ksum", "rescue_level", "realized_rescue_level", "MMS_penalty", "mmsW_survival", "mmsR_survival", "gmean_T_spawner", "sd_T_spawner", "min_T_spawner", "final_T_spawner", "percent5_spawner", "percent50_spawner", "percent95_spawner", "percent_extinction", "percent_rescued", "gmean_T_fry", "gmean_T_smolt") 
  
  ## varying parameters
  parms[,1] <- Kegg                               # record egg K values (Kegg)
  
  cand.K.fry <- c(Ksum, Ksum*0.8, Ksum*0.6, Ksum*0.4, Ksum*0.2) # possible summer carrying capacity values (equal spacing around specified Ksum value)
  parms[,2] <- safeSample(cand.K.fry, reps)       # randomly selected Ksum values (Ksum)
  
  cand.rescue.level <- c(0, rescue_level*0.7, rescue_level, rescue_level*1.3) # possible rescue level values
  parms[,3] <- safeSample(cand.rescue.level, reps)# randomly selected rescue level values (rescue_level)
  
  parms[,5] <- mmsPenalty                         # record penalty value
  
  parms[,6] <- meanS_mmsW                         # record wild early survival values
  
  parms[,7] <- mmsPenalty*parms[,"mmsW_survival"] # calculate values of rescued mms (MMS penalty * mmsW); (mmsR)  
  
  ##########################################################################
  ################################# Model ##################################
  ##########################################################################
  input$submit
  isolate({
  for (sim in 1:reps){
    # variable parameter values (vary between sims)
    Ksum <- parms[sim,"Ksum"]                        # summer capacity values
    rescue_level <- parms[sim,"rescue_level"]        # rescue level values
    meanS_mmsW <- parms[sim,"mmsW_survival"]         # wild early marine survival values
    meanS_mmsR <- parms[sim,"mmsR_survival"]         # rescued early marine survival values
    
    alpha_mmsW <- beta.mom(meanS_mmsW, sdS_mmsW)[1]  # alpha parm calculated from mean,sd for mmsW
    beta_mmsW <- beta.mom(meanS_mmsW, sdS_mmsW)[2]   # beta parm calculated from mean,sd for mmsW
    alpha_mmsR <- beta.mom(meanS_mmsR, sdS_mmsR)[1]  # alpha parm calculated from mean,sd for mmsR
    beta_mmsR <- beta.mom(meanS_mmsR, sdS_mmsR)[2]   # beta parm calculated from mean,sd for mmsR
    alpha_msm <- beta.mom(meanS_msm, sdS_msm)[1]     # alpha parm calculated from mean,sd for msm
    beta_msm <- beta.mom(meanS_msm, sdS_msm)[2]      # beta parm calculated from mean,sd for msm
    
    for (gen in 1:gens) {
      
      ## calculate number of total fry based on previous spawner abundance
      if(gen == 1){
        data[gen,1,sim] <- Rep(initial, fecund, sex_ratio, interceptRep, Kegg, determ)  # use initial abundance for first generation
      } else {
        data[gen,1,sim] <- Rep(data[gen-1,10,sim], fecund, sex_ratio, interceptRep, Kegg, determ) # must come from previous row spawners
      }
########################################################################################     
      
      ## rescued fry and rescued parr *short-term only*
      if (Rlength == 0) {res_FP <- CapST(data[gen,1,sim], rescue_level, facility_cap, encounter_prob, interceptCap, determ)   # output: Rfry and Rparr
      data[gen,6,sim] <- res_FP[1]      # rescued fry
      data[gen,7,sim] <- res_FP[2]      # rescued parr
      
      ## rescued fry and rescued smolt *long-term only*  
      } else {res_FP <- Cap(data[gen,1,sim], rescue_level, facility_cap, encounter_prob, interceptCap, determ)
      data[gen,6,sim] <- res_FP[1]  # rescued fry
      data[gen,8,sim] <- res_FP[2]  # rescued smolt 
      }
      
      ## wild fry and wild parr from total fry and rescue level
      wild_FP <- Sum(data[gen,1,sim], rescue_level, facility_cap, encounter_prob, interceptSum, Ksum, determ)
      data[gen,2,sim] <- wild_FP[1]  # wild fry
      data[gen,3,sim] <- wild_FP[2]  # wild parr
      
      ## all smolt from all parr 
      if (Rlength == 0) {all_PS <- WinST(data[gen,3,sim],data[gen,7,sim],interceptWin, Kwin, determ) # short term
      data[gen,4,sim] <- all_PS[1]  # wild smolt
      data[gen,8,sim] <- all_PS[2]  # rescued smolt
      
      } else {all_PS <- Win(data[gen,3,sim], interceptWin, Kwin, determ) # long term
      data[gen,4,sim] <- all_PS[1]  # wild smolt
      }
      
      ## wild adult from wild smolt
      data[gen,5,sim] <- mmsW(data[gen,4,sim], alpha_mmsW, beta_mmsW, determ)
      
      ## rescued adult from rescued smolt
      data[gen,9,sim] <- mmsR(data[gen,8,sim], alpha_mmsR, beta_mmsR, determ)
      
      ## total spawners from wild + rescued adult
      data[gen,10,sim] <- msm(data[gen,5,sim], data[gen,9,sim], alpha_msm, beta_msm, determ)
      
########################################################################################   
    }
    
    ## outputs
    parms[sim,4] <- mean(data[,6,sim])             # mean realized rescue level
    parms[sim,8] <- exp(mean(log(data[,10,sim])))  # geometric mean total spawner for each simulation
    parms[sim,9] <- sd(data[,10,sim])              # standard deviation of total spawners for each simulation
    parms[sim,10] <- min(data[,10,sim])            # minimum number of total spawners for each simulation
    parms[sim,11] <- data[gens,10,sim]             # final total spawners for each simulation
    parms[sim,12] <- quantile(data[,10,sim], 0.05) # 5th percentile of total spawners
    parms[sim,13] <- quantile(data[,10,sim], 0.50) # 50th percentile of total spanwers
    parms[sim,14] <- quantile(data[,10,sim], 0.95) # 95th percentile of total spawners
    parms[sim,15] <- sum(data[,10,sim]<extn_thresh)/gen # calculate % of generations when spawner abundance falls below extinction threshold
    parms[sim,16] <- mean(data[,6,sim]/data[,1,sim], na.rm = T) # average percent of fry rescued
    parms[sim,17] <- exp(mean(log(data[,1,sim])))               # geometric mean of total fry for each sim
    parms[sim,18] <- exp(mean(log(data[,4,sim]+data[,8,sim])))  # geometric mean of total smolt for each sim
    
  }
  
    ##########################################################################
    ########################### Additional Outputs ###########################
    ##########################################################################
    
    parms <- data.frame(parms)
    
  
  comboData <- list(data = data, parms = parms)
  return(comboData)
  })
  })

  
  usertimeData <- reactive({
    
    ## make matrix to hold time series data
    userData <- modelData()$data[,,modelData()$parms[,"Ksum"] == input$Ksum & modelData()$parms[,"MMS_penalty"] == input$mmsPenalty & modelData()$parms[,"rescue_level"] == input$rescue_level]  # select only matrices with user-specified parameters
    timeseries <- matrix(nrow = input$gens, ncol = 8) # third dimension of userData array is length of timeseries matrix
    colnames(timeseries) <- c("Generation", "percent2.5_spawner", "median_spawner", "percent97.5_spawner", "percent2.5_smolt", "median_smolt", "percent97.5_smolt", "extinction")
    
    ## calculate time series data
    timeseries[,1] <- 1:input$gens   # generation (time)
    
    timeseries[,2] <- apply(userData[,10,], 1, quantile, probs = 0.025)  # lower quartile spawner abundance (25%)
    timeseries[,3] <- apply(userData[,10,], 1, quantile, probs = 0.5)   # median spawner abundance
    timeseries[,4] <- apply(userData[,10,], 1, quantile, probs = 0.975)  # upper quartile spawner abundance (75%)
    
    timeseries[,5] <- apply(userData[,4,] + userData[,8,], 1, quantile, probs = 0.025)  # lower quartile smolt abundance (25%)
    timeseries[,6] <- apply(userData[,4,] + userData[,8,], 1, quantile, probs = 0.5)   # median total smolt abundance
    timeseries[,7] <- apply(userData[,4,] + userData[,8,], 1, quantile, probs = 0.975)  # upper quartile smolt abundance (75%)
    
    timeseries[,8] <- apply(userData[,10,], 1, extinctRisk, threshold = input$extn_thresh, gens = dim(userData)[3])  # extinction risk

    return(data.frame(timeseries))
  })
  
  output$timeTable <- renderTable({
    input$submit
    isolate({
      head(usertimeData())
    })
  })
  
  ## Spawner Timeseries Plot
  output$spawner_timeseries_plot <- renderPlot({
    input$submit
    isolate({
      ggplot(usertimeData(), aes(x = Generation)) +
        geom_ribbon(aes(ymin = percent2.5_spawner, ymax = percent97.5_spawner), alpha = 0.2, fill = "blue") +
        geom_line(aes(y = median_spawner)) +
        #geom_line(aes(y = percent5_spawner), linetype = 2) +
        #geom_line(aes(y = percent95_spawner), linetype = 2) +
        geom_hline(yintercept = input$extn_thresh, color = "red") +
        labs(x = "Generation (3 years/gen)", y = "Spawner Abundance") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    })
  })
  
  ## Smolt Timeseries Plot
  output$smolt_timeseries_plot <- renderPlot({
    input$submit
    isolate({
      ggplot(usertimeData(), aes(x = Generation)) +
        geom_ribbon(aes(ymin = percent2.5_smolt, ymax = percent97.5_smolt), alpha = 0.2, fill = "blue") +
        geom_line(aes(y = median_smolt)) +
        #geom_line(aes(y = percent5_smolt), linetype = 2) +
        #geom_line(aes(y = percent95_smolt), linetype = 2) +
        labs(x = "Generation (3 years/gen)", y = "Smolt Abundance") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    })
  })
  
  ## Extinction Timeseries Plot
  output$extinction_timeseries_plot <- renderPlot({
    input$submit
    isolate({
      ggplot(usertimeData(), aes(x = Generation)) +
        geom_line(aes(y = extinction, color = "red")) +
        labs(x = "Generation (3 years/gen)", y = "Percent Extinction") +
        ylim(0,1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    })
  })

  
  output$table <- renderTable({
    input$submit
    isolate({
    modelData()$parms[1:15,]
    })
  })
  
  #output$table2 <- renderTable({
    #input$submit
    #isolate({
    #data.frame(modelData()$data[,,1])
    #})
  #})
  

  Parameter = c("Number of Simulations", "Initial Abundance", "Extinction Threshold", "Fecundity", "Sex Ratio", "Egg to Fry Survival", "Oversummer Survival", "Overwinter Survival", "Early Marine Survival", "Early Marine Variation", "Late Marine Survival", "Late Marine Variation", "Spawning Capacity (# eggs)", "Number of Fry Rescued", "Captivity Survival", "Rescue Penalty")
  
  output$selected_parameters <- renderTable({
    input$submit
    isolate({
      data.frame(Parameter, Value = c(input$reps, input$initial, input$extn_thresh, input$fecund, input$sex_ratio, input$interceptRep, input$interceptSum, input$interceptWin, input$meanS_mmsW, input$sdS_mmsW, input$meanS_msm, input$sdS_msm, input$Kegg, input$rescue_level, input$interceptCap, input$mmsPenalty))
    })
    })
    
    output$downloadResults <- downloadHandler(
      input$submit,
      filename = function() {
        paste("Rescue_results_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(modelData()$parms, file)
      }
    )
  

  
  ########################  Geometric mean spawners #########################
  ################## against drought (color = rescue level) #################
  output$point_abund_by_drought <- renderPlot({
  input$submit
  isolate({
    ggplot(modelData()$parms, aes(x = Ksum, y = gmean_T_spawner, color = as.factor(rescue_level)) ) +
      geom_jitter(alpha = 0.3) +
      geom_smooth(method = "loess") +
      labs(x = "Maximum Summer Rearing Capacity", y = "Geo Mean Spawner Abundance", color = "Rescue Level") +
      scale_x_continuous(breaks = seq(min(modelData()$parms$Ksum), max(modelData()$parms$Ksum), by = max(modelData()$parms$Ksum*0.2))) +  ## this is the line causing an error "by" must be length 1
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  })
  })
  
  
  ############################  Extinction Risk #############################
  ################## against drought (color = rescue level) #################
  
  output$point_extinct_by_drought <- renderPlot({
    input$submit
    isolate({
      ggplot(modelData()$parms, aes(x = as.factor(Ksum), y = percent_extinction, color = as.factor(rescue_level)) ) +
        geom_boxplot(fill="transparent") +
        labs(x = "Maximum Summer Rearing Capacity", y = "Extinction Risk", color = "Rescue Level") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    })
  })
  
  
}


shinyApp(ui = ui, server = server)

