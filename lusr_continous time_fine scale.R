library(deSolve)

# parameters
OD <- 18500000 # cell number in 1mL when OD = 1
volume <- 0.128 # uL, medium volume
Ninitial <- volume*OD # initial population size of both population
freq <- c(0.1,0.9) # initial frequency of leu- strain

dilution <- 2^(-10) # dilution factor
transfer <- 3 # transfer times
t1 <- 300 #  the first sererial dilution
t2 <- t1*2 #  the second sererial dilution
timef <- t1*3 # total time period

# leu- strain ecological parameters
Rm1 <- c(0.1459, 0.1912, 0.1989, 0.1708) # all parameters used for Rmax1
K1 <- c(51.2176, 63.4227, 33.7455, 29.6791) # for Km1
Y1 <- c(0.0329, 0.0327, 0.0289, 0.0311) # yield of leu- strain for leucine
Y1 <- 10000000000*Y1

# trp- strain ecological parameters
Rm2 <- c(0.17883, 0.179907, 0.170852, 0.173637, 0.173451)  # for Rmax2
K2 <- c(3.3646, 3.4371, 2.9538, 2.8635, 2.3002)  # for K2
Y2 <- c(0.2119, 0.2073, 0.2073, 0.2134, 0.1907) # yield of trp- strain for tryptophan
Y2 <- 10000000000*Y2

# file scale nutrient levels
Trp <- c(0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15) #Specific nutrient for population 1 - different concentration of trp (mg/L), fold 0.02 ~ 0.5
Leu <- c(0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15) #% % Specific nutrient for population 2 - different concentration of leu (mg/L),fold 0.02 ~ 0.25

# change unit to for X to mg/L
Trp <- 73.4*Trp  #concetration of nutients (mg/L)
Leu <- 367*Leu # concentration of leucine (mg/L) 

Clusr <- 1 # coefficient for LUSR

# store results - final frequency of leu strain
ffreqLeu <- array(0, dim = c(length(Leu), length(Trp), length(freq)))

for (i in 1:length(Leu)) {
  
  for (j in 1:length(Trp)) {
    
    for (m in 1:length(freq)) {
      
      # ODE functions
      lusr <- function(t, y, parms){
        with(as.list(c(y, parms)), {
          
          # reset lusr coefficient if nutrient goes to 0.1
          if (X >= 0.1){
            Clusr2 <- Clusr
          } else {
            Clusr2 <- 0
          }
          if (Y >= 0.1){
            Clusr1 <- Clusr
          } else {
            Clusr1 <- 0
          }
          
          r1 = Rmax1 * X / (Km1 + X)
          r2 = Rmax2 * Y / (Km2 +Y)
          # double check that growth rate can't be negative 
          if (r1 < 0) {r1 = 0}
          if (r2 < 0) {r2 = 0}
          # print(c(t, r1, r2))
          
          # population growth
          dL = L * r1   # leu strain
          dP = P * r2   # trp strain
          
          # amino acid depletion
          dX =  - dL/Y1leu - Clusr2 * dP/Y2leu # leucine
          dY =  - dP/Y2trp - Clusr1 * dL/Y1trp # tryptophan
          
          return(list(c(dL, dP, dX, dY)))
        })
      }
      
      # initial values
      yini <- c(L = Ninitial * freq[m], P = Ninitial * (1 - freq[m]), X = Leu[i], Y = Trp[j])
      
      # solve ODE
      times <- seq(from = 0, to = timef, by = 1)
      
      parms <- c(Rmax1 = Rm1[1], Rmax2 = Rm2[1],
                 Km1 = K1[1], Km2 = K2[1],
                 Y1leu = Y1[1], Y2leu = Y1[1],
                 Y1trp = Y2[1], Y2trp = Y2[1])
      
      Serialdilution <- data.frame( 
        var =  c("L", "P", "X", "Y", "L", "P", "X", "Y"), # name of the variable, 
        time = c(t1, t1, t1, t1, t2, t2, t2, t2), # at which times do the events occur,
        value = c(dilution, dilution, Leu[i], Trp[j], dilution, dilution, Leu[i], Trp[j]), 
        method = c("mult", "mult", "rep", "rep", "mult", "mult", "rep", "rep"))
      
      out2 <- ode(times = times, y = yini, 
                  func = lusr, parms = parms, 
                  atol = 0,
                  method="ode45",
                  events = list(data = Serialdilution))
      out2_df <- as.data.frame(out2)
      
      # plot 
      # plot(out2, select = "X") #leucine
      # plot(out2, select = "Y") #tryptophan
      # plot(out2, select = "L") #leu strain
      # plot(out2, select = "P")
      
      # store final frequency of Leu strain
      ffreqLeu[(length(Leu)+1)-i,j,m] <- out2[nrow(out2),2]/(out2[nrow(out2),2] + out2[nrow(out2),3])
    }
  }
}

# organize final frequency of leu strain
f1 <- as.data.frame(ffreqLeu[ , , 1]) # final frequency of leu strain when start with 0.1 initial frequency
f9 <- as.data.frame(ffreqLeu[ , , 2]) # final frequency of leu strain when start with 0.9 initial frequency

# plot final frequency
library("plotly") 
plot_ly(z = as.matrix(f1[rev(rownames(f1)), ]), 
        type = "heatmap",
        zmin=0, zmax=1) 
plot_ly(z = as.matrix(f9[rev(rownames(f9)), ]), type = "heatmap",
        zmin=0, zmax=1) 

# final dynamic for coexistence
finaldynamic <- matrix(0, nrow = dim(f1)[1], ncol = dim(f1)[2])

# encode final dynamic
# 1: coexistence
# 2: trp- wins
# 3: leu- wins
for (i in 1:dim(f1)[1]) {
  for (j in 1:dim(f1)[2]) {
    
    if(f1[i,j] > 0.1 && f9[i, j] < 0.9) { # coexistence
      finaldynamic[i, j] <- 1
    } else if(f1[i,j] < 0.1 && f9[i, j] < 0.9) { # trp- wins
      finaldynamic[i, j] <- 2
    } else if(f1[i,j] > 0.1 && f9[i, j] > 0.9) { # leu- wins
      finaldynamic[i, j] <- 3
    } else {
      finaldynamic[i, j] <- 4
    }
  }
}

# # store results in excel
# library(openxlsx) # Load on each new session
# library(xlsx)
# write.csv(finaldynamic, "finefinal_lusr0.csv")
# write.csv(finaldynamic, "finefinal_lusr05.csv")
# write.csv(finaldynamic, "finefinal_lusr1.csv")
