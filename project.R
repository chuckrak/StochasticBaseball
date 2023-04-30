setwd("/Users/chuckrak/math156/finalProject")
library(tidyverse)
pitch_data <- read.csv('pitches.csv') 
player_names <- read.csv('player_names.csv')
atbat_data <- read.csv('atbats.csv')
whiff_data <- read.csv('stats.csv')

trimmed_pitch_data <- pitch_data %>% select(ab_id, b_count, s_count, outs, 
                                            pitch_num, pitch_type, event_num, 
                                            code, type)
trimmed_ab_data <- atbat_data %>% select(ab_id, batter_id, event, o)

#filter bad codes / errors in data
league_pitch_data <- left_join(trimmed_ab_data, trimmed_pitch_data, by = "ab_id") %>% 
  filter(!((code == "D" & 
              (event == "Grounded Into DP" | 
                 event == "Groundout" | 
                 event == "Forceout" | 
                 event == "Sac Fly")) | 
             code == "Z"))

league_pitch_data %>% group_by(event,code) %>% summarize(n()) %>% View()


swingmiss_data <- whiff_data %>% 
  select(last_name, first_name, player_id, b_swinging_strike, b_total_swinging_strike)
  
swingmiss_data <- whiff_data %>% group_by(player_id, last_name, first_name) %>% 
  summarize(total_swing_miss = sum(b_swinging_strike), 
            total_swing = sum(b_total_swinging_strike), 
            seasons = n()) %>%
  mutate(swingmiss_rate = total_swing_miss / total_swing)

woba_data <- whiff_data %>% 
  mutate(wobaxpa = woba * b_total_pa) %>%
  group_by(player_id, last_name, first_name) %>% 
  summarize(total_pa = sum(b_total_pa), 
            seasons = n(),
            wobaxpa = sum(wobaxpa)) %>%
  mutate(woba = wobaxpa / total_pa)





top_swingmiss_player_ids <- swingmiss_data %>% 
  filter(swingmiss_rate > quantile(swingmiss_data$swingmiss_rate, .95)) %>%
  select(player_id)

top_woba_player_ids <- woba_data %>% 
  filter(woba > quantile(woba_data$woba, .95)) %>%
  select(player_id)




get_player_id <- function(fname, lname){
  player <- player_names %>% filter(first_name == fname & last_name == lname)
  return(player$id)
}

get_player_name <- function(pid){
  player <- player_names %>% filter(id == pid)
  return(c(player$first_name, player$last_name))
}

get_player_data <- function(player_id){
  pbp_data <- league_pitch_data %>% filter(batter_id == player_id)
  return(pbp_data)
}

counts <- c("0-0", "0-1", "0-2", "1-0", "1-1", "1-2", "2-0", "2-1", "2-2", "3-0", "3-1", "3-2")

#state indices
#counts-1-12
#strikeout - 13
#walk - 14
#field out - 15
#single - 16
#double 17
#triple 18
#home run 19
#sacrifice fly 20
#hit by pitch 21

strikeout_state <- 13
walk_state <- 14
fieldout_state <- 15
single_state <- 16
double_state <- 17
triple_state <- 18
homerun_state <- 19
sacrificefly_state <- 20
hitbypitch_state <- 21

generate_markov_chain <- function(player_data){
  data <- player_data %>% 
    filter(!(event == "Runner Out" | 
               event == "Catcher Interference" | 
               event == "Intent Walk" | 
               event == "Sac Bunt" | 
               code == "" | 	
               event == "Sacrifice Bunt DP" | 
               event == "Batter Interference" | 
               event == "Fielders Choice"))  
  output <- c()
  for(ball_count in c(0:3)){
    for(strike_count in c(0:2)){
      # 0-0, 0-1, 0-2, 1-0, 1-1, 1-2, 2-0, 2-1, 2-2, 3-0, 3-1, 3-2
      
      #print("Starting with current ball and strike amounts")
      #print(ball_count)
      #print(strike_count)
      state_number <- strike_count + 3 * ball_count + 1
      
      #total obs in state
      count <- data %>% 
        filter(b_count == ball_count & s_count == strike_count) %>%
        summarize(n()) %>%
        as.integer()
      
      #print("number of observations")
      #print(count)
      
      #state goes up one number
      count_strike <- data %>% 
        filter(b_count == ball_count & 
                 s_count == strike_count &
                 (code == "W" |
                    code == "S" |
                    code == "C" |
                    code == "T" |
                    code == "Q" | 
                    code == "M" | 
                    code == "F" |
                    code == "R" |
                    code == "L")) %>% 
        summarize(n()) %>% 
        as.integer()
      #state goes up three numbers
      count_ball <- data %>% 
        filter(b_count == ball_count &
                 s_count == strike_count &
                 (code == "B" | 
                    code == "*B" | 
                    code == "P" |
                    code == "V")) %>%
        summarize(n()) %>%
        as.integer()
      count_fieldout <- data %>% 
        filter(b_count == ball_count & s_count == strike_count & 
                 (code == "X" | code == "E") & 
                 (event == "Double Play" | 
                    event == "Grounded Into DP" |
                    event == "Fielders Choice Out" |
                    event == "Groundout" | 
                    event == "Forceout" | 
                    event == "Flyout" | 
                    event == "Bunt Groundout" | 
                    event == "Bunt Lineout" | 
                    event == "Bunt Pop Out" | 
                    event == "Pop Out" | 
                    event == "Lineout" | 
                    event == "Triple Play")) %>%
        summarize(n()) %>%
        as.integer()
      count_error <- data %>% 
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 (code == "D" | 
                    code == "E") &
                 event == "Field Error") %>% 
        summarize(n()) %>%
        as.integer()      
      # errors to be treated as outs, as is standard
      count_out <- count_fieldout + count_error
      count_sacrifice <- data %>%
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 code == "E" & 
                 (event == "Sac Fly" | 
                    event == "Sac Fly DP")) %>% 
        summarize(n()) %>% 
        as.integer()
      count_hitbypitch <- data %>%
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 code == "H") %>% 
        summarize(n()) %>%
        as.integer()
      count_single <- data %>%
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 (code == "D" |
                    code == "E" |
                    type == "X") & 
                 event == "Single") %>% 
        summarize(n()) %>%
        as.integer()
      count_double <- data %>% 
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 (code == "D" | 
                    code == "E" |
                    type == "X") & 
                 event == "Double") %>% 
        summarize(n()) %>% 
        as.integer()
      count_triple <- data %>%
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 (code == "D" | 
                    code == "E" | 
                    type == "X") & 
                 event == "Triple") %>% 
        summarize(n()) %>%
        as.integer()
      count_homerun <- data %>%
        filter(b_count == ball_count & 
                 s_count == strike_count & 
                 (code == "D" | 
                    code == "E" |
                    type == "X") & 
                 event == "Home Run") %>% 
        summarize(n()) %>% 
        as.integer()
    
      
      #stores chance of going from on state in a count to another
      count_transition_vector <- c(rep(0, 20))
      #only non-zero can values be count + 1 (another strike) or count + 3 (another ball)
      
      if(ball_count == 3){
        count_transition_vector[walk_state] <- count_ball
      } else{
        count_transition_vector[state_number + 3] <- count_ball
      }
      if(strike_count == 2) {
        count_strike <- data %>% 
          filter(b_count == ball_count & 
                   s_count == strike_count & 
                   (code == "W" | 
                      code == "S" | 
                      code == "C" | 
                      code == "T" | 
                      code == "Q" | 
                      code == "M")) %>% 
          summarize(n()) %>% 
          as.integer()
        count_foul <- data %>% 
          filter(b_count == ball_count & 
                   s_count == strike_count & 
                   (code == "F" | 
                      code == "R")) %>% 
          summarize(n()) %>%
          as.integer()
        count_transition_vector[state_number] <- count_foul
        count_transition_vector[strikeout_state] <- count_strike
      } else{
        count_transition_vector[state_number + 1] <- count_strike
      }
    
      count_transition_vector[single_state] <- count_single
      count_transition_vector[double_state] <- count_double
      count_transition_vector[triple_state] <- count_triple
      count_transition_vector[homerun_state] <- count_homerun
      count_transition_vector[sacrificefly_state] <- count_sacrifice
      count_transition_vector[hitbypitch_state] <- count_hitbypitch
      
      # convert to probabilities by dividing by total number of events
      count_transition_vector <- count_transition_vector / count
      count_transition_vector[fieldout_state] <- 1-sum(count_transition_vector)
      output <- c(output, count_transition_vector)
    }
  }
  
  last_9_rows <- c(rep(0, 189))
  last_9_rows[13] <- 1
  last_9_rows[35] <- 1
  last_9_rows[57] <- 1
  last_9_rows[79] <- 1
  last_9_rows[101] <- 1
  last_9_rows[123] <- 1
  last_9_rows[145] <- 1
  last_9_rows[167] <- 1
  last_9_rows[189] <- 1
  output <- c(output, last_9_rows)
  return(output)
}

calculate_absorption_time <- function(markov_chain){
  transition_mat <- matrix(markov_chain, nrow=21, byrow=TRUE)
  Q <- transition_mat[1:12, 1:12]
  F1 <- solve(diag(12) - Q) %*% matrix(c(rep(1,12)), nrow=12, byrow=TRUE)
  return(data.frame("count" = counts,"time_to_abs" = F1))
}

calculate_count_woba <- function(markov_chain){
  transition_mat <- matrix(markov_chain, nrow=21, byrow=TRUE)
  woba <- c(rep(0, 12))
  visited <- c(rep(FALSE,12))
  
  visited[12] <- TRUE
  queue <- c(3,2)
  # While there are nodes left to visit...
  while(length(queue) > 0) {
    b_count <- queue[1]
    s_count <- queue[2]
    #remove first 2
    queue <- queue[-1][-1]
    state_number <- s_count + 3 * b_count + 1
    #expected value from at bat ending events at this count
    E_absorption <- 0.89 * transition_mat[state_number, 16] + 
      1.27 * transition_mat[state_number, 17] + 
      1.62 * transition_mat[state_number, 18] +
      2.10 * transition_mat[state_number, 19] +
      0.72 * transition_mat[state_number, 21]
    
    if(b_count == 3 && s_count == 2){
      #walk p * walk value
      E_walk <- transition_mat[state_number, 14] * 0.69
      #Formula for infinite series sum a/1-r to account for foul ball probability.
      woba[state_number] <- 
        (E_absorption + E_walk) / (1 - transition_mat[state_number, state_number])
    } else if(b_count == 3){
      E_walk <- transition_mat[state_number, 14] * 0.69
      woba[state_number] <- E_walk + 
        E_absorption +
        transition_mat[state_number, state_number + 1] * woba[state_number + 1]
      
    } else if(s_count == 2){
      E_ball <- transition_mat[state_number, state_number + 3] * woba[state_number + 3]
      #Formula for infinite series sum a/1-r to account for foul ball probability.
      woba[state_number] <- 
        (E_absorption + E_ball) / (1 - transition_mat[state_number, state_number])
    } else {
      woba[state_number] <- E_absorption + 
        transition_mat[state_number,state_number+3] * woba[state_number + 3] + 
        transition_mat[state_number, state_number + 1] * woba[state_number + 1]
    }
    
    
  
    if(b_count > 0 && !visited[3*(b_count - 1) + s_count + 1]){
      visited[3*(b_count - 1) + s_count + 1] <- TRUE
      queue <- c(queue, b_count-1, s_count)
    }
    
    if(s_count > 0 && !visited[3*b_count + s_count]){
      visited[3*b_count + s_count] <- TRUE
      queue <- c(queue, b_count, s_count - 1)
    }
    
    
    
  }
  return (data.frame("count" = counts, woba))
}


get_player_values <- function(fname, lname){
  player_name <- paste(fname, lname)
  player_id <- get_player_id(fname, lname)
  player_data <- get_player_data(player_id)
  player_MC <- generate_markov_chain(player_data)
  player_abs_time <- calculate_absorption_time(player_MC)
  player_count_woba <- calculate_count_woba(player_MC)
  names_col <- c(rep(player_name, 12))
  values <- player_abs_time %>% left_join(player_count_woba, by="count")
  values <- cbind(values, "hitter" = names_col)
  return(values)
}



trout_values <- get_player_values("Mike", "Trout")

league_markov_chain <- generate_markov_chain(league_pitch_data) 
league_abs_time <- calculate_absorption_time(league_markov_chain)
league_count_woba <- calculate_count_woba(league_markov_chain)

league_values <- left_join(league_abs_time, league_count_woba, by="count") %>% 
  cbind("hitter" = c(rep("League Average", 12)))

all_vals <- rbind(league_values, trout_values)

ggplot(data = all_vals, aes(x=count, y=woba, fill=hitter)) +
  geom_bar(stat='identity', position = 'dodge')


aggressive_hitters_df <- data.frame()
for(player_id in top_swingmiss_player_ids$player_id){
  player_name <- get_player_name(player_id)
  first <- player_name[1]
  last <- player_name[2]
  vals <- get_player_values(first,last)
  aggressive_hitters_df <- rbind(aggressive_hitters_df, vals)
}
aggressive_hitters_df <- aggressive_hitters_df %>% 
  group_by(count) %>%
  summarize(time_to_abs = mean(time_to_abs),
            woba = mean(woba)) %>%
  cbind("hitter" = c(rep("Aggressive", 12)))

aggressive_league_df <- rbind(aggressive_hitters_df, league_values)

aggressive_hitters_ahead_counts <- aggressive_league_df %>%
  filter(count=="2-0" | count == "3-0" | count == "3-1")

ggplot(data = aggressive_hitters_ahead_counts, aes(x=count, y=woba, fill=hitter)) +
  geom_bar(stat='identity', position = 'dodge')



woba_leaders_df <- data.frame()
for(player_id in top_woba_player_ids$player_id){
  player_name <- get_player_name(player_id)
  first <- player_name[1]
  last <- player_name[2]
  vals <- get_player_values(first,last)
  woba_leaders_df <- rbind(woba_leaders_df, vals)
}

woba_leaders_df <- woba_leaders_df %>% 
  group_by(count) %>%
  summarize(time_to_abs = mean(time_to_abs),
            woba = mean(woba)) %>% 
  cbind("hitter" = c(rep("Top 5% WOBA", 12)))


woba_league_df <- rbind(woba_leaders_df, league_values)



ggplot(data = woba_league_df, aes(x=count, y=woba, fill=hitter)) +
  geom_bar(stat='identity', position = 'dodge')


#sample 50 hitters, compare 0-2 abs_time to 0-2 woba, see what happens
# do same for other counts where hitters are behind.

player_sample <- woba_data$player_id %>% sample(50, replace=FALSE)
woba_abs_time <- data.frame()

for(player_id in player_sample){
  player_name <- get_player_name(player_id)
  first <- player_name[1]
  last <- player_name[2]
  vals <- get_player_values(first,last)
  woba_abs_time <- rbind(woba_abs_time, vals)
}

woba_abs_time %>% 
  filter(count == "0-2") %>%
  ggplot(aes(x=time_to_abs, y=woba)) + geom_point()
  
woba_abs_time %>% 
  group_by(count) %>%
  summarize(r2 = cor(time_to_abs, woba)^2)




