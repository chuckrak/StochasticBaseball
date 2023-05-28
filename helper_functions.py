import numpy as np
import pandas as pd
import project

FIELDOUT_EVENTS = ['Double Play', 'Grounded Into DP', 'Fielders Choice Out', 'Groundout', 'Forceout', 'Flyout', 'Bunt Groundout', 'Bunt Lineout', 'Bunt Pop Out', 'Pop Out', 'Lineout', 'Triple Play']
FIELDOUT_CODES = ['X', 'E']
STRIKE_CODES = ['S', 'C', 'T', 'Q', 'M', 'F', 'R', 'L']
BALL_CODES = ['B', '*B', 'P', 'V']
BAD_EVENTS = ['Runner Out', 'Catcher Interference', 'Intent Walk', 'Sac Bunt', 'Sacrifice Bunt DP', 'Batter Interference', 'Fielders Choice']
BAD_CODES = ['']
NONOUT_CODES = ['D', 'E']


states_dict = {"K": 12, "BB": 13, "FO": 14, "1B": 15, "2B": 16, "3B": 17, "HR": 18, "SF": 19, "HBP": 20}

def get_player_id(player_name):
    player_name = player_name.split(' ')
    first_name = player_name[0]
    last_name = player_name[1]
    player_id = player_data[player_data['first_name'] == first_name][player_data['last_name'] == last_name]['player_id'].tolist()
    return player_id[0]

def get_player_name(player_id):
    player_name = player_data[player_data['player_id'] == player_id][['first_name', 'last_name']]
    return player_name['first_name'].tolist()[0] + ' ' + player_name['last_name'].tolist()[0]

def get_player_data(player_id):
    player_data = pitch_data[savant_data_agg.index.get_level_values('batter_id') == player_id]
    return player_data




def generate_markov_chain(player_id):
    data = get_player_data(player_id)
    data = data[~(data['event'].isin(BAD_EVENTS) | data['code'].isin(BAD_CODES))]
    output = []
    for ball_count in range(3):
        for strike_count in range(2):
            # state number, observations in state, transition probabilities
            cur_state = strike_count + 3 * ball_count
            n_observations = len(data[(data['b_count'] == ball_count) & (data['s_count'] == strike_count)])
            cur_data = data[(data['b_count'] == ball_count) & (data['s_count'] == strike_count)]


            n_strike = len(cur_data[(cur_data['code'].isin(STRIKE_CODES))])
            n_ball = len(cur_data[(cur_data['code'].isin(BALL_CODES))])
            
            
            n_fieldout = len(cur_data[(cur_data['code'].isin(FIELDOUT_CODES)) & (cur_data['event'].isin(FIELDOUT_EVENTS))])
            n_error = len(cur_data[(cur_data['code'].isin(NONOUT_CODES)) & cur_data['event'] == "Field Error"])
            n_out = n_fieldout + n_error

            n_sf = len(cur_data[(cur_data['code'] == "E") & (cur_data['event'] == 'Sac Fly' or cur_data['event'] == 'Sac Fly DP')])
            
            n_hbp = len(cur_data[(cur_data['code'] == 'H')])

            n_single = len(cur_data[(cur_data['code'].isin(NONOUT_CODES)) & (cur_data['event'] == 'Single')])
            n_double = len(cur_data[(cur_data['code'].isin(NONOUT_CODES)) & (cur_data['event'] == 'Double')])
            n_triple = len(cur_data[(cur_data['code'].isin(NONOUT_CODES)) & (cur_data['event'] == 'Triple')])
            n_hr = len(cur_data[(cur_data['code'].isin(NONOUT_CODES)) & (cur_data['event'] == 'Home Run')])

            obs_count = [0] * 20
            
            if ball_count == 3:
                obs_count[states_dict['BB']] = n_ball / n_observations
            else:
                obs_count[cur_state + 3] = n_ball / n_observations
            
            if strike_count == 2:
                nonfoul_strikes = ["W", "S", "C", "T", "Q", "M"]
                n_strike = len(cur_data[(cur_data['code'].isin(nonfoul_strikes))])
                obs_count[states_dict['K']] = n_strike / n_observations
            else:
                obs_count[cur_state + 1] = n_strike / n_observations

            obs_count[states_dict['1B']] = n_single / n_observations
            obs_count[states_dict['2B']] = n_double / n_observations
            obs_count[states_dict['3B']] = n_triple / n_observations
            obs_count[states_dict['HR']] = n_hr / n_observations
            obs_count[states_dict['SF']] = n_sf / n_observations
            obs_count[states_dict['HBP']] = n_hbp / n_observations
            total_nonout = sum(obs_count)
            obs_count[states_dict['FO']] = 1 - total_nonout
            output.append([obs_count])
    for i in range(9):
        row_to_add = [0] * 21
        row_to_add[i + 12] = 1
        output.append([row_to_add])

    # want to return output as a stochastic matrix
    return np.array(output).reshape(21, 21)

def calculate_absorption_time(transition_matrix):
    #Get the 12x12 matrix of transient states
    Q = transition_matrix[:12, :12]
    F1 = np.linalg.inv(np.identity(12) - Q)
    F2 = np.ones((12, 1))
    F = np.matmul(F1, F2)
    #return data.frame("count"=counts, "time_to_abs"= F)
    #want this ^
    counts = ["0-0", "0-1", "0-2", "1-0", "1-1", "1-2", "2-0", "2-1", "2-2", "3-0", "3-1", "3-2"]

    return 




            




