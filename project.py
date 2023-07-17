import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque

# Make dataframes
pitch_data = pd.read_csv("pitches.csv")
atbat_data = pd.read_csv("atbats.csv")
player_data = pd.read_csv("player_names.csv")
savant_data = pd.read_csv("savant_stats.csv")

# Trim pitch data

# want this:
"""
trimmed_pitch_data <- pitch_data %>% 
    select(ab_id, b_count, s_count, outs, 
    pitch_num, pitch_type, event_num, 
    code, type)
"""

pitch_data = pitch_data[
    [
        "ab_id",
        "b_count",
        "s_count",
        "outs",
        "pitch_num",
        "pitch_type",
        "event_num",
        "code",
        "type",
    ]
]

# Trim atbat data
atbat_data = atbat_data[["ab_id", "batter_id", "event", "o"]]


# Join atbat & pitch, filter
# want this
"""
league_pitch_data <- left_join(trimmed_ab_data, trimmed_pitch_data, by = "ab_id") %>% 
  filter(!((code == "D" & 
              (event == "Grounded Into DP" | 
                 event == "Groundout" | 
                 event == "Forceout" | 
                 event == "Sac Fly")) | 
             code == "Z"))
"""

league_pitch_data = pd.merge(atbat_data, pitch_data, on="ab_id", how="left")
league_pitch_data = league_pitch_data[
    ~(
        (league_pitch_data["code"] == "D")
        & (
            (league_pitch_data["event"] == "Grounded Into DP")
            | (league_pitch_data["event"] == "Groundout")
            | (league_pitch_data["event"] == "Forceout")
            | (league_pitch_data["event"] == "Sac Fly")
        )
        | (league_pitch_data["code"] == "Z")
    )
]

# Trim savant data
# out_zone_swing,iz_contact_percent,in_zone_swing_miss,in_zone_swing,whiff_percent,
savant_data = savant_data[
    [
        "last_name",
        "first_name",
        "player_id",
        "b_swinging_strike",
        "b_total_swinging_strike",
        "oz_swing_percent",
        "out_zone_swing_miss",
        "out_zone_swing",
        "iz_contact_percent",
        "in_zone_swing_miss",
        "in_zone_swing",
        "whiff_percent",
    ]
]


# Calculate swing and miss percent
# want this:
"""
"""
savant_data["scaled_whiff_pct"] = (
    savant_data["whiff_percent"] * savant_data["b_total_pa"]
)
savant_data["scaled_woba"] = savant_data["woba"] * savant_data["b_total_pa"]
savant_data["scaled_chase_rate"] = (
    savant_data["oz_swing_percent"] * savant_data["b_total_pa"]
)


savant_data_agg = savant_data.groupby(["player_id", "last_name", "first_name"]).agg(
    b_swinging_strike=("b_swinging_strike", "sum"),
    b_total_swinging_strike=("b_total_swinging_strike", "sum"),
    oz_swing_percent=("oz_swing_percent", "mean"),
    out_zone_swing_miss=("out_zone_swing_miss", "sum"),
    out_zone_swing=("out_zone_swing", "sum"),
    iz_contact_percent=("iz_contact_percent", "mean"),
    in_zone_swing_miss=("in_zone_swing_miss", "sum"),
    in_zone_swing=("in_zone_swing", "sum"),
    whiff_percent=("whiff_percent", "mean"),
    seasons=("player_id", "count"),
    b_total_pa=("b_total_pa", "sum"),
    b_total_ab=("b_ab", "sum"),
    scaled_whiff_pct=("scaled_whiff_pct", "sum"),
    scaled_woba=("scaled_woba", "sum"),
)


savant_data_agg["swingmiss_rate"] = (
    savant_data_agg["b_swinging_strike"] / savant_data_agg["b_total_swinging_strike"]
)
savant_data_agg["whiff_rate"] = (
    savant_data_agg["scaled_whiff_pct"] / savant_data_agg["b_total_pa"]
)
savant_data_agg["woba"] = savant_data_agg["scaled_woba"] / savant_data_agg["b_total_pa"]


# Calculate woba over all seasons


# Get highest swing miss players
"""
top_swingmiss_player_ids <- swingmiss_data %>% 
  filter(swingmiss_rate > quantile(swingmiss_data$swingmiss_rate, .95)) %>%
  select(player_id)
"""

top_swingmiss_player_ids = (
    savant_data_agg[
        savant_data_agg["swingmiss_rate"]
        > savant_data_agg["swingmiss_rate"].quantile(0.95)
    ]
    .index.get_level_values("player_id")
    .tolist()
)
top_woba_player_ids = (
    savant_data_agg[savant_data_agg["woba"] > savant_data_agg["woba"].quantile(0.95)]
    .index.get_level_values("player_id")
    .tolist()
)


FIELDOUT_EVENTS = [
    "Double Play",
    "Grounded Into DP",
    "Fielders Choice Out",
    "Groundout",
    "Forceout",
    "Flyout",
    "Bunt Groundout",
    "Bunt Lineout",
    "Bunt Pop Out",
    "Pop Out",
    "Lineout",
    "Triple Play",
]
FIELDOUT_CODES = ["X", "E"]
STRIKE_CODES = ["S", "C", "T", "Q", "M", "F", "R", "L"]
BALL_CODES = ["B", "*B", "P", "V"]
BAD_EVENTS = [
    "Runner Out",
    "Catcher Interference",
    "Intent Walk",
    "Sac Bunt",
    "Sacrifice Bunt DP",
    "Batter Interference",
    "Fielders Choice",
]
BAD_CODES = [""]
NONOUT_CODES = ["D", "E"]


states_dict = {
    "K": 12,
    "BB": 13,
    "FO": 14,
    "1B": 15,
    "2B": 16,
    "3B": 17,
    "HR": 18,
    "SF": 19,
    "HBP": 20,
}


def get_player_id(player_name):
    first_name, last_name = player_name.split(" ")
    player_id = player_data[player_data["first_name"] == first_name][
        player_data["last_name"] == last_name
    ]["player_id"].tolist()
    return player_id[0]


def get_player_name(player_id):
    player_name = player_data[player_data["player_id"] == player_id][
        ["first_name", "last_name"]
    ]
    return (
        player_name["first_name"].tolist()[0]
        + " "
        + player_name["last_name"].tolist()[0]
    )


def get_player_data(player_id):
    player_data = pitch_data[
        savant_data_agg.index.get_level_values("batter_id") == player_id
    ]
    return player_data


def generate_markov_chain(player_id):
    data = get_player_data(player_id)
    data = data[~(data["event"].isin(BAD_EVENTS) | data["code"].isin(BAD_CODES))]
    output = []
    for ball_count in range(3):
        for strike_count in range(2):
            # state number, observations in state, transition probabilities
            cur_state = strike_count + 3 * ball_count
            n_observations = len(
                data[
                    (data["b_count"] == ball_count) & (data["s_count"] == strike_count)
                ]
            )
            cur_data = data[
                (data["b_count"] == ball_count) & (data["s_count"] == strike_count)
            ]

            n_strike = len(cur_data[(cur_data["code"].isin(STRIKE_CODES))])
            n_ball = len(cur_data[(cur_data["code"].isin(BALL_CODES))])

            n_fieldout = len(
                cur_data[
                    (cur_data["code"].isin(FIELDOUT_CODES))
                    & (cur_data["event"].isin(FIELDOUT_EVENTS))
                ]
            )
            n_error = len(
                cur_data[
                    (cur_data["code"].isin(NONOUT_CODES)) & cur_data["event"]
                    == "Field Error"
                ]
            )
            n_out = n_fieldout + n_error

            n_sf = len(
                cur_data[
                    (cur_data["code"] == "E")
                    & (
                        cur_data["event"] == "Sac Fly"
                        or cur_data["event"] == "Sac Fly DP"
                    )
                ]
            )

            n_hbp = len(cur_data[(cur_data["code"] == "H")])

            n_single = len(
                cur_data[
                    (cur_data["code"].isin(NONOUT_CODES))
                    & (cur_data["event"] == "Single")
                ]
            )
            n_double = len(
                cur_data[
                    (cur_data["code"].isin(NONOUT_CODES))
                    & (cur_data["event"] == "Double")
                ]
            )
            n_triple = len(
                cur_data[
                    (cur_data["code"].isin(NONOUT_CODES))
                    & (cur_data["event"] == "Triple")
                ]
            )
            n_hr = len(
                cur_data[
                    (cur_data["code"].isin(NONOUT_CODES))
                    & (cur_data["event"] == "Home Run")
                ]
            )

            obs_count = [0] * 20

            if ball_count == 3:
                obs_count[states_dict["BB"]] = n_ball / n_observations
            else:
                obs_count[cur_state + 3] = n_ball / n_observations

            if strike_count == 2:
                nonfoul_strikes = ["W", "S", "C", "T", "Q", "M"]
                n_strike = len(cur_data[(cur_data["code"].isin(nonfoul_strikes))])
                obs_count[states_dict["K"]] = n_strike / n_observations
            else:
                obs_count[cur_state + 1] = n_strike / n_observations

            obs_count[states_dict["1B"]] = n_single / n_observations
            obs_count[states_dict["2B"]] = n_double / n_observations
            obs_count[states_dict["3B"]] = n_triple / n_observations
            obs_count[states_dict["HR"]] = n_hr / n_observations
            obs_count[states_dict["SF"]] = n_sf / n_observations
            obs_count[states_dict["HBP"]] = n_hbp / n_observations
            total_nonout = sum(obs_count)
            obs_count[states_dict["FO"]] = 1 - total_nonout
            output.append([obs_count])
    for i in range(9):
        row_to_add = [0] * 21
        row_to_add[i + 12] = 1
        output.append([row_to_add])

    # want to return output as a stochastic matrix
    return np.array(output).reshape(21, 21)


def calculate_absorption_time(transition_matrix):
    # Get the 12x12 matrix of transient states
    Q = transition_matrix[:12, :12]
    F1 = np.linalg.inv(np.identity(12) - Q)
    F2 = np.ones((12, 1))
    F = np.matmul(F1, F2)
    counts = [
        "0-0",
        "0-1",
        "0-2",
        "1-0",
        "1-1",
        "1-2",
        "2-0",
        "2-1",
        "2-2",
        "3-0",
        "3-1",
        "3-2",
    ]

    return pd.DataFrame({"count": counts, "time_to_abs": F.flatten()})


#TODO: Don't need to dynamically adjust order of graph, can be fixed from the start just need to do Topo sort
def calculate_count_woba(transition_matrix):
    visited = [0] * 12
    visited[11] = 1
    e_woba = [0] * 12
    queue = deque([(3, 2)])
    while queue:
        cur_state = queue.popleft()
        ball, strike = cur_state
        state_num = ball * 3 + strike
        visited[state_num] = 1
        #        E_absorption <- 0.89 * transition_mat[state_number, 16] +
        #   1.27 * transition_mat[state_number, 17] +
        #   1.62 * transition_mat[state_number, 18] +
        #   2.10 * transition_mat[state_number, 19] +
        #   0.72 * transition_mat[state_number, 21]
        #   want this ^
        e_woba_count = (
            0.89 * transition_matrix[state_num][16]
            + 1.27 * transition_matrix[state_num][17]
            + 1.62 * transition_matrix[state_num][18]
            + 2.10 * transition_matrix[state_num][19]
            + 0.72 * transition_matrix[state_num][21]
        )

        if ball == 3 and strike == 2:
            e_walk = 0.69 * transition_matrix[state_num][14]
            # Formula for infinite series sum a/1-r to account for foul ball probability.
            e_woba_count = (e_woba_count + e_walk) / (1 - transition_matrix[state_num][state_num])
            e_woba[state_num] = e_woba_count
        elif ball == 3:
            e_walk = 0.69 * transition_matrix[state_num][14]
            e_strike = transition_matrix[state_num][state_num + 1] * e_woba[state_num + 1]
            e_woba[state_num] = e_woba_count + e_walk  + e_strike
        elif strike == 2:
            e_ball = transition_matrix[state_num][state_num + 3] * e_woba[state_num + 3]
            e_woba[state_num] = (e_woba_count + e_ball) / (1 - transition_matrix[state_num][state_num])
        else:
            e_ball = transition_matrix[state_num][state_num + 3] * e_woba[state_num + 3]
            e_strike = transition_matrix[state_num][state_num + 1] * e_woba[state_num + 1]
            e_woba[state_num] = e_woba_count + e_ball + e_strike
        
        if ball > 0 and not visited[3*ball  + strike - 2]:
            queue.append((ball, strike - 1))
            
        

        for i in range(12):
            if (
                transition_matrix[cur_state[0] * 3 + cur_state[1]][i] > 0
                and visited[i] == 0
            ):
                queue.append((i // 3, i % 3))


def calculate_woba_and_abs_time(data, name):
    transition_matrix = generate_markov_chain(data)
    absorption_time = calculate_absorption_time(transition_matrix)
    woba = calculate_count_woba(transition_matrix)
    values = absorption_time.merge(woba, on="count")
    values["hitter"] = name



def get_player_values(fname, lname):
    player_name = fname + " " + lname
    player_id = get_player_id(player_name)
    player_data = get_player_data(player_id)
    return calculate_woba_and_abs_time(player_data, player_name)



trout_values = get_player_values("Mike", "Trout")

league_transition_mat = generate_markov_chain(league_pitch_data)
league_values = calculate_woba_and_abs_time(league_transition_mat, "League Average")


values = pd.concat([trout_values, league_values])

#now want this:
# ggplot(data = all_vals, aes(x=count, y=woba, fill=hitter)) +
#  geom_bar(stat='identity', position = 'dodge') + 
#  ggtitle("Expected wOBA From Each Count: Mike Trout vs. League Average")

plt.bar(values["count"], values["woba"], values["hitter"])
plt.show()