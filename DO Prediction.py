import numpy as np
from numba import njit
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import mean_absolute_error
import os
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# --- Helper function to explicitly calculate DO values ---
# JIT-compiled function for simulating DO levels for a 24-hour period
@njit
def simulate_day_DO(n_steps, n_sim, tempC, irr, press, salinity, avg_temp, last_DO, constants, sigma):
    predictions = np.full((n_steps, n_sim), np.nan) # Store predicted DO values
    prev_DO_matrix = np.full((n_steps, n_sim), np.nan) # Store the previous DO used for each step

    # Loop over Monte Carlo simulations
    for sim in range(n_sim):
        # Sample once per simulation
        K, theta, R, A, p = constants[sim]

        # Add random noise to the initial DO measurement
        prev_DO = last_DO + np.random.normal(0, sigma)

        # Loop over time steps for the day
        for j in range(n_steps):
            # For the first step, use the first row of measurements (which corresponds to the last measurement of the previous day)
            if j == 0:
                temp = tempC[0]
                irrad = irr[0]
                press_j = press[0]
                sal_j = salinity[0]
                avgT = avg_temp[0]
            # For subsequent steps, use the current step's measurements
            else:
                temp = tempC[j - 1]
                irrad = irr[j - 1]
                press_j = press[j - 1]
                sal_j = salinity[j - 1]
                avgT = avg_temp[j - 1]

            # Convert temperature to Kelvin
            T_k = temp + 273.15

            # Compute oxygen solubility terms
            S1 = 157570.1 / T_k
            S2 = -6.6423080E7 / T_k**2
            S3 = 1.2438E10 / T_k**3
            S4 = -8.621949E11 / T_k**4
            S5 = -sal_j * (0.017674 - 10.754/T_k + 2140.7/T_k**2)

            # Compute DO saturation concentration with salinity corrections
            DO_sat_sal = np.exp(-139.34411 + S1 + S2 + S3 + S4 + S5)

            # Pressure correction terms
            alpha = 0.000975 - 0.00001426*T_k + 0.00000006436*T_k**2
            beta = np.exp(11.8571 - (3840.7 / T_k) - (216961 / T_k**2))
            gamma = ((1 - beta / press_j) / (1 - beta)) * ((1 - alpha * press_j) / (1 - alpha))

            # Final DO saturation value with salinity and pressure corrections
            DO_sat = DO_sat_sal * press_j * gamma

            # Store the previous DO
            prev_DO_matrix[j, sim] = prev_DO

            # Compute Gross Primary Production (GPP)
            GPP = A * (irrad ** p) if irrad > 0 else 0

            # Compute DO using regression formula from Grace et al. 2015
            DO = prev_DO + GPP - R * (theta ** (temp - avgT)) + K * (DO_sat - prev_DO) * (1.0241 ** (temp - avgT))
            
            # Bound DO to a reasonable range
            predictions[j, sim] = min(max(DO, 0), 30)

            # Update prev_DO for the next timestep
            prev_DO = DO

    return predictions, prev_DO_matrix

# ----- Load & Preprocess Data -----
def load_data(measurement_data, metabolism_data):
    # Load measurement data
    measurements = pd.read_csv(measurement_data)
    measurements['datetime'] = pd.to_datetime(measurements['Date'] + ' ' + measurements['Time'])

    # Sort data chronologically
    measurements = measurements.sort_values('datetime').reset_index(drop=True)

    # Load metabolism data with start time
    metabolism = pd.read_csv(metabolism_data)
    metabolism['period_start'] = pd.to_datetime(metabolism['Start Date'] + ' ' + metabolism['Start Time'])

    # Define end time for each metabolism period (24 hours after start)
    metabolism['period_end'] = metabolism['period_start'] + pd.Timedelta(days=1)

    # For clarity, define a "period_id"
    metabolism['period_id'] = metabolism['period_start'].dt.strftime('%Y-%m-%d_%H')

    # Assign each measurement row to a period based on datetime range
    period_labels = []
    for _, row in measurements.iterrows():
        dt = row['datetime']

        # find matching metabolism row
        match = metabolism[(dt >= metabolism['period_start']) & (dt < metabolism['period_end'])]
        if not match.empty:
            period_labels.append(match.iloc[0]['period_id'])
        else:
            period_labels.append(None)

    measurements['period_id'] = period_labels

    # Merge measurement and metabolism data on "period_id"
    merged = pd.merge(measurements, metabolism, how='left', on='period_id')

    return merged

# ----- Predict DO Levels -----
def predict_do_variable_constants(data, n_simulations=1000):
    # Create results directory and name results file
    results_dir = os.path.join(os.path.dirname(__file__), 'Results')
    os.makedirs(results_dir, exist_ok=True)

    # Sort data by period and datetime
    data = data.sort_values(['period_id', 'datetime'])

    # Compute average water temperature for each period
    data['avg_water_temperature'] = data.groupby('period_id')['tempC'].transform('mean')

    # Get unique periods in chronological order
    unique_periods = data['period_id'].dropna().unique()

    # Ensure enough periods exist
    if len(unique_periods) < 2:
        print('Not enough periods of data.')
        return

    # Define output file paths
    predictions_path = os.path.join(results_dir, 'DO_predictions_montecarlo.csv')
    r2_path = os.path.join(results_dir, 'daily_r2_montecarlo.csv')

    if os.path.exists(predictions_path):
        os.remove(predictions_path)
    daily_r2_list = []

    # Loop through periods
    for i in range(1, len(unique_periods)):
        prev_period = unique_periods[i - 1]
        curr_period = unique_periods[i]

        prev_day = data[data['period_id'] == prev_period]
        curr_day = data[data['period_id'] == curr_period].copy()

        if prev_day.empty or curr_day.empty:
            continue

        # Get the last DO measurement from the previous day to use as the starting point for predictions
        last_row = prev_day.iloc[-1]
        last_DO_meas = last_row['DO.meas']

        # Get mean and sd of metablism parameters
        prev_K_mean = prev_day['K.mean'].iloc[0]
        prev_theta_mean = prev_day['theta.mean'].iloc[0]
        prev_R_mean = prev_day['R.mean'].iloc[0]
        prev_A_mean = prev_day['A.mean'].iloc[0]
        prev_p_mean = prev_day['p.mean'].iloc[0]

        prev_K_sd = prev_day['K.sd'].iloc[0]
        prev_theta_sd = prev_day['theta.sd'].iloc[0]
        prev_R_sd = prev_day['R.sd'].iloc[0]
        prev_A_sd = prev_day['A.sd'].iloc[0]
        prev_p_sd = prev_day['p.sd'].iloc[0]

        n_steps = len(curr_day)

        # Measurement noise standard deviation (assumed)
        tau_mean = 1000
        sigma = np.sqrt(1 / tau_mean)

        # Prepare input arrays for the simulation
        tempC_vals = curr_day['tempC'].values
        irr_vals = curr_day['I'].values
        atm_p_vals = curr_day['atmo.pressure'].values
        sal_vals = curr_day['salinity'].values
        avg_temp_vals = curr_day['avg_water_temperature'].values

        # Include last time measurement from previous day as the first input for the simulation
        last_temp = last_row['tempC']
        last_irr = last_row['I']
        last_press = last_row['atmo.pressure']
        last_sal = last_row['salinity']
        last_avg_temp = prev_day['avg_water_temperature'].iloc[0]

        tempC_input = np.concatenate(([last_temp], tempC_vals))
        irr_input = np.concatenate(([last_irr], irr_vals))
        atm_p_input = np.concatenate(([last_press], atm_p_vals))
        sal_input = np.concatenate(([last_sal], sal_vals))
        avgT_input = np.concatenate(([last_avg_temp], avg_temp_vals))

        # Sample constants for all simulations at once
        constants = np.zeros((n_simulations, 5))
        constants[:, 0] = np.random.normal(prev_K_mean, prev_K_sd, size=n_simulations) / (86400 / 600)
        constants[:, 1] = np.random.normal(prev_theta_mean, prev_theta_sd, size=n_simulations)
        constants[:, 2] = np.random.normal(prev_R_mean, prev_R_sd, size=n_simulations)
        constants[:, 3] = np.random.normal(prev_A_mean, prev_A_sd, size=n_simulations)
        constants[:, 4] = np.random.normal(prev_p_mean, prev_p_sd, size=n_simulations)

        # Run Monte Carlo simulations for the current day
        predictions, prev_DO_matrix = simulate_day_DO(
            n_steps=n_steps,
            n_sim=n_simulations,
            tempC=tempC_input,
            irr=irr_input,
            press=atm_p_input,
            salinity=sal_input,
            avg_temp=avgT_input,
            last_DO=last_DO_meas,
            constants=constants,
            sigma=sigma
        )

        # Store statistics of predictions
        curr_day['DO_predicted_median'] = np.nanmedian(predictions, axis=1)
        curr_day['DO_predicted_std'] = np.nanstd(predictions, axis=1)
        curr_day['prev_DO_used'] = np.nanmedian(prev_DO_matrix, axis=1)

        # Save predictions
        curr_day.to_csv(predictions_path, mode='a', header=not os.path.exists(predictions_path), index=False)

        # Filter to valid rows for comparison (non-NaN, non-zero)
        valid = curr_day.dropna(subset=['DO.meas', 'DO_predicted_median'])
        valid = valid[(valid['DO.meas'] != 0) & (valid['DO_predicted_median'] != 0)]
        if not valid.empty:
            # Compute RMSE
            rmse = root_mean_squared_error(valid['DO.meas'], valid['DO_predicted_median'])

            # --- Font size settings ---
            title_fontsize = 20
            label_fontsize = 16
            legend_fontsize = 16
            x_tick_fontsize = 16
            y_tick_fontsize = 16

            # Plot prediction vs measured DO for the current day
            plt.figure(figsize=(10, 6))
            plt.plot(valid['datetime'], valid['DO.meas'], label='Actual DO', color='red')
            plt.plot(valid['datetime'], valid['DO_predicted_median'], label='Predicted DO (median)', color='black', dashes=(6, 2))
            plt.fill_between(
                valid['datetime'],
                valid['DO_predicted_median'] - 1,
                valid['DO_predicted_median'] + 1,
                color='gray',
                alpha=0.3,
                label='±1 mg/L DO'
            )

            # Extract start datetime from curr_period and compute end datetime (+24 hours)
            period_start = datetime.strptime(curr_period, "%Y-%m-%d_%H")
            period_end = period_start + timedelta(hours=24)

            # Format the range string for the x-axis label
            time_range_str = f"{period_start.strftime('%Y-%m-%d %H:%M')} to {period_end.strftime('%Y-%m-%d %H:%M')}"

            # Also compute the next day for the plot title
            next_day_str = (period_start + timedelta(days=1)).strftime("%Y-%m-%d")

            # Plot formatting
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            plt.xlabel(f'Time of measurement ({time_range_str})', fontsize=label_fontsize)
            plt.ylabel('Dissolved Oxygen (mg/L)', fontsize=label_fontsize)
            plt.title(f'DO Prediction for {next_day_str}', fontsize=title_fontsize)
            plt.legend(fontsize=legend_fontsize)
            plt.xticks(rotation=45, fontsize=x_tick_fontsize)
            plt.yticks(fontsize=y_tick_fontsize)
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, f'DO_comparison_{curr_period}.png'), dpi=400, bbox_inches='tight', pad_inches=0.02)
            plt.close()

            # Compute R², MAE, and store results
            r2 = r2_score(valid['DO.meas'], valid['DO_predicted_median'])
            mae = mean_absolute_error(valid['DO.meas'], valid['DO_predicted_median'])
            n = len(valid)

            # Compute minimum DO values and their difference
            min_measured_do = np.min(curr_day["DO.meas"].values)
            min_predicted_do = np.min(curr_day["DO_predicted_median"].values)
            min_do_diff = min_measured_do - min_predicted_do
            
            daily_result = {
                'period_id': curr_period,
                'r2': r2,
                'rmse': rmse,
                'mae': mae,
                "Min_DO_diff": min_do_diff,
                "Min_measured_DO": min_measured_do,
                "Min_predicted_DO": min_predicted_do
            }
            daily_r2_list.append(daily_result)

            print(f'{curr_period} — R²: {r2:.4f}, RMSE: {rmse:.4f}, MAE: {mae:.4f}, n: {n}')

    results_df = pd.DataFrame(daily_r2_list)
    results_df.to_csv(r2_path, index=False)

# ----- Main Execution -----
def main():
    # Define file paths
    measurement_data = 'Measurement Datasets/Final_Format_P7T.csv'
    metabolism_data = 'Metabolism Datasets/Final_Format_P7T_Metabolism.csv'

    # Load and preprocess data
    data = load_data(measurement_data, metabolism_data)

    # Predict DO levels
    predict_do_variable_constants(data)
    
if __name__ == "__main__":
    main()