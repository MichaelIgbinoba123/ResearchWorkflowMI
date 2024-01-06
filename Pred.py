import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from keras.models import load_model

# Load the trained model
model = load_model('best_model.h5')  # Replace with the actual path to your trained model file

# Load the scaler used for normalization
scaler = StandardScaler()
scaler.mean_ = np.array([...])  # Replace with the mean values used for scaling during training
scaler.scale_ = np.array([...])  # Replace with the scale values used for scaling during training

# Prompt for user input
print("Please provide input values for the following variables:")
Latitude = float(input("Latitude (deg): "))
Eye_Temp = float(input("Eye Temp (C): "))
Eye_Diameter = float(input("Eye Diameter (km): "))
Symmetry = float(input("Symmetry: "))
WV_Temp = float(input("WV Temp (C): "))
Eyewall_Temp = float(input("Eyewall Temperature (C): "))
Eye_Diameter_CDO_Ratio = float(input("Eye Diameter/CDO Diameter Ratio: "))

# Create an input array for prediction
input_data = np.array([[Latitude, Eye_Temp, Eye_Diameter, Symmetry, WV_Temp, Eyewall_Temp, Eye_Diameter_CDO_Ratio]])

# Normalize input data using the saved scaler
normalized_input_data = scaler.transform(input_data)

# Make a prediction
predicted_pressure_deficit = model.predict(normalized_input_data)

print(f"Predicted Pressure Deficit (hPa): {predicted_pressure_deficit[0][0]}")