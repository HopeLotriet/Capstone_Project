from flask import Flask, render_template, request, jsonify
import joblib
import os
import pandas as pd

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/search")
def search():
    return render_template("search.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

# Get the absolute path to the model file
model_file_path = os.path.abspath('RFModel.pkl')

# Load the pre-trained model
model = joblib.load(model_file_path)

#Load the pre-trained Random Forest Model
#model = joblib.load('RandomForestModel.pkl')
# print(model)

# Sample feature names (adjust to match your model's feature names)
feature_names = ['HBD', 'HBA', 'MolecularWeight', 'LogP', 'TPSA', 'Num_Rotatable_Bonds', 'SAS']


@app.route('/predict_lipinski', methods=['POST'])
def predict_lipinski():

    # Extract data from the POST request and preprocess it as needed
    data = request.get_json()

    # Use the loaded model to make predictions
    predictions = model.predict(data)

    # Return the predictions as a JSON response
    return jsonify(predictions.tolist())


# @app.route('/predict_lipinski', methods=['POST'])
# def predict_lipinski():
    
#    try:
#         # Get the input data from the request
#         input_data = request.json

#          # Debugging: Print the received input data
#         print("Received input data:", input_data)        

#         # Ensure that the input data contains all required features
#         if not all(feature in input_data for feature in feature_names):
#             return jsonify({'error': 'Missing or incorrect input features'}), 400

#         # Preprocess the input data if needed (e.g., convert to DataFrame)
#         input_df = pd.DataFrame([input_data], columns=feature_names)

#         # Debugging: Print the input DataFrame
#         print("Input DataFrame:", input_df)

#         # Use your Random Forest Model for prediction
#         prediction = model.predict(input_df)

#          # Debugging: Print the prediction
#         print("Prediction:", prediction)

#         # Interpret the prediction (0 for False, 1 for True, if 'PassesLipinski' is a binary label)
#         prediction_result = 'Passes Lipinski' if prediction[0] == 1 else 'Fails Lipinski'

#         # Return the prediction as JSON
#         return jsonify({'prediction': prediction_result})
#    except Exception as e:
#         return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)


