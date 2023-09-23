from flask import Flask, render_template, request, jsonify
import joblib
import pandas as pd #For data preprocessing

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


#Load the pre-trained Random Forest Model
model = joblib.load('Random_Forest_Model.pkl')

# Sample feature names (adjust to match your model's feature names)
feature_names = ['HBD', 'HBA', 'MolecularWeight', 'LogP', 'TPSA', 'Num_Rotatable_Bonds', 'SAS']


@app.route('/predict_lipinski', methods=['POST'])
def predict_lipinski():
   try:
        # Get the input data from the request
        input_data = request.json

        # Ensure that the input data contains all required features
        if not all(feature in input_data for feature in feature_names):
            return jsonify({'error': 'Missing or incorrect input features'}), 400

        # Preprocess the input data if needed (e.g., convert to DataFrame)
        input_df = pd.DataFrame([input_data], columns=feature_names)

        # Use your Random Forest Model for prediction
        prediction = model.predict(input_df)

        # Interpret the prediction (0 for False, 1 for True, if 'PassesLipinski' is a binary label)
        prediction_result = 'Passes Lipinski' if prediction[0] == 1 else 'Fails Lipinski'

        # Return the prediction as JSON
        return jsonify({'prediction': prediction_result})
   except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)


