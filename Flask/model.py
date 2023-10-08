from flask import Flask, render_template, request, jsonify
import joblib
import os
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.preprocessing import LabelEncoder


app = Flask(__name__)



# @app.route("/")
# def index():
#     return render_template("index.html")


@app.route("/search")
def search():
    return render_template("search.html")



# Get the absolute path to the model file
model_file_path = os.path.abspath('Best_RFModel.pkl')

# Load the pre-trained model
model = joblib.load(model_file_path)



@app.route('/predict', methods=['POST'])
def predict():
    # Get the input data from the form
    input_features = request.form.to_dict()
    
    # Convert input features to a DataFrame
    data = pd.DataFrame([input_features])
    
    # Make predictions
    prediction = model.predict(data)
    
    # Determine if the molecular properties pass the Lipinski Rule
    result = "Pass" if prediction[0] == 1 else "Fail"
    
    # Render the result in the HTML file
    return jsonify({'prediction': result})

if __name__ == '__main__':
    app.run(debug=True)
