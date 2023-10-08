from flask import Flask, render_template, request, jsonify, send_from_directory
import joblib
import os
import pandas as pd
# from sklearn import preprocessing
# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler
# from sklearn.ensemble import RandomForestClassifier
# from sklearn import metrics
# from sklearn.preprocessing import LabelEncoder


app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/predict")
def predict():
    return render_template("predict.html")

@app.route("/classification")
def classification():
    return render_template("classification.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

@app.route("/guideline")
def guideline():
    return render_template("guideline.html")

# @app.route('/')
# def serve_index():
#     return send_from_directory('public', 'index.html')

# Get the absolute path to the model file
model_file_path = os.path.abspath('Best_RFModel.pkl')

# Load the pre-trained model
model = joblib.load(model_file_path)



@app.route('/predict', methods=['POST'])
def check_rule():
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

@app.route('/classify', methods=['POST'])
def classify_drug():
    # Get the input data from the form
    input_features = request.form.to_dict()
    
    # Convert input features to a DataFrame
    data = pd.DataFrame([input_features])
    
    # Make predictions
    prediction = model.predict(data)
    
    # Determine if the molecular properties pass the Lipinski Rule
    properties_passed = sum(prediction)
    total_properties = len(prediction)
    pass_percentage = (properties_passed / total_properties) * 100
    
    # Classify the drug based on pass percentage
    if pass_percentage >= 90: 
        classification = "Excellent Drug"
    elif pass_percentage >= 70: 
        classification = "Good Drug"
    elif pass_percentage >= 50: 
        classification = "Moderate Drug"
    elif pass_percentage >= 20: 
        classification = "Bad Drug"
    else: 
        classification = "Very Bad Drug"
    
    # Return the prediction result and classification as JSON
    return jsonify({'pass_percentage': pass_percentage, 'classification': classification})

if __name__ == '__main__':
    app.run(debug=True)
