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
from rdkit import Chem
from edkit.Chem import Draw
import sqlite3  # Import SQLite library

app = Flask(__name__)

# Establish a connection to the SQLite database
db_connection = sqlite3.connect("lipinski_properties.db")


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

#data preprocessing
def preprocess_data(data):

    df = pd.DataFrame(data)

    le = preprocessing.LabelEncoder()
    df['PassesLipinski'] = le.fit_transform(df['PassesLipinski'])

#Encode the 6 properties

# Create a LabelEncoder object
    label_encoder = LabelEncoder()

# Fit and transform the labels
    for column in ['HBD', 'HBA', 'MolecularWeight', 'LogP', 'TPSA', 'Num_Rotatable_Bonds', 'SAS']:
        df[column] = label_encoder.fit_transform(df[column])

    X = df[['HBD', 'HBA', 'MolecularWeight','LogP','TPSA','Num_Rotatable_Bonds','SAS']] #Features
    y = df['PassesLipinski'] #Target variable

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, shuffle=True, stratify=y)

    print("Train data shape:", X_train.shape)
    print("Train labels shape:", y_train.shape)
    print("Test data shape:", X_test.shape)
    print("Test labels shape:", y_test.shape)

    return df



@app.route('/predict_lipinski', methods=['POST'])
def predict_lipinski():
   try:
       # Extract data from the POST request and preprocess it
       data = request.get_json()
       preprocessed_data = preprocess_data(data)

       # Check if the model is fitted
       if hasattr(model, "classes_"):
           # Use the loaded and fitted model to make predictions
           predictions = model.predict(preprocessed_data)

           # Return the predictions as a JSON response
           return jsonify(predictions.tolist())
       else:
           return jsonify({"error": "Model is not fitted. Call 'fit' before making predictions."})

   except Exception as e:
        return jsonify({"error": str(e)})



# API endpoint to retrieve compound information by ID
@app.route('/get_compound/<int:compound_id>', methods=['GET'])
def get_compound(compound_id):
    try:
        # Execute a SQL query to fetch compound data by ID
        cursor = db_connection.cursor()
        cursor.execute("SELECT * FROM compounds WHERE id = ?", (compound_id,))
        compound_data = cursor.fetchone()

        if compound_data:
            # Convert the compound data to a dictionary and return it as JSON
            compound_info = {
                "id": compound_data[0],
                "Name": compound_data[1],
                "SMILES": compound_data[2],
                "Formula": compound_data[3],
                "HBD": compound_data[4],
                "HBA": compound_data[5],
                "MolecularWeight": compound_data[6],
                "LogP": compound_data[7],
                "TPSA": compound_data[8],
                "Num_Rotatable_Bonds": compound_data[9],
                "PassesLipinski": compound_data[10],
                "SAS": compound_data[11],
                "Structure": compound_data[12],
                
            }
            return jsonify(compound_info)
        else:
            return jsonify({"error": "Compound not found."}), 404

    except Exception as e:
        return jsonify({"error": str(e)})
    
    
@app.route('/get_logp/<int:compound_id>', methods=['GET'])
def get_logp(compound_id):
    try:
        # Execute a SQL query to fetch the LogP property of a compound by ID
        cursor = db_connection.cursor()
        cursor.execute("SELECT LogP FROM compounds WHERE id = ?", (compound_id,))
        logp = cursor.fetchone()

        if logp:
            return jsonify({"LogP": logp[0]})
        else:
            return jsonify({"error": "Compound not found."}), 404

    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/calculate_properties', methods=['POST'])
def calculate_properties():
    try:
        # Extract data from the POST request and preprocess it
        data = request.get_json()
        preprocessed_data = preprocess_data(data)

        # Calculate the 6 feature properties using the same preprocessing steps
        calculated_properties = preprocessed_data.mean().to_dict()

        return jsonify(calculated_properties)

    except Exception as e:
        return jsonify({"error": str(e)})


# Define a function for substructure search
def substructure_search(smiles_query):
    try:
        # Initialize a list to store matching compound data
        matching_compounds = []

        # Execute a SQL query to fetch compounds by SMILES and Name
        cursor = db_connection.cursor()
        cursor.execute("SELECT * FROM compounds WHERE SMILES LIKE ? OR Name LIKE ?",
                       ('%' + smiles_query + '%', '%' + smiles_query + '%'))

        for row in cursor.fetchall():
            compound_data = {
                "id": row[0],
                "Name": row[1],
                "SMILES": row[2],
                "Formula": row[3],
                "HBD": row[4],
                "HBA": row[5],
                "MolecularWeight": row[6],
                "LogP": row[7],
                "TPSA": row[8],
                "Num_Rotatable_Bonds": row[9],
                "PassesLipinski": row[10],
                "SAS": row[11],
                "Structure": row[12],
            }
            matching_compounds.append(compound_data)

        return matching_compounds

    except Exception as e:
        return {"error": str(e)}

@app.route('/search_compounds', methods=['GET'])
def search_compounds():
    try:
        # Get the search query from the request
        search_query = request.args.get('query')

        # Perform substructure search using RDKit
        matching_compounds = substructure_search(search_query)

        # Return the matching compounds as JSON
        return jsonify(matching_compounds)

    except Exception as e:
        return jsonify({"error": str(e)})

if __name__ == '__main__':
    app.run(debug=True)




