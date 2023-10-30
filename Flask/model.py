from flask import Flask, render_template, request, jsonify, g
# import requests
from flask_cors import CORS
import joblib
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
import traceback
import csv
import sqlite3
import logging


app = Flask(__name__)
CORS(app)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect('lipinski_properties.db')
    return db

conn = sqlite3.connect('lipinski_properties.db', check_same_thread=False)
cursor = conn.cursor()

    # cursor.execute('''
    #     CREATE TABLE IF NOT EXISTS Compounds (
    #         id INTEGER PRIMARY KEY AUTOINCREMENT,
    #         Name TEXT,
    #         Formula TEXT,
    #         HBD REAL,
    #         HBA REAL,
    #         MolecularWeight REAL,
    #         LogP REAL,
    #         TPSA REAL,
    #         Num_Rotatable_Bonds REAL,
    #         SAS REAL,
    #         PassesLipinski BOOLEAN
    #     )
    # ''')
    # conn.commit()

    # Read data from the CSV file and insert it into the database
with open('Capstone.csv', 'r') as csvfile:
    csv_reader = csv.DictReader(csvfile)
    for row in csv_reader:
            # Convert 'False' to 0 and 'True' to 1
        PassesLipinski = 1 if row['PassesLipinski'] == 'True' else 0
        cursor.execute('SELECT * FROM Compounds WHERE Name = ?', (row['Name'],))
        existing_record = cursor.fetchone()
        if existing_record is None: 

         cursor.execute('''INSERT INTO Compounds (Name, Formula, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, PassesLipinski) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (row['Name'], row['Formula'], float(row['HBD']), float(row['HBA']), float(row['MolecularWeight']),
                float(row['LogP']), float(row['TPSA']), float(
                 row['Num_Rotatable_Bonds']),
                float(row['SAS']), row['PassesLipinski']))

    conn.commit()

    # Close the connection outside the for loop
    # conn.close()

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


@app.route("/search")
def search():
    return render_template("search.html")

# @app.route("/database")
# def database():
#     return render_template("database.html")

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



# Set up logging configuration
# logging.basicConfig(level=logging.ERROR)  # Set the logging level as needed


import requests
from flask import jsonify, request

@app.route('/calculate_properties', methods=['POST'])
def calculate_properties():
    result = {}

    try:
        data = request.get_json()
        chemical_formula = data.get('chemical_formula')
        molecular_weight = data.get('molecular_weight')

        # Input validation
        if not chemical_formula and not molecular_weight:
            result['error'] = 'Invalid input. Please enter a valid chemical formula or molecular weight.'
            return jsonify(result), 400

        if chemical_formula:
            # Convert chemical formula to SMILES using CIR
            smiles = formula_to_smiles(chemical_formula)
            if smiles is None:
                result['error'] = f'Invalid chemical formula: {chemical_formula}'
                return jsonify(result), 400

            mol = MolFromSmiles(smiles)
        else:
            mol = MolFromSmiles('C')
            mol.SetProp('_MolWt', molecular_weight)

        # Rest of your code for property calculations...
        # (logp, hbd, hba, tpsa, sas_score, num_rotatable_bonds, is_lipinski_compliant)

        result = {
            'molecular_weight': mol.GetProp('_MolWt'),
            'logP': logp,
            'hydrogen_bond_donors': hbd,
            'hydrogen_bond_acceptors': hba,
            'tpsa': tpsa,
            'sas_score': sas_score,
            'num_rotatable_bonds': num_rotatable_bonds,
            'is_lipinski_compliant': is_lipinski_compliant
        }

    except Exception as e:
        # Handle exceptions appropriately
        result['error'] = 'Error occurred during calculation: ' + str(e)
        return jsonify(result), 400

    return jsonify(result)

def formula_to_smiles(formula):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{formula}/smiles"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip('"')
    else:
        return None


#....

def format_data(data):
    formatted_data = []
    for row in data:
        formatted_row = []
        for item in row:
            if isinstance(item, (int, float)):
                formatted_row.append(round(item, 2))
            else:
                formatted_row.append(item)
        formatted_data.append(formatted_row)
    return formatted_data


@app.route("/database")
def database():
    
        conn = sqlite3.connect('lipinski_properties.db')
        cursor = conn.cursor()
        # Replace 'your_table_name' with your actual table name
        cursor.execute("SELECT * FROM Compounds")
        data = cursor.fetchall()
        conn.close()  # Close the connection after fetching data

        # formatted_data = [(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]) for row in data]
    
    # Pass the formatted data to the template
        # print(formatted_data)
        formatted_data = format_data(data)
        return render_template("database.html", data=formatted_data)
        # return render_template("database.html", data=data)
#     except Exception as e:
#         print(e)
#         return render_template("database.html", data=None)


@app.post('/search')
def search_database():
    print('SEARCHING...')
    query = request.args.get('query')
    # Get the search query from the request parameters

    # Execute the SQL query with placeholders for partial matching in Name or Formula
    cursor.execute("SELECT * FROM Compounds WHERE Name LIKE ? OR Formula LIKE ?",
                   ('%' + query + '%', '%' + query + '%'))
    
    search_results = cursor.fetchall()
    # columns = [desc[0] for desc in cursor.description]
    # results = [dict(zip(columns, row)) for row in search_results]
    print(search_results)
    return jsonify(search_results)


# .....

    # try:
    #     query = request.args.get('query')
    #     # Get the search query from the request parameters

    #     # Execute the SQL query with placeholders for partial matching in Name or Formula
    #     cursor.execute("SELECT * FROM Compounds WHERE Name LIKE ? OR Formula LIKE ?",
    #                    ('%' + query + '%', '%' + query + '%'))

    #     # Fetch the search results
    #     search_results = cursor.fetchall()

    #     # Return search results as JSON with HTTP status 200 (OK)

    #     return jsonify(search_results), 200

    # except Exception as e:
    # Handle errors and return an error response as JSON with HTTP status 500 (Internal Server Error)
    # error_message = str(e)  # Get the error message as a string
    # print("Error:", error_message)  # Log the error for debugging purposes
    # return jsonify(error="An error occurred while searching the database."), 500

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

if __name__ == "__main__":
    app.run(debug=True)
