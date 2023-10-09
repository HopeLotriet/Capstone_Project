from flask import Flask, render_template, request, jsonify, send_from_directory
import joblib
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski
import traceback
import csv
import sqlite3
# from sklearn import preprocessing
# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler
# from sklearn.ensemble import RandomForestClassifier
# from sklearn import metrics
# from sklearn.preprocessing import LabelEncoder


app = Flask(__name__)


if not os.path.exists('lipinski_properties.db'):
    conn = sqlite3.connect('lipinski_properties.db')
    cursor = conn.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Compounds (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            Name TEXT,
            Formula TEXT,
            HBD REAL,
            HBA REAL,
            MolecularWeight REAL,
            LogP REAL,
            TPSA REAL,
            Num_Rotatable_Bonds REAL,
            SAS REAL,
            PassesLipinski BOOLEAN
        )
    ''')
    conn.commit()

    # Read data from the CSV file and insert it into the database
    with open('Capstone.csv', 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        for row in csv_reader:
            # Convert 'False' to 0 and 'True' to 1
            PassesLipinski = 1 if row['PassesLipinski'] == 'True' else 0

            cursor.execute('''INSERT INTO Compounds (Name, Formula, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, PassesLipinski) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', 
                           (row['Name'], row['Formula'], float(row['HBD']), float(row['HBA']), float(row['MolecularWeight']),
                            float(row['LogP']), float(row['TPSA']), float(row['Num_Rotatable_Bonds']),
                            float(row['SAS']), row['PassesLipinski']))
           
        conn.commit()

    # Close the connection outside the for loop
    conn.close()

else:
    conn = sqlite3.connect('lipinski_properties.db')
    cursor = conn.cursor()



 # If the database file already exists, just connect to it



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



@app.route('/calculate_properties', methods=['POST'])
def calculate_properties():
    data = request.get_json()
    chemical_formula = data.get('chemical_formula')
    molecular_weight = data.get('molecular_weight')

    if not chemical_formula and not molecular_weight:
        return jsonify(error='Invalid input. Please enter a valid chemical formula or molecular weight.'), 400

    try:
        if chemical_formula:
            mol = Chem.MolFromSmiles(chemical_formula)
        else:
            mol = Chem.MolFromSmiles('C')
            mol.SetProp('_MolWt', molecular_weight)

        if mol is None:
            return jsonify(error='Invalid input. Please enter a valid chemical formula or molecular weight.'), 400

        logp = Crippen.MolLogP(mol)
        hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
        hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
        tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
        sas_score = Chem.CalcSAS(mol)
        num_rotatable_bonds = Chem.CalcNumRotatableBonds(mol)

        is_lipinski_compliant = all([
            Lipinski.NumHDonors(mol) <= 5,
            Lipinski.NumHAcceptors(mol) <= 10,
            logp <= 5,
            mol.GetNumHeavyAtoms() <= 500
        ])

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

        return jsonify(result)

    except Exception as e:
        print(e)  # Log the specific error message
        return jsonify(error='Error occurred during calculation.'), 500
    
@app.route("/database")
def database():
    try:
        conn = sqlite3.connect('lipinski_properties.db')
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM Compounds")  # Replace 'your_table_name' with your actual table name
        data = cursor.fetchall()
        conn.close()  # Close the connection after fetching data
        return render_template("database.html", data=data)
    except Exception as e:
        print(e)
        return render_template("database.html", data=None)



if __name__ == "__main__":
    app.run(debug=True)
