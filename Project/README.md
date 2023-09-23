Information based on the project

app.py: This is the main Flask application file where you define your API endpoints, routes, and interact with your Random Forest Model.

random_forest_model.pkl: This is your pre-trained Random Forest Model saved as a pickle file.

import_csv_to_db.py: This script is used to import your dataset into the SQLite database, as explained in previous messages. You can run this script once to populate your database.

templates/: This folder can contain HTML templates if you are building a web interface for your application using Flask. You can create HTML files for different pages of your web app.

static/: This folder can contain static assets such as CSS files, JavaScript files, images, etc., used in your web interface.


# Send a POST request to predict Lipinski's properties
curl -X POST -H "Content-Type: application/json" -d '{"HBD": 2, "HBA": 3, "MolecularWeight": 200, "LogP": 1.5, "TPSA": 50, "Num_Rotatable_Bonds": 4, "SAS": 80}' http://localhost:5000/predict_lipinski
