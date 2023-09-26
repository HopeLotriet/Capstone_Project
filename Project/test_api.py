import requests

data = [{
    "HBD": 2,
    "HBA": 3,
    "MolecularWeight": 200,
    "LogP": 1.5,
    "TPSA": 50,
    "Num_Rotatable_Bonds": 4,
    "SAS": 80
}]

# Send a POST request to predict Lipinski's properties
response = requests.post("http://localhost:5000/predict_lipinski", json=data)
print(response)
if response.status_code == 200:
    result = response.json()
    print("Prediction:", result.get("prediction"))
else:
    print("Error:", response.status_code)
