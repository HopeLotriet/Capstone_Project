import csv
import sqlite3
from rdkit import Chem

# Connect to the SQLite database
conn = sqlite3.connect("lipinski_properties.db")
cursor = conn.cursor()

# Drop the table if it exists
cursor.execute('''
    DROP TABLE IF EXISTS compounds
''')

# Create a table to store the data
cursor.execute('''
    CREATE TABLE compounds (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        Name TEXT NOT NULL,
        SMILES TEXT,
        Formula TEXT,
        HBD REAL,
        HBA REAL,
        MolecularWeight REAL,
        LogP REAL,
        TPSA REAL,
        Num_Rotatable_Bonds REAL,
        PassesLipinski BOOLEAN,
        SAS REAL,
        Structure TEXT
    )
''')

# Read data from the CSV file and insert it into the database
with open('Capstone_Dataset.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)

    # Skip the header row
    next(csv_reader)

    for row in csv_reader:
        # Convert 'True' or 'False' strings to boolean values for the PassesLipinski column
        passes_lipinski = True if row[9].strip().lower() == 'true' else False

        # Parse the SMILES string to an RDKit molecule object
        smiles = row[1]
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            mol_str = Chem.MolToMolBlock(mol)
        else:
            mol_str = None

        # Insert data into the database
        cursor.execute('''
            INSERT INTO compounds (
                Name, SMILES, Formula, HBD, HBA, MolecularWeight,
                LogP, TPSA, Num_Rotatable_Bonds, PassesLipinski, SAS, Structure
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            row[0],  # Name
            smiles,  # SMILES
            row[2],  # Formula
            float(row[3]),  # HBD (converted to float)
            float(row[4]),  # HBA
            float(row[5]),  # MolecularWeight
            float(row[6]),  # LogP
            float(row[7]),  # TPSA
            float(row[8]),  # Num_Rotatable_Bonds
            passes_lipinski,  # PassesLipinski
            float(row[10]),  # SAS
            mol_str  # Structure
        ))

# Commit the changes and close the database connection
conn.commit()
conn.close()

print("Data imported successfully.")
