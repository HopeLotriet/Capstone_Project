create table if not exists lipinski_properties (
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
);