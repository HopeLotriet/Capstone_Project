import express from 'express';
import * as sqlite from 'sqlite';
import sqlite3 from 'sqlite3';
import fs from 'fs';
import csvParser from 'csv-parser';

const app = express();

app.use(express.static('public'))
app.use(express.json());

async function insertDataFromCSV() {
    const db = await sqlite.open({
        filename: './lipinski_properties.db',
        driver: sqlite3.Database
    });

    try {
        await db.migrate(); // Perform migrations

        fs.createReadStream('Capstone_Dataset.csv')
            .pipe(csvParser())
            .on('data', async (row) => {
                try {
                    // Insert data into the lipinski_properties table
                    await db.run('INSERT INTO lipinski_properties (id, Name, SMILES, Formula, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, PassesLipinski, SAS, Structure) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                        row.id, row.Name, row.SMILES, row.Formula, row.HBD, row.HBA, row.MolecularWeight, row.LogP, row.TPSA, row.Num_Rotatable_Bonds, row.PassesLipinski, row.SAS, row.Structure);
                } catch (error) {
                    console.error(error);
                }
            })
            .on('end', async () => {
               // console.log('Data inserted successfully.');
                await db.close(); // Close the database connection after the migration and data insertion
            });
    } catch (error) {
        console.error('Error during migration:', error);
        await db.close(); // Ensure the database connection is closed in case of an error
    }
}

insertDataFromCSV();



app.get('/get_compound', async function (req, res) {

    try {
        const getCompound = await getAllCompounds();
        res.json({
            getCompound: getCompound
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve compound Id.' });
    }
});

export async function getAllCompounds() {
    const db = await sqlite.open({
        filename: './lipinski_properties.db',
        driver: sqlite3.Database
    });
    const result = await db.all(`select * from lipinski_properties`);
    await db.close(); // Close the database connection after retrieving data
    return result;
}

// export async function getAllCompounds() {
//     const result = await db.all(`select * from lipinski_properties`);
//     return result;
// }

///////////

app.get('/get_logp', async function (req, res) {

    try {
        const getLogP = await getLogPs();
        res.json({
            getLogP: getLogP
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve LogP for compound.' });
    }
});

export async function getLogPs() {
    const result = await db.all(`select * from lipinski_properties`);
    return result;
}
////////////

app.post('/predict_lipinski', async function (req, res) {
    try {
        const Name = req.body.Name
        const HBA = req.body.HBA;
        const HBD = req.body.HBD;
        const MolecularWeight = req.body.MolecularWeight;
        const LogP = req.body.LogP;
        const TPSA = req.body.TPSA;
        const Num_Rotatable_Bonds = req.body.Num_Rotatable_Bonds;
        const SAS = req.body.SAS;
        const passesLipinski = req.body.passesLipinski;

        await predictLipinski(Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski);

        res.json({
            status: 'Sucess',
            message: `The Lipinski Properties for ${Name} is ${passesLipinski}`
        })

    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to predict the Lipinski Properties.' });
    }
});

export async function predictLipinski(Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski) {
    const sql = `predict_lipinski lipinski_properties set Name = ?, HBD = ?, HBA = ?, MolecularWeight = ?, LogP = ?, TPSA = ?, Num_Rotatable_Bonds = ?, SAS = ? where passesLipinski = ?`;
    await db.run(sql, [Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski])
}

////////

app.post('/calculate_molecular_weight', async function (req, res) {
    try {
        const Name = req.body.Name;
        const MolecularWeight = req.body.MolecularWeight;

        await calculateMolecularWeight(Name, MolecularWeight);

        res.json({
            status: 'Sucess',
            message: `Calculated Molecular Weight for ${Name}`
        })

    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to calculate the Molecular Weight.' });
    }
});

export async function calculateMolecularWeight(Name, MolecularWeight) {
    const sql = `calculate_molecular_weight from lipinski_properties set MolecularWeight = ?, where Name = ?`;
    await db.run(sql, [Name, MolecularWeight])
}

///////

app.post('/search_compound', async function (req, res) {
    try {
        const id = req.body.id;
        const Name = req.body.Name;

        await searchCompound(id, Name);

        res.json({
            status: 'Sucess',
            message: `Compound Search successful for ${Name}`
        })

    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Search Failed.' });
    }
});

export async function searchCompound(id, Name) {
    const sql = `search_compound from lipinski_properties set Name = ?, where id = ?`;
    await db.run(sql, [id, Name])
}



const PORT = process.env.PORT || 5000;
app.listen(PORT, () => console.log(`Server started on port: ${PORT}`))