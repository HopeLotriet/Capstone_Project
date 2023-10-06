import express from 'express';
import * as sqlite from 'sqlite';
import sqlite3 from 'sqlite3';
import fs from 'fs';
import csvParser from 'csv-parser';

console.log("Current Working Directory: ", process.cwd());

const app = express();

app.use(express.static('public'));
app.use(express.json());

async function initializeDatabase() {
    try {
        const db = await sqlite.open({
            filename: './lipinski_properties.db',
            driver: sqlite3.Database
        });
        await db.migrate(); // Perform migrations
        await insertDataFromCSV(db);
        return db;
    } catch (error) {
        console.error('Error initializing database:', error);
        throw error;
    }
}

// async function insertDataFromCSV(db) {
//     const data = [];
//     return new Promise((resolve, reject) => {
//         fs.createReadStream('Project/Express API/Capstone_Dataset.csv')
//             .pipe(csvParser())
//             .on('data', (row) => {
//                 data.push(row);
//             })
//             .on('end', async () => {
//                 try {
//                     await db.run('BEGIN');
//                     for (const row of data) {
//                         await db.run('INSERT INTO lipinski_properties (id, Name, SMILES, Formula, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, PassesLipinski, SAS, Structure) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
//                             row.id, row.Name, row.SMILES, row.Formula, row.HBD, row.HBA, row.MolecularWeight, row.LogP, row.TPSA, row.Num_Rotatable_Bonds, row.PassesLipinski, row.SAS, row.Structure);
//                     }
//                     await db.run('COMMIT');
//                     console.log('Data inserted successfully.');
//                     resolve();
//                 } catch (error) {
//                     await db.run('ROLLBACK');
//                     console.error('Error inserting data:', error);
//                     reject(error);
//                 }
//             });
//     });
// }

let dbInstance; // Database connection instance

(async () => {
    try {
        dbInstance = await initializeDatabase();
        console.log('Database initialized successfully.');
    } catch (error) {
        console.error('Failed to initialize database:', error);
    }
})();

// Error handling middleware
app.use((err, req, res, next) => {
    console.error(err.stack);
    res.status(500).json({ error: 'Something went wrong!' });
});



app.get('/get_compound', async function (req, res) {
    try {
        const getCompound = await getAllCompounds(dbInstance);
        res.json({
            getCompound: getCompound
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve compound data.' });
    }
});

// export async function getAllCompounds() {
//     const db = await sqlite.open({
//         filename: './lipinski_properties.db',
//         driver: sqlite3.Database
//     });
//     const result = await db.all(`select * from lipinski_properties`);
//     await db.close(); // Close the database connection after retrieving data
//     return result;
// }

// export async function getAllCompounds() {
//     const result = await db.all(`select * from lipinski_properties`);
//     return result;
// }

///////////

app.get('/get_logp', async function (req, res) {
    try {
        const getLogP = await getLogPs(dbInstance);
        res.json({
            getLogP: getLogP
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve LogP data.' });
    }
});

// export async function getLogPs() {
//     const result = await db.all(`select * from lipinski_properties`);
//     return result;
// }
////////////

app.post('/predict_lipinski', async function (req, res) {
    try {
        const { Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski } = req.body;
        await predictLipinski(dbInstance, Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski);
        res.json({
            status: 'Success',
            message: `The Lipinski Properties for ${Name} have been updated.`
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to predict the Lipinski Properties.' });
    }
});


// export async function predictLipinski(Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski) {
//     const sql = `predict_lipinski lipinski_properties set Name = ?, HBD = ?, HBA = ?, MolecularWeight = ?, LogP = ?, TPSA = ?, Num_Rotatable_Bonds = ?, SAS = ? where passesLipinski = ?`;
//     await db.run(sql, [Name, HBD, HBA, MolecularWeight, LogP, TPSA, Num_Rotatable_Bonds, SAS, passesLipinski])
// }

////////

app.post('/calculate_molecular_weight', async function (req, res) {
    try {
        const { Name, MolecularWeight } = req.body;
        await calculateMolecularWeight(dbInstance, Name, MolecularWeight);
        res.json({
            status: 'Success',
            message: `Calculated Molecular Weight for ${Name}.`
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to calculate the Molecular Weight.' });
    }
});

// export async function calculateMolecularWeight(Name, MolecularWeight) {
//     const sql = `calculate_molecular_weight from lipinski_properties set MolecularWeight = ?, where Name = ?`;
//     await db.run(sql, [Name, MolecularWeight])
// }

///////

app.post('/search_compound', async function (req, res) {
    try {
        const { id, Name } = req.body;
        await searchCompound(dbInstance, id, Name);
        res.json({
            status: 'Success',
            message: `Compound Search successful for ${Name}.`
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Search Failed.' });
    }
});

// export async function searchCompound(id, Name) {
//     const sql = `search_compound from lipinski_properties set Name = ?, where id = ?`;
//     await db.run(sql, [id, Name])
// }

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => console.log(`Server started on port: ${PORT}`))