import express from 'express';
import * as sqlite from 'sqlite';
import sqlite3 from 'sqlite3';

const app = express ();

app.use(express.static('public'))

app.use(express.json());

const  db = await sqlite.open({
    filename:  './lipinski_properties.db',
    driver:  sqlite3.Database
});

await db.migrate()

app.get('/get_compound', async function (req, res) {

    try {
        const getCompound = await getAllCompounds();
        res.json ({
            getCompound : getCompound
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve compound Id.'});
    }
});

export async function getAllCompounds () {
    const result = await db.all(`select * from lipinski_properties`);
    return result;
}

app.get('/get_logp', async function (req, res) {

    try {
        const getLogP = await getLogPs();
        res.json ({
            getLogP : getLogP
        });
    } catch (error) {
        console.error(error);
        res.status(500).json({ error: 'Failed to retrieve compound Id.'});
    }
});

export async function getLogPs () {
    const result = await db.all(`select * from lipinski_properties`);
    return result;
}



const PORT = process.env.PORT || 5000;
app.listen(PORT, () => console.log (`Server started on port: ${PORT}`))