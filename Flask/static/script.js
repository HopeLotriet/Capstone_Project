function init() {
    document.addEventListener('alpine:init', () => {
        Alpine.data('predictApp', () => {
            return {
                searchQuery: '',
                chemicalFormula: '',
                molecularWeight: '',
                searchResults: [],
                predictionResult: null,
                async predictProperties() {
                    try {
                        const response = await fetch('/predict_lipinski', {
                            method: 'POST',
                            headers: {
                                'Content-Type': 'application/json'
                            },
                            body: JSON.stringify({
                                Name: this.compoundName,
                                // Add other input fields here if needed
                            })
                        });

                        const data = await response.json();
                        this.predictionResult = data; // Update the reactive variable
                    } catch (error) {
                        console.error('Error:', error);
                        // Handle error
                    }
                },
                async calculateProperties() {
                    const data = {};
                    if (this.chemicalFormula) {
                        data.chemical_formula = this.chemicalFormula;
                    } else if (!isNaN(this.molecularWeight) && this.molecularWeight > 0) {
                        data.molecular_weight = this.molecularWeight;
                    } else {
                        this.resultDiv = 'Invalid input. Please enter a valid chemical formula or molecular weight.';
                        return;
                    }

                    try {
                        const response = await fetch('/calculate_properties', {
                            method: 'POST',
                            headers: {
                                'Content-Type': 'application/json'
                            },
                            body: JSON.stringify(data)
                        });

                        if (response.status === 200) {
                            const responseData = await response.json();
                            this.resultDiv = `Molecular Weight: ${responseData.molecular_weight} g/mol\nLogP: ${responseData.logP}\nHydrogen Bond Donors: ${responseData.hydrogen_bond_donors}\nHydrogen Bond Acceptors: ${responseData.hydrogen_bond_acceptors}\nTPSA: ${responseData.tpsa} Å²\nSAS Score: ${responseData.sas_score}\nNumber of Rotatable Bonds: ${responseData.num_rotatable_bonds}\n\nLipinski's Rule of Five Compliance: ${responseData.is_lipinski_compliant ? 'Yes' : 'No'}`;
                        } else {
                            this.resultDiv = 'Error occurred during calculation.';
                        }
                    } catch (error) {
                        console.error('Error:', error);
                        this.resultDiv = 'Error occurred during calculation.';
                    }
                },
                async searchDatabase() {
                    this.searchResults = []; // Clear previous search results

                    try {
                        const response = await fetch(`/search?query=${this.searchQuery}`, {
                            method: 'GET',
                            headers: {
                                'Content-Type': 'application/json'
                            }
                        });

                        if (response.status === 200) {
                            const searchData = await response.json();
                            if (searchData.length > 0) {
                                // Render search results in the searchResults array
                                this.searchResults = searchData;
                            } else {
                                this.searchResults.push({ id: null, Name: 'No results found.', Formula: null });
                            }
                        } else {
                            this.searchResults.push({ id: null, Name: 'Error occurred during search.', Formula: null });
                        }
                    } catch (error) {
                        console.error('Error:', error);
                        this.searchResults.push({ id: null, Name: 'Error occurred during search.', Formula: null });
                    }
                }
            };
        });
    });
}

init();
