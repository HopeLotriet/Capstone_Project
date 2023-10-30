function init() {
    document.addEventListener('alpine:init', () => {
        Alpine.data('predictApp', () => {
            return {
                searchQuery: '',
                chemicalFormula: '',
                molecularWeight: '',
                searchResults: [],
                predictionResult: null,

             
async calculateProperties() {
    const data = {
        chemical_formula: this.chemicalFormula,
        molecular_weight: this.molecularWeight
    };

    try {
        const response = await axios.post('/calculate_properties', data, {
            headers: {
                'Content-Type': 'application/json'
            }
        });

        console.log(response);
        if (response.status === 200) {
            const responseData = response.data;
            this.resultDiv = `Molecular Weight: ${responseData.molecular_weight} g/mol\nLogP: ${responseData.logP}\n...`; // Update UI with the response data
        } else {
            console.error('Error occurred during calculation:', response.statusText);
            this.resultDiv = 'Error occurred during calculation.';
        }
    } catch (error) {
        console.error('Error:', error);
        this.resultDiv = 'Error occurred during calculation.';
    }
},


                async searchDatabase() {
                    this.searchResults = []; // Clear previous search results

                    axios.post(`/search?query=${this.searchQuery}`).then(res=>{
                        console.log(res.data);
                        this.searchResults = res.data;

                         // Toggle visibility of the table based on search results
                         let searchResults = document.getElementById('searchResults');
                         searchResults.style.display = this.searchResults.length > 0 ? 'table' : 'none';
                     }).catch(error => {
                         console.error('Error occurred during search:', error);
                         this.searchResults.push({ id: null, Name: 'Error occurred during search.', Formula: null });
 
                         // Hide the table if there are no search results
                         let searchResults = document.getElementById('searchResults');
                         searchResults.style.display = 'none';
                    });
                    // try {
                    //     const response = await fetch(`/search?query=${this.searchQuery}`, {
                    //         method: 'POST',
                    //         headers: {
                    //             'Content-Type': 'application/json'
                    //         },
                    //         body: {
                    //             query: this.searchQuery
                    //         }
                    //     });

                    //     if (response.ok) {
                    //         const searchData = await response.json(); // Parse the response as JSON
                    //         if (searchData.length > 0) {
                    //             this.searchResults = searchData;
                    //         } else {
                    //             this.searchResults.push({ id: null, Name: 'No results found.', Formula: null });
                    //         }
                    //     } else {
                    //         console.error('Error occurred during search:', response.statusText);
                    //         this.searchResults.push({ id: null, Name: 'Error occurred during search.', Formula: null });
                    //     }
                    // } catch (error) {
                    //     console.error('Error:', error);
                    //     this.searchResults.push({ id: null, Name: 'Error occurred during search.', Formula: null });
                    // }
                }
            };
        });
    });
}

init();
