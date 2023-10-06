document.addEventListener('alpine:init', () => {
    Alpine.data('searchApp', () => {
        return {
            // Initialize ChemDoodle sketcher
            chemDoodleCanvas: null,
            selectedCompound: {},

            initChemDoodleCanvas() {
                this.chemDoodleCanvas = new ChemDoodle.SketcherCanvas('drawingCanvas', 400, 400);
                this.chemDoodleCanvas.setIsotopesVisible(true);

                // Add event listener for structure search
                document.getElementById('searchButton').addEventListener('click', () => {
                    const molecule = this.chemDoodleCanvas.getMoleculeJSON();
                    if (molecule) {
                        this.clearSearchResults(); // Clear previous results
                        this.searchChemicalStructure(molecule);
                    } else {
                        console.error('Invalid chemical structure.');
                    }
                });

                // Add event listener for clearing the canvas
                document.getElementById('clearButton').addEventListener('click', () => {
                    this.clearChemDoodleCanvas();
                    this.clearSearchResults();
                });
            },

            clearSearchResults() {
                const resultsDiv = document.getElementById('searchResults');
                resultsDiv.innerHTML = '';
            },

            clearChemDoodleCanvas() {
                if (this.chemDoodleCanvas) {
                    this.chemDoodleCanvas.clear();
                }
            },

            saveStructure() {
                const molfile = this.chemDoodleCanvas.getMoleculeMolfile();
                // Send the 'molfile' to your server for saving or download
            },

            undo() {
                if (this.chemDoodleCanvas) {
                    this.chemDoodleCanvas.undo();
                }
            },

            redo() {
                if (this.chemDoodleCanvas) {
                    this.chemDoodleCanvas.redo();
                }
            },

            zoomIn() {
                if (this.chemDoodleCanvas) {
                    this.chemDoodleCanvas.zoomIn();
                }
            },

            zoomOut() {
                if (this.chemDoodleCanvas) {
                    this.chemDoodleCanvas.zoomOut();
                }
            },

            exportAsImage() {
                if (this.chemDoodleCanvas) {
                    const imageData = this.chemDoodleCanvas.toDataURL('image/png');
                    // You can handle the image data (e.g., display, save, or download)
                }
            },

            displayCompoundDetails(compound) {
                this.selectedCompound = compound;
            },

            // Function to send a chemical structure to your search API
            searchChemicalStructure(molecule) {
                // Make an AJAX POST request to your search API endpoint
                fetch('/search-api', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({ molecule }),
                })
                .then(response => response.json())
                .then(data => {
                    // Handle the response data (e.g., display results in #searchResults)
                    this.displaySearchResults(data.results);
                })
                .catch(error => {
                    console.error('Error:', error);
                });
            },

            displaySearchResults(results) {
                const resultsDiv = document.getElementById('searchResults');
                const ul = document.createElement('ul');

                if (results.length === 0) {
                    const li = document.createElement('li');
                    li.textContent = 'No results found.';
                    ul.appendChild(li);
                } else {
                    results.forEach(result => {
                        const li = document.createElement('li');
                        const link = document.createElement('a');
                        link.textContent = result.Name;
                        link.href = '#'; // Add a link to display more details
                        link.addEventListener('click', () => {
                            this.displayCompoundDetails(result);
                        });
                        li.appendChild(link);
                        ul.appendChild(li);
                    });
                }

                resultsDiv.appendChild(ul);
            },

            // You can add more methods and properties here as needed
        };
    });
});
