document.addEventListener('alpine:init', () => {
    Alpine.data('predictApp', () => {
        return {
            compoundName: '',
            predictionResult: null,
            async predictProperties() {
                try {
                    const response = await fetch('http://localhost:5000/predict_lipinski', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json'
                        },
                        body: JSON.stringify({
                            Name: this.compoundName,
                            // Add other input fields here
                        })
                    });

                    const data = await response.json();
                    this.predictionResult = data; // Update the reactive variable
                } catch (error) {
                    console.error('Error:', error);
                    // Handle error
                }
            }
        };
    });
});
