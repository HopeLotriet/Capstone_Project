function search() {
    return {
        searchQuery: "",
        predictionResults: null,
        showError: false,

        async search() {
            try {
                const response = await fetchPredictionAPI(this.searchQuery);
                if (!response.ok) {
                    throw new Error("Network response was not ok");
                }
                const data = await response.json();
                this.predictionResults = data;
                this.showError = false; // Reset error state
            } catch (error) {
                console.error("Error fetching data:", error);
                this.predictionResults = null;
                this.showError = true;
            }
        },
    };
}

async function fetchPredictionAPI(query) {
    // Replace with your actual API endpoint
    const apiUrl = "https://your-chemical-prediction-api.com/predict";
    const requestBody = { compound: query };

    return await fetch(apiUrl, {
        method: "POST",
        headers: {
            "Content-Type": "application/json",
        },
        body: JSON.stringify(requestBody),
    });
}
