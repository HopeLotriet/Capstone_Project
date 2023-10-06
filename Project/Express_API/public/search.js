function searchApp() {
    return {
        searchQuery: "",
        searchResults: [],
        resultsPerPage: 3,
        currentPage: 1,
        totalPages: 0,

        async search() {
            this.currentPage = 1;
            await this.updateResults();
        },

        async updateResults() {
            try {
                const response = await fetch(`https://api.example.com/search?query=${this.searchQuery}`);
                if (!response.ok) {
                    throw new Error("Network response was not ok");
                }
                const data = await response.json();
                this.totalPages = Math.ceil(data.length / this.resultsPerPage);
                this.searchResults = data.slice(
                    (this.currentPage - 1) * this.resultsPerPage,
                    this.currentPage * this.resultsPerPage
                );
            } catch (error) {
                console.error("Error fetching data:", error);
            }
        },

        changePage(page) {
            this.currentPage = page;
            this.updateResults();
        },
    };
}
