// pagination.js
// const jsonData = {{ data | tojson | safe }};
// const jsonData = JSON.parse(document.getElementById('data').textContent);

// Retrieve JSON data from the 'data-data' attribute of the 'jsonData' div element
const jsonDataElement = document.getElementById('jsonData');
const jsonData = JSON.parse(jsonDataElement.getAttribute('data-data'));

function pagination() {
    return {
        currentPage: 1,
        itemsPerPage: 100,
        data: jsonData,

        get paginatedData() {
            const start = (this.currentPage - 1) * this.itemsPerPage;
            const end = start + this.itemsPerPage;
            return this.data.slice(start, end);
        },

        get totalPages() {
            return Math.ceil(this.data.length / this.itemsPerPage);
        },

        prevPage() {
            if (this.currentPage > 1) {
                this.currentPage--;
            }
        },

        nextPage() {
            if (this.currentPage < this.totalPages) {
                this.currentPage++;
            }
        }
    };
}
