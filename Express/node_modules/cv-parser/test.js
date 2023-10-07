const parser = require("./src");
require("dotenv").config();

const file = process.env.FILE;

parser(file)
  .then((data) => console.log(data))
  .catch((error) => {
    console.log(error);
  });
