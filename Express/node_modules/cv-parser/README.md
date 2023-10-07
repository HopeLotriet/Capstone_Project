# Resume Parser

Extracts resume information

## Installation
```
  npm install --save doc-resume-parser
```

## Usage

```
const parser = require("doc-resume-parser");
require("dotenv").config();

const file = process.env.FILE;

parser(file)
  .then((data) => console.log(data))
  .catch((error) => {
    console.log(error);
  });
```

