const path = require("path");
const textract = require("textract");
const mime = require("mime-types");
const promisify = require("../");

class Extract {
  constructor(file) {
    this.path = file;
    this.mime = mime.lookup(file);
    this.ext = mime.extension(this.mime);
    this.raw = "";
    this.name = path.basename(file);
  }

  extractTextFile = async () => {
    const data = await promisify(textract.fromFileWithPath, this.path, {
      preserveLineBreaks: true,
    });
    const cleaned = this.cleanTextByRows(data);
    this.raw = cleaned;
    return this;
  };

  cleanTextByRows = (data) => {
    const result =
      data
        .split("\n")
        .map((row) => row.replace(/\r?\n|\r|\t|\n/g, "").trim())
        .join("\n") + "\n{end}";

    return result;
  };
}

module.exports = Extract;
