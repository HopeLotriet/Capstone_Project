const ParseBoy = require("./ParseBoy");
const Extract = require("./libs/processing");
const objParseBoy = new ParseBoy();

const onFileReady = (preppedFile) => objParseBoy.parseFile(preppedFile);

const parser = async (file) => {
  const processing = new Extract(file);
  const data = await processing.extractTextFile();
  return onFileReady(data);
};

module.exports = parser;
