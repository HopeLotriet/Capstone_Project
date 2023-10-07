const {
  parseDictionaryRegular,
  parseDictionaryProfiles,
  parseDictionaryTitles,
  parseDictionaryInline,
} = require("./libs/parser");

const Resume = require("./Resume");

/**
 *
 * @constructor
 */
function ParseBoy() {}

/**
 *
 * @param PreparedFile
 * @param cbGetResume
 */
ParseBoy.prototype.parseFile = function(preppedFile) {
  const rawFileData = preppedFile.raw;
  const resume = new Resume();
  let rows = rawFileData.split("\n");

  parseDictionaryRegular(rawFileData, resume);
  rows.forEach((row, i) => {
    row = parseDictionaryProfiles(row, resume);
    parseDictionaryTitles(resume, rows, i);
    parseDictionaryInline(resume, row);
  });
  return resume;
};

module.exports = ParseBoy;
