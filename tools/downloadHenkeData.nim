import std / [httpclient, os]
import shell

const basePath = "https://henke.lbl.gov/optical_constants/sf.tar.gz"
const outPath = currentSourcePath().parentDir.parentDir / "resources/henke_form_factors"

let client = newHttpClient()
createDir(outPath)
writeFile(outPath / "sf.tar.gz", client.getContent(basePath))

shell:
  one:
    cd ($outPath)
    tar xzf sf.tar.gz
