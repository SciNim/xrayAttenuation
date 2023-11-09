import std / [strutils, httpclient, asyncdispatch, sequtils, parseutils, os]

#[
Downloads all the NIST data files for the atomic scattering factors.
]#

const basePath = "https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=$#&Formula=&gtype=4&range=U&lower=200&upper=700&density="
const outPath = currentSourcePath().parentDir.parentDir / "resources/nist_form_factors"

proc downloadFile(element: int): Future[(int, string)] {.async.} =
  var client = newAsyncHttpClient()
  let data = await client.getContent(basePath % [element.intToStr(2)])
  return (element, data)

proc extractBetween(s: string, start, stop: string): string =
  let idxStart = s.find(start)
  let idxStop = s.find(stop, start = idxStart + start.len)
  echo "from : ", (idxStart, idxStop)
  result = s[idxStart + start.len ..< idxStop]

proc cleanElementName(s: string): string =
  result = s.multiReplace(
    ("<BR>", " "), ("<i>", ""), ("</i>", ""), ("<I>", ""), ("</I>", ""), ("\r", ""), ("\n", ""),
    ("=", "")
  ).replace("  ", " ")

proc cleanData(s: string): seq[string] =
  result = s.strip.splitLines.mapIt(it.strip.split(" ").filterIt(it.len > 0).join(","))

proc afterDash(s: string): string =
  const header = "___"
  var buf = ""
  let start = s.find(header)
  let idx = parseWhile(s, buf, validChars = {'_'}, start = start)
  result = s[start + idx .. ^1]

proc extractData(html: string): tuple[element: string, data: seq[string]] =
  let elementNode = html.extractBetween("<b>", ")")
  let element = elementNode.extractBetween(" ", "&#160;").strip
  # `nm` from the last column, `λ`
  let data = html.extractBetween("nm", "</BODY>")
  result[0] = element #html.extractBetween("<B>", "</B>").cleanElementName()
  result[1] = data.cleanData()

proc writeData(tup: tuple[element: string, data: seq[string]]) =
  createDir(outpath)
  var elementName = tup.element
  # special case for "Carbon" as it also contains "Graphite" in the name
  if "Carbon" in elementName:
    doAssert "," in elementName
    elementName = "data_element_Carbon_Z_6.csv"
  echo "Data: ", tup.data
  let columns = ["E [keV]", "f1 [e atom⁻¹]", "f2 [e atom⁻¹]", "μ/ρ [cm² g⁻¹]", "σ/ρ [cm² g⁻¹]", "μ/ρ total [cm² g⁻¹]", "μ/ρ_K [cm² g⁻¹]", "λ [nm]"]
  var outfile = open(outpath / "data_element_" & elementName & ".csv", fmWrite)
  outfile.write(columns.join(",") & "\n")
  for line in tup.data:
    let lineSpl = line.split(",")
    doAssert lineSpl.len == 8
    outfile.write(line & "\n")
  outfile.close()

proc main(download = false, extract = false) =

  if download:
    # 1. iterate all elements and download files
    var futs = newSeq[Future[(int, string)]]()
    for element in 1 ..< 92:
      futs.add downloadFile(element)

    echo "INFO: Downloading all files..."
    while futs.anyIt(not it.finished()):
      poll()

    echo "INFO: Downloading done! Writing to ", outpath
    var files = newSeq[(int, string)]()
    for fut in futs:
      files.add waitFor(fut)

    for f in files:
      let (element, data) = (f[0], f[1])
      writeFile(outpath / "nist_element_" & $element & ".txt", data)

  if extract:
    for f in walkFiles(outpath / "*.txt"):
      echo "looking at: ", f
      f.readFile.extractData.writeData()

when isMainModule:
  import cligen
  dispatch main
