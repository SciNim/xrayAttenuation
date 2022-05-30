import std / [strutils, httpclient, asyncdispatch, sequtils, parseutils, os]

const basePath = "https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z$#.html"
const outPath = currentSourcePath().parentDir.parentDir / "resources/nist_mass_attenuation"

proc downloadFile(element: int): Future[string] {.async.} =
  var client = newAsyncHttpClient()
  return await client.getContent(basePath % [element.intToStr(2)])

proc extractBetween(s: string, start, stop: string): string =
  let idxStart = s.find(start)
  let idxStop = s.find(stop, start = idxStart + start.len)
  result = s[idxStart + start.len ..< idxStop]

proc cleanElementName(s: string): string =
  result = s.multiReplace(
    ("<BR>", " "), ("<i>", ""), ("</i>", ""), ("<I>", ""), ("</I>", ""), ("\r", ""), ("\n", ""),
    ("=", "")
  ).replace("  ", " ")

proc cleanData(s: string): seq[string] =
  result = s.strip.splitLines.mapIt(it.strip.split(" ").filterIt(it.len > 0).join("\t"))

proc afterDash(s: string): string =
  const header = "___"
  var buf = ""
  let start = s.find(header)
  let idx = parseWhile(s, buf, validChars = {'_'}, start = start)
  result = s[start + idx .. ^1]

proc extractData(html: string): tuple[element: string, data: seq[string]] =
  let preNode = html.extractBetween("<PRE>", "</PRE>")
  let data = preNode.afterDash().afterDash()
  result[0] = html.extractBetween("<B>", "</B>").cleanElementName()
  result[1] = data.cleanData()

proc writeData(tup: tuple[element: string, data: seq[string]]) =
  let header = "Lines\tEnergy[MeV]\tμ/ρ\tμ_en/ρ\n"
  createDir(outpath)
  var outfile = open(outpath / "data_element_" & tup.element.replace(" ", "_") & ".csv", fmWrite)
  outfile.write(header)
  for line in tup.data:
    let lineSpl = line.split("\t")
    if lineSpl.len == 3:
      outfile.write("-\t" & line & "\n") # prepend an empty tab if no fluorescence line
    else:
      doAssert lineSpl.len == 4 # includes fluorescence line
      outfile.write(line & "\n")
  outfile.close()

# 1. iterate all elements and download files
var futs = newSeq[Future[string]]()
for element in 1 ..< 92:
  futs.add downloadFile(element)

echo "INFO: Downloading all files..."
while futs.anyIt(not it.finished()):
  poll()

echo "INFO: Downloading done! Writing to ", outpath
var files = newSeq[string]()
for fut in futs:
  files.add waitFor(fut)

for f in files:
  f.extractData.writeData()
