import std / [macros, macrocache]

## These tables are filled in the main program where the types are generated
const ElementTable = CacheTable"Elements"
const ElementSymbolTable = CacheTable"ElementSymbols"
const ElementChemToNameTable = CacheTable"ElementChemToName"
const ElementSeq = CacheSeq"ElementSeq"

proc getTypedescField(n: NimNode): NimNode =
  case n.typeKind
  of ntyTypeDesc:
    result = n[1]
  else:
    result = n

proc getName*(t: NimNode): NimNode =
  let typInst = t.getTypeInst
  let inner = getTypedescField(typInst)
  case inner.kind
  of nnkBracketExpr:
    # look up based on int
    expectKind(inner[1], nnkIntLit)
    result = ElementSeq[inner[1].intVal.int - 1]
  of nnkSym:
    result = inner
  else:
    error("Invalid AST found in `lookupInverseName`: " & $inner.treerepr)

macro lookupInverseName*(t: typed): untyped =
  result = getName(t)

proc lookupNameFromChemSymbol*(t: string): NimNode =
  result = ElementChemToNameTable[t]

macro lookupChemSymbol*(t: typed): untyped =
  let name = getName(t)
  result = newLit(ElementSymbolTable[name.strVal].strVal)

proc genTypeClass*(e: var seq[NimNode]): NimNode =
  ## Helper to generate a "type class" (using `|`) of multiple
  ## types, because for _reasons_ the Nim AST for that is nested infix
  ## calls apparently.
  if e.len == 2:
    result = nnkInfix.newTree(ident"|", e[0], e[1])
  else:
    let el = e.pop
    result = nnkInfix.newTree(ident"|",
                              genTypeClass(e),
                              el)

proc getArgStr*(n: NimNode): string =
  case n.kind
  of nnkIdent, nnkSym,
     nnkStrLit, nnkTripleStrLit, nnkRStrLit: result = n.strVal
  of nnkOpenSymChoice, nnkClosedSymChoice: result = n[0].strVal
  else:
    error("Invalid node for argument `" & $(n.repr) & " in aes macro!")
