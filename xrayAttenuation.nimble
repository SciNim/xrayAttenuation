# Package

version       = "0.4.3"
author        = "Vindaar"
description   = "Library for X-ray reflectivity and transmission / absorption through matter"
license       = "MIT"

# Dependencies

requires "nim >= 1.6.0"
requires "ggplotnim >= 0.5.1"
requires "unchained >= 0.1.9"
requires "numericalnim"
requires "cligen"

task ciDeps, "Install CI dependencies":
  exec "nimble install nimhdf5"

task test, "Run the tests":
  exec "nim c -r tests/tests.nim"
  exec "nim c -r playground/mean_free_path.nim"
  # only compile for now
  ## Does not compile right now due to some weird nimhdf5 regression!
  #exec "nim c playground/llnl_telescope_reflectivity.nim"
