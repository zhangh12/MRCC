
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage('MRCC ', packageVersion('MRCC'), '  For help type ?mrcc ')
    packageStartupMessage('The most frequently updated version can be downloaded from https://github.com/zhangh12/ARTP2')
  }
}
