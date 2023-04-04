
.onAttach = function(libname, pkgname)
{
    temp = packageDescription("sklarsomega")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version, " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "Copyright (c) 2018-2023 John Hughes\n", sep = "")
    msg = paste(msg, 'For citation information, type citation("sklarsomega").\n', sep = "")
    msg = paste(msg, 'Type help(package = sklarsomega) to get started.\n', sep = "")
    packageStartupMessage(msg)
}

