# https://aryoda.github.io/tutorials/tryCatchLog/tryCatchLog-intro-slides.html#/easiest-way-to-add-logging-of-all-conditions
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking
# options("tryCatchLog.write.error.dump.file" = TRUE) # dump for post-mortem analysis

flog.appender(appender.file("my_app.log"))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)