# From https://stackoverflow.com/questions/59537482/how-to-get-line-number-of-a-function-call-in-r
caller_info <- function(fmtstring = NULL, level = 1) # https://stackoverflow.com/q/59537482/684229
{
    x <- .traceback(x = level + 1)

    i <- 1
    repeat { # loop for subexpressions case; find the first one with source reference
        srcref <- getSrcref(x[[i]])
        if (is.null(srcref)) {
            if (i < length(x)) {
                i <- i + 1
                next
            } else {
                warning("caller_info(): not found\n")
                return(NULL)
            }
        }
        srcloc <- list(fun = getSrcref(x[[i + 1]]), file = getSrcFilename(x[[i]]), line = getSrcLocation(x[[i]]))
        break
    }

    if (is.null(fmtstring)) {
        return(srcloc)
    }

    fmtstring <- sub("%f", paste0(srcloc$fun, collapse = ""), fmtstring)
    fmtstring <- sub("%s", srcloc$file, fmtstring)
    fmtstring <- sub("%l", srcloc$line, fmtstring)
    fmtstring
}

write_csv_append_md <- function(df_table, file, version_code = NULL, SAVE_RDS = FALSE, ...) {
    docker_container <- Sys.getenv("CONTAINER_VERSION")
    str_func_call <- caller_info("Function: %f Location: %s#%l")
    if (is.null(version_code)) {
        version_code <- as.numeric(Sys.time())
    }


    str_git_branch <- system("git branch --show-current", intern = TRUE)[[1]]
    str_git_commit <- system("git rev-parse HEAD", intern = TRUE)[[1]]
    str_git_remote <- system("git config --get remote.origin.url", intern = TRUE)[[1]]

    str_git_info <- sprintf("%s@%s", str_git_branch, str_git_commit)
    str_git_repository <- sprintf("%s", str_git_remote)

    df_table$PROVENANCE_VERSION_CODE <- version_code
    df_table$PROVENANCE_DOCKER_IMAGE <- docker_container
    df_table$PROVENANCE_FUNCTION_CALL <- str_func_call
    df_table$PROVENANCE_GIT_REPOSITORY <- str_git_repository
    df_table$PROVENANCE_GIT_COMMIT <- str_git_info
    if (SAVE_RDS) {
        saveRDS(df_table, file = file, ...)
    } else {
        write.csv(df_table, file = file, ...)
    }
    return(list(df_table = df_table, file = file, version_code = version_code, docker_container = docker_container, str_func_call = str_func_call, git_info = str_git_info, git_repo = str_git_repository))
}
