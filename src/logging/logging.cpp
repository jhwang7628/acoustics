#include <ctime>
#include "logging.h"
#include "utils/macros.h"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

char BUF__[1024];
static const char TIME_FORMAT[] = "%Y-%m-%d %H:%M:%S";

////////////////////////////////////////////////////////////////////////
LoggingManager LoggingManager::instance_;

/*!
 * Set the current logging level. If an invalid parameter is specified,
 * the current logging level is unchanged, and \c ERROR_RETURN is returned.
 *
 * \param level The logging level
 * \return If set the level successfully, \c SUCC_RETURN is returned,
 *         otherwise, \c ERROR_RETURN is returned.
 */
int LoggingManager::set_logging_level(LOGGING_LEVEL level)
{
    if ( level < 0 || level > LOG_NONE )
        return ERROR_RETURN;
    log_level_ = level;
    return SUCC_RETURN;
}

/*!
 * After calling this method, the logging information is redirected to 
 * the specified log file.
 *
 * If you didn't call this method, than the log information is output
 * to console as default.
 *
 * If the given file name is \c NULL, then the logging information would
 * be redirected to standard output \c stdout.
 *
 * \param file The log file
 * \return If set the level successfully, \c SUCC_RETURN is returned,
 *         otherwise, \c ERROR_RETURN is returned.
 */
int LoggingManager::set_logfile(const char* file, const char * mode)
{ 
    if ( !file && out_ != stdout )
    {
        fclose(out_);
        out_ = stdout;
        return SUCC_RETURN;
    }

    FILE* fp = fopen(file, mode);
    if ( !fp ) return ERROR_RETURN;

    out_ = fp;
    return SUCC_RETURN;
}

LoggingManager::LoggingManager():log_level_(LOG_NONE), out_(stdout)
{ 
}

LoggingManager::~LoggingManager()
{
    if ( out_ != stdout )
        fclose(out_);
}

void LoggingManager::logging_debug(const char * str)
{    
    if ( log_level_ > LOG_DEBUG ) return;
    time_t t;
    time(&t);
    struct tm* ct = localtime(&t);
    
    strftime(TIME__, 64, TIME_FORMAT, ct);
    fprintf(out_, "%s [DEBUG]   %s\n", TIME__, str);
}

void LoggingManager::logging_info(const char * str)
{    
    if ( log_level_ > LOG_INFO ) return;
    time_t t;
    time(&t);
    struct tm* ct = localtime(&t);

    strftime(TIME__, 64, TIME_FORMAT, ct);
    fprintf(out_, "%s [INFO]    %s\n", TIME__, str);
}

void LoggingManager::logging_warning(const char * str)
{    
    if ( log_level_ > LOG_WARNING ) return;
    time_t t;
    time(&t);
    struct tm* ct = localtime(&t);

    strftime(TIME__, 64, TIME_FORMAT, ct);
    fprintf(out_, "%s [WARNING] %s\n", TIME__, str);
}

void LoggingManager::logging_error(const char * str)
{    
    if ( log_level_ > LOG_ERROR ) return;
    time_t t;
    time(&t);
    struct tm* ct = localtime(&t);

    strftime(TIME__, 64, TIME_FORMAT, ct);
    fprintf(out_, "%s [ERROR]   %s\n", TIME__, str);
}

#ifdef USE_NAMESPACE
} // end namespace
#endif
