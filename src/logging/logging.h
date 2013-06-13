#ifndef CARBINE_LOGGER_H
#   define CARBINE_LOGGER_H

#include <cstdio>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

//! Level of logging
enum LOGGING_LEVEL
{    
    //! Log info and debug message
    LOG_DEBUG = 0,
    //! Log only the info message
    LOG_INFO,
    //! Log info, debug, and warning message
    LOG_WARNING,
    //! Log all message
    LOG_ERROR,
    //! Do not log anything
    LOG_NONE
};

//! A simplified logger
/*!
 * I'm trying to implement it more like the \e Python
 * logging library since I think it's more convenient than
 * the ACE logging does.
 */
class LoggingManager
{
    public:        
        static LoggingManager instance_;

    private:
        LOGGING_LEVEL   log_level_;
        FILE*           out_;

        char            TIME__[64];   // time string buffer

    private:
        LoggingManager();

    public:
        // Destructor
        virtual ~LoggingManager();

        //! Retrieve the singleton instance.
        static LoggingManager& instance() { return instance_; }

        //! set the logging level
        int set_logging_level(LOGGING_LEVEL);
        //! set the log file that would be wrote
        int set_logfile(const char*, const char *); 
        // logging information
        void logging_debug  (const char *);
        void logging_info   (const char *);
        void logging_warning(const char *);
        void logging_error  (const char *);
        // ! current logging level
        LOGGING_LEVEL logging_level() { return log_level_; }        
};

extern char BUF__[];

#define LOGGING_DEBUG(...)       \
    if ( LoggingManager::instance().logging_level() <= LOG_DEBUG ) \
    { \
        sprintf(BUF__, __VA_ARGS__);                          \
        LoggingManager::instance().logging_debug(BUF__);      \
    }

#define LOGGING_INFO(...)        \
    if ( LoggingManager::instance().logging_level() <= LOG_INFO )  \
    { \
        sprintf(BUF__, __VA_ARGS__);                          \
        LoggingManager::instance().logging_info(BUF__);       \
    }

#define LOGGING_WARNING(...)     \
    if ( LoggingManager::instance().logging_level() <= LOG_WARNING ) \
    { \
        sprintf(BUF__, __VA_ARGS__);                            \
        LoggingManager::instance().logging_warning(BUF__);      \
    }

#define LOGGING_ERROR(...)       \
    if ( LoggingManager::instance().logging_level() <= LOG_ERROR ) \
    { \
        sprintf(BUF__, __VA_ARGS__);                          \
        LoggingManager::instance().logging_error(BUF__);      \
    }

#ifdef USE_NAMESPACE
} // end namespace
#endif

#endif
