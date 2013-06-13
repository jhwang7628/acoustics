#ifndef PRINT_MSG_H
#   define PRINT_MSG_H

#define TERM_COLOR_BEGIN     "\033["
#define TERM_COLOR_END       "\033[0m"

#define TERM_COLOR_BOLD      1
#define TERM_COLOR_BLINK     5

#define TERM_COLOR_RED       31
#define TERM_COLOR_GREEN     32
#define TERM_COLOR_YELLOW    33
#define TERM_COLOR_BLUE      34
#define TERM_COLOR_MAGENTA   35
#define TERM_COLOR_WHITE     37

#define TERM_COLOR_BLACK_BG  40
#define TERM_COLOR_RED_BG    41
#define TERM_COLOR_BLUE_BG   44

#define PRINT_WARNING(...)                                  \
        {                                                   \
            fprintf(stderr, "%s%d;%dmWARNING:%s %s%dm",     \
                    TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                    TERM_COLOR_YELLOW, TERM_COLOR_END,      \
                    TERM_COLOR_BEGIN, TERM_COLOR_WHITE);    \
            fprintf(stderr, __VA_ARGS__);                   \
            fprintf(stderr, TERM_COLOR_END);                \
        }

#define PRINT_ERROR(...)                                    \
        {                                                   \
            fprintf(stderr, "%s%d;%dmERROR:%s %s%dm",       \
                    TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                    TERM_COLOR_RED, TERM_COLOR_END,         \
                    TERM_COLOR_BEGIN, TERM_COLOR_YELLOW);   \
            fprintf(stderr, __VA_ARGS__);                   \
            fprintf(stderr, TERM_COLOR_END);                \
        }

#define PRINT_MSG(...)                                     \
        {                                                  \
            printf("%s%d;%dmMESSAGE:%s %s%dm",             \
                   TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                   TERM_COLOR_GREEN, TERM_COLOR_END,       \
                   TERM_COLOR_BEGIN, TERM_COLOR_WHITE);    \
            printf(__VA_ARGS__);                           \
            printf(TERM_COLOR_END);                        \
        }

#endif
