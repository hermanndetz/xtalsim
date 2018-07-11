#ifndef EASYLOGGINGPPINTERNAL_H
#define EASYLOGGINGPPINTERNAL_H

#include "easylogging++.h"

namespace el {

namespace base {
/// @brief Namespace containing constants used internally.
namespace consts {
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

static const base::type::char_t* kInfoLevelLogValue     =   ELPP_LITERAL("INFO");
static const base::type::char_t* kDebugLevelLogValue    =   ELPP_LITERAL("DEBUG");
static const base::type::char_t* kWarningLevelLogValue  =   ELPP_LITERAL("WARNING");
static const base::type::char_t* kErrorLevelLogValue    =   ELPP_LITERAL("ERROR");
static const base::type::char_t* kFatalLevelLogValue    =   ELPP_LITERAL("FATAL");
static const base::type::char_t* kVerboseLevelLogValue  =
  ELPP_LITERAL("VERBOSE"); // will become VERBOSE-x where x = verbose level
static const base::type::char_t* kTraceLevelLogValue    =   ELPP_LITERAL("TRACE");
static const base::type::char_t* kInfoLevelShortLogValue     =   ELPP_LITERAL("I");
static const base::type::char_t* kDebugLevelShortLogValue    =   ELPP_LITERAL("D");
static const base::type::char_t* kWarningLevelShortLogValue  =   ELPP_LITERAL("W");
static const base::type::char_t* kErrorLevelShortLogValue    =   ELPP_LITERAL("E");
static const base::type::char_t* kFatalLevelShortLogValue    =   ELPP_LITERAL("F");
static const base::type::char_t* kVerboseLevelShortLogValue  =   ELPP_LITERAL("V");
static const base::type::char_t* kTraceLevelShortLogValue    =   ELPP_LITERAL("T");
// Format specifiers - These are used to define log format
static const base::type::char_t* kAppNameFormatSpecifier          =      ELPP_LITERAL("%app");
static const base::type::char_t* kLoggerIdFormatSpecifier         =      ELPP_LITERAL("%logger");
static const base::type::char_t* kThreadIdFormatSpecifier         =      ELPP_LITERAL("%thread");
static const base::type::char_t* kSeverityLevelFormatSpecifier    =      ELPP_LITERAL("%level");
static const base::type::char_t* kSeverityLevelShortFormatSpecifier    =      ELPP_LITERAL("%levshort");
static const base::type::char_t* kDateTimeFormatSpecifier         =      ELPP_LITERAL("%datetime");
static const base::type::char_t* kLogFileFormatSpecifier          =      ELPP_LITERAL("%file");
static const base::type::char_t* kLogFileBaseFormatSpecifier      =      ELPP_LITERAL("%fbase");
static const base::type::char_t* kLogLineFormatSpecifier          =      ELPP_LITERAL("%line");
static const base::type::char_t* kLogLocationFormatSpecifier      =      ELPP_LITERAL("%loc");
static const base::type::char_t* kLogFunctionFormatSpecifier      =      ELPP_LITERAL("%func");
static const base::type::char_t* kCurrentUserFormatSpecifier      =      ELPP_LITERAL("%user");
static const base::type::char_t* kCurrentHostFormatSpecifier      =      ELPP_LITERAL("%host");
static const base::type::char_t* kMessageFormatSpecifier          =      ELPP_LITERAL("%msg");
static const base::type::char_t* kVerboseLevelFormatSpecifier     =      ELPP_LITERAL("%vlevel");
static const char* kDateTimeFormatSpecifierForFilename            =      "%datetime";
// Date/time
static const char* kDays[7]                         =      { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday" };
static const char* kDaysAbbrev[7]                   =      { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
static const char* kMonths[12]                      =      { "January", "February", "March", "Apri", "May", "June", "July", "August",
                                                             "September", "October", "November", "December"
                                                           };
static const char* kMonthsAbbrev[12]                =      { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
static const char* kDefaultDateTimeFormat           =      "%Y-%M-%d %H:%m:%s,%g";
static const char* kDefaultDateTimeFormatInFilename =      "%Y-%M-%d_%H-%m";

static const char* kAm                              =      "AM";
static const char* kPm                              =      "PM";

#ifdef ELPP_DEFAULT_PERFORMANCE_LOGGER
static const char* kPerformanceLoggerId                    =      ELPP_DEFAULT_PERFORMANCE_LOGGER;
#else
static const char* kPerformanceLoggerId                    =      "performance";
#endif

static const char* kNullPointer                            =      "nullptr";

static const char* kUnknownUser                            =      "user";
static const char* kUnknownHost                            =      "unknown-host";
#if defined(ELPP_DEFAULT_LOG_FILE)
static const char* kDefaultLogFile                         =      ELPP_DEFAULT_LOG_FILE;
#else
#  if ELPP_OS_UNIX
#      if ELPP_OS_ANDROID
static const char* kDefaultLogFile                         =      "logs/myeasylog.log";
#      else
static const char* kDefaultLogFile                         =      "logs/myeasylog.log";
#      endif  // ELPP_OS_ANDROID
#  elif ELPP_OS_WINDOWS
static const char* kDefaultLogFile                         =      "logs\\myeasylog.log";
#  endif  // ELPP_OS_UNIX
#endif  // defined(ELPP_DEFAULT_LOG_FILE)
#if !defined(ELPP_DISABLE_LOG_FILE_FROM_ARG)
static const char* kDefaultLogFileParam                    =      "--default-log-file";
#endif  // !defined(ELPP_DISABLE_LOG_FILE_FROM_ARG)
#if defined(ELPP_LOGGING_FLAGS_FROM_ARG)
static const char* kLoggingFlagsParam                      =      "--logging-flags";
#endif  // defined(ELPP_LOGGING_FLAGS_FROM_ARG)

static const char* kValidLoggerIdSymbols                   =
  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._ ";
static const char* kConfigurationComment                   =      "##";
static const char* kConfigurationLevel                     =      "*";
static const char* kConfigurationLoggerId                  =      "--";

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
}  // namespace consts
}  // namespace base
}  // namespace el

#endif // EASYLOGGINGPPINTERNAL_H
