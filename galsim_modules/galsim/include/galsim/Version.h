
// This file is auto-generated by SCons.  Do not edit.
#define GALSIM_MAJOR 2
#define GALSIM_MINOR 2
#define GALSIM_REVISION 4

#include <string>
#include <sstream>

namespace galsim {
    // Compiled versions of the above #define values.
    extern int major_version();
    extern int minor_version();
    extern int revision();

    // Returns string of the form "1.4.2"
    extern std::string version();

    // Checks if the compiled library version matches the #define values in this header file.
    inline bool check_version() {
        // Same code as version(), but inline, so we get the above values to compare
        // to the values compiled into the library.
        std::ostringstream oss;
        oss << GALSIM_MAJOR << '.' << GALSIM_MINOR << '.' << GALSIM_REVISION;
        return oss.str() == version();
    }
}