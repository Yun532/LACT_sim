#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <sys/types.h>
#include <stdint.h>
#include <sys/stat.h>
#include <time.h>
#if 0
#include <shmem.h>
#endif

using std::cout;

int main()
{
   cout << "Variable type      Typeid    Sizeof\n";
   cout << "  bool             " << typeid(bool).name()              << "         " << sizeof(bool) << "\n";
   cout << "  char             " << typeid(char).name()              << "         " << sizeof(char) << "\n";
   cout << "  signed char      " << typeid(signed char).name()       << "         " << sizeof(signed char) << "\n";
   cout << "  unsigned char    " << typeid(unsigned char).name()     << "         " << sizeof(unsigned char) << "\n";
   cout << "  wchar_t          " << typeid(wchar_t).name()           << "         " << sizeof(wchar_t) << "\n";
//   cout << "  signed wchar_t   " << typeid(signed wchar_t).name()    << "         " << sizeof(signed wchar_t) << "\n";
//   cout << "  unsigned wchar_t " << typeid(unsigned wchar_t).name()  << "         " << sizeof(unsigned wchar_t) << "\n";
   cout << "  short            " << typeid(short).name()             << "         " << sizeof(short) << "\n";
   cout << "  signed short     " << typeid(signed short).name()      << "         " << sizeof(signed short) << "\n";
   cout << "  unsigned short   " << typeid(unsigned short).name()    << "         " << sizeof(unsigned short) << "\n";
   cout << "  int              " << typeid(int).name()               << "         " << sizeof(int) << "\n";
   cout << "  signed int       " << typeid(signed int).name()        << "         " << sizeof(signed int) << "\n";
   cout << "  unsigned int     " << typeid(unsigned int).name()      << "         " << sizeof(unsigned int) << "\n";
   cout << "  long             " << typeid(long).name()              << "         " << sizeof(long) << "\n";
   cout << "  signed long      " << typeid(signed long).name()       << "         " << sizeof(signed long) << "\n";
   cout << "  unsigned long    " << typeid(unsigned long).name()     << "         " << sizeof(unsigned long) << "\n";
   cout << "  ssize_t          " << typeid(ssize_t).name()           << "         " << sizeof(ssize_t) << "\n";
   cout << "  size_t           " << typeid(size_t).name()            << "         " << sizeof(size_t) << "\n";
   cout << "  off_t            " << typeid(off_t).name()             << "         " << sizeof(off_t) << "\n";
   cout << "  ptrdiff_t        " << typeid(ptrdiff_t).name()         << "         " << sizeof(ptrdiff_t) << "\n";
   cout << "  time_t           " << typeid(time_t).name()            << "         " << sizeof(time_t) << "\n";
   cout << "  quad_t           " << typeid(quad_t).name()            << "         " << sizeof(quad_t) << "\n";
   cout << "  loff_t           " << typeid(loff_t).name()            << "         " << sizeof(loff_t) << "\n";
   cout << "  ino_t            " << typeid(ino_t).name()             << "         " << sizeof(ino_t) << "\n";
   cout << "  dev_t            " << typeid(dev_t).name()             << "         " << sizeof(dev_t) << "\n";
   cout << "  gid_t            " << typeid(gid_t).name()             << "         " << sizeof(gid_t) << "\n";
   cout << "  mode_t           " << typeid(mode_t).name()            << "         " << sizeof(mode_t) << "\n";
   cout << "  uid_t            " << typeid(uid_t).name()             << "         " << sizeof(uid_t) << "\n";
   cout << "  int8_t           " << typeid(int8_t).name()            << "         " << sizeof(int8_t) << "\n";
   cout << "  uint8_t          " << typeid(uint8_t).name()           << "         " << sizeof(uint8_t) << "\n";
   cout << "  int16_t          " << typeid(int16_t).name()           << "         " << sizeof(int16_t) << "\n";
   cout << "  uint16_t         " << typeid(uint16_t).name()          << "         " << sizeof(uint16_t) << "\n";
   cout << "  int32_t          " << typeid(int32_t).name()           << "         " << sizeof(int32_t) << "\n";
   cout << "  uint32_t         " << typeid(uint32_t).name()          << "         " << sizeof(uint32_t) << "\n";
   cout << "  int64_t          " << typeid(int64_t).name()           << "         " << sizeof(int64_t) << "\n";
   cout << "  uint64_t         " << typeid(uint64_t).name()          << "         " << sizeof(uint64_t) << "\n";
   cout << "  intmax_t         " << typeid(intmax_t).name()          << "         " << sizeof(intmax_t) << "\n";
   cout << "  uintmax_t        " << typeid(uintmax_t).name()         << "         " << sizeof(uintmax_t) << "\n";
   cout << "  long long        " << typeid(long long).name()         << "         " << sizeof(long long) << "\n";
   cout << "  uns. long long   " << typeid(unsigned long long).name()<< "         " << sizeof(unsigned long long) << "\n";
   cout << "  float            " << typeid(float).name()             << "         " << sizeof(float) << "\n";
   cout << "  double           " << typeid(double).name()            << "         " << sizeof(double) << "\n";
   cout << "  long double      " << typeid(long double).name()       << "         " << sizeof(long double) << "\n";
   cout << "  void *           " << typeid(void *).name()            <<  "        " << sizeof(void *) << "\n";

   return 0;
}
