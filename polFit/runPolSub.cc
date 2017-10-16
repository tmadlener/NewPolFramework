#include "general/ArgParser.h"

#include "polSub.C"

#include <string>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto genFile = parser.getOptionVal<std::string>("--file");

  polSub(genFile);

  return 0;
}

#endif
