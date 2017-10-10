#include "general/ArgParser.h"

#include "polPPD.C"

#include <string>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto sigfile = parser.getOptionVal<std::string>("--sigfile");
  const auto reffile = parser.getOptionVal<std::string>("--reffile");
  const auto sigtree = parser.getOptionVal<std::string>("--sigtree");
  const auto reftree = parser.getOptionVal<std::string>("--reftree");
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "polPPD.root");
  const auto reflth = parser.getOptionVal<double>("--reflth", 0.0);

  polPPD(sigfile, reffile, sigtree, reftree, outfile, reflth);

  return 0;
}

#endif
